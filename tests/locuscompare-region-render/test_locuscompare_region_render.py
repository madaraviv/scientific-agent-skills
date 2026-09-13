"""CLI contract, pre-fetched-input loaders, the offline synthetic demo end to end,
and the sibling-skill resolution the K-Dense layout relies on.

The composer (`locuscompare_region_render.py`) imports three sibling skills from
`skills/<sibling>/scripts/`. Tests that need it skip, with the composer's own
ImportError as the reason, when the siblings are neither installed in this
checkout nor pointed at by LOCUSCOMPARE_SIBLING_SKILLS_ROOT (see conftest.py).
The resolution tests themselves run everywhere: they build a throwaway skills/
tree in tmp_path so the absent-sibling path is exercised deterministically.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

import skill_contract  # noqa: E402

import cli as locuscompare  # noqa: E402

# Every argparse script under scripts/ answers --help in a subprocess (cli.py).
CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)

try:
    import locuscompare_region_render as composer  # noqa: E402
    SIBLING_ERROR = None
except ImportError as _e:  # a sibling skill is absent
    composer = None
    SIBLING_ERROR = str(_e)

needs_siblings = pytest.mark.skipif(SIBLING_ERROR is not None, reason=SIBLING_ERROR or "")

EXAMPLES = SKILL_ROOT / "scripts" / "examples"


# ----------------------- argparse / config contract (no siblings needed) -----------------------


def test_cli_requires_input_or_demo(tmp_path: Path):
    with pytest.raises(SystemExit) as excinfo:
        locuscompare.main(["--output", str(tmp_path)])
    assert excinfo.value.code == 2


@needs_siblings  # cli._run imports the composer before validating the config
def test_cli_creates_output_directory_on_invalid_config(tmp_path: Path):
    """--output is created up-front, even when the config later fails validation."""
    out = tmp_path / "sub1" / "sub2"
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps({"schema_version": "1.0"}))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(out)])
    assert out.is_dir()
    assert rc == 2


@needs_siblings  # cli._run imports the composer before validating the config
def test_cli_rejects_missing_lead(tmp_path: Path, capsys):
    config = {
        "exposure": {"trait_label": "x", "fetch": {"source": "eqtl_catalogue", "dataset_id": "X"}},
        "outcome": {"trait_label": "y", "fetch": {"source": "gwas_catalog", "accession": "Y"}},
    }
    config_path = tmp_path / "c.json"
    config_path.write_text(json.dumps(config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(tmp_path / "out")])
    assert rc == 2
    assert "lead" in capsys.readouterr().err.lower()


@needs_siblings  # cli._run imports the composer before validating the config
def test_cli_rejects_unsupported_exposure_source(tmp_path: Path, capsys):
    config = {
        "lead": {"variant_id": "1_1_A_T", "chromosome": "1", "position_bp": 1, "window_bp": 1000000},
        "exposure": {"trait_label": "x", "fetch": {"source": "made_up_source", "dataset_id": "X"}},
        "outcome": {"trait_label": "y", "fetch": {"source": "gwas_catalog", "accession": "Y"}},
    }
    config_path = tmp_path / "c.json"
    config_path.write_text(json.dumps(config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(tmp_path / "out")])
    assert rc == 2
    assert "made_up_source" in capsys.readouterr().err


@needs_siblings  # cli._run imports the composer before validating the config
def test_cli_rejects_missing_exposure_input_vector(tmp_path: Path, capsys):
    """Exposure must declare either `fetch:` or `sumstats_path:`; missing both is an error."""
    config = {
        "lead": {"variant_id": "1_1_A_T", "chromosome": "1", "position_bp": 1, "window_bp": 1000000},
        "exposure": {"trait_label": "x"},
        "outcome": {"trait_label": "y", "fetch": {"source": "gwas_catalog", "accession": "Y"}},
    }
    config_path = tmp_path / "c.json"
    config_path.write_text(json.dumps(config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(tmp_path / "out")])
    assert rc == 2
    err = capsys.readouterr().err
    assert "sumstats_path" in err and "fetch" in err


def test_load_config_yaml(tmp_path: Path):
    yaml_path = tmp_path / "c.yaml"
    yaml_path.write_text("schema_version: '1.0'\nfoo: bar\n")
    assert locuscompare._load_config(yaml_path) == {"schema_version": "1.0", "foo": "bar"}


def test_load_config_json(tmp_path: Path):
    json_path = tmp_path / "c.json"
    json_path.write_text('{"schema_version": "1.0", "foo": "bar"}')
    assert locuscompare._load_config(json_path) == {"schema_version": "1.0", "foo": "bar"}


def test_load_config_unknown_extension(tmp_path: Path):
    bad = tmp_path / "c.toml"
    bad.write_text("[section]\nfoo = 'bar'\n")
    with pytest.raises(ValueError, match="unsupported config extension"):
        locuscompare._load_config(bad)


def test_format_provenance_prefix_empty():
    assert locuscompare._format_provenance_prefix({}) == ""


def test_format_provenance_prefix_ot_release():
    assert locuscompare._format_provenance_prefix({"provenance": {"ot_release": "26.03"}}) == "OT release: 26.03 | "


def test_format_provenance_prefix_gwas_lookup():
    assert locuscompare._format_provenance_prefix(
        {"provenance": {"gwas_lookup_run_dir": "runs/gl/"}}
    ) == "gwas-lookup chain: runs/gl/ | "


# ----------------------- sibling-skill resolution -----------------------


def _stage_skill_copy(root: Path) -> Path:
    """Copy this skill into `root/skills/locuscompare-region-render/` so that the
    composer's `parents[2]` is a skills/ directory this test controls."""
    dst = root / "skills" / "locuscompare-region-render"
    shutil.copytree(SKILL_ROOT / "scripts", dst / "scripts")
    return dst


def _import_composer_in_subprocess(skill_dir: Path, *, python_path: str = "") -> subprocess.CompletedProcess:
    code = (
        "import sys; sys.path.insert(0, sys.argv[1]); "
        "import locuscompare_region_render as m; "
        "print('UKB_PPP_AVAILABLE=%s' % m.UKB_PPP_AVAILABLE)"
    )
    env = {"PYTHONDONTWRITEBYTECODE": "1", "PATH": "/usr/bin:/bin"}
    if python_path:
        env["PYTHONPATH"] = python_path
    return subprocess.run(
        [sys.executable, "-c", code, str(skill_dir / "scripts")],
        capture_output=True, text=True, timeout=120, env=env,
    )


def test_missing_required_sibling_raises_an_import_error_naming_the_skill(tmp_path: Path):
    """No sibling installed beside a copy of the skill: the failure is an ImportError
    whose message names the sibling skill and the path it was expected at, and it
    is not a bare ModuleNotFoundError traceback."""
    skill_dir = _stage_skill_copy(tmp_path)
    result = _import_composer_in_subprocess(skill_dir)
    assert result.returncode != 0
    last = result.stderr.strip().splitlines()[-1]
    assert last.startswith("ImportError: locuscompare-region-render needs the sibling skill")
    assert "'eqtl-catalogue-region-fetch'" in last
    assert str(tmp_path / "skills" / "eqtl-catalogue-region-fetch" / "scripts") in last
    assert "ModuleNotFoundError" not in result.stderr


@needs_siblings
def test_siblings_installed_beside_the_skill_are_found_without_pythonpath(tmp_path: Path):
    """The K-Dense layout: skills side by side under skills/, code under scripts/.
    Stage the three required siblings next to a copy of this skill and import the
    composer with an empty PYTHONPATH; resolution must come from the layout alone."""
    skill_dir = _stage_skill_copy(tmp_path)
    for module_name, skill_name in composer._REQUIRED_SIBLING_SKILLS.items():
        source = Path(sys.modules[module_name].__file__).resolve().parent
        target = tmp_path / "skills" / skill_name / "scripts"
        if not target.exists():
            shutil.copytree(source, target)
    result = _import_composer_in_subprocess(skill_dir)
    assert result.returncode == 0, result.stderr
    assert "UKB_PPP_AVAILABLE=" in result.stdout


@needs_siblings
def test_composer_loads_without_ukb_ppp_and_says_so():
    """ukb-ppp-region-fetch is optional: absent, the module still imports and the
    pQTL names resolve to None / a stand-in exception; present, they are real."""
    if composer.UKB_PPP_AVAILABLE:
        assert composer.UKBPPPClient is not None
        pytest.skip("ukb-ppp-region-fetch is installed; the absent path is exercised elsewhere")
    assert composer.UKBPPPClient is None
    assert composer.UKBPPPRegionResult is None
    assert issubclass(composer.UKBPPPAccessError, Exception)
    assert composer.dispatch_exposure_kind("UKB_PPP_EUR_SORT1") == composer.EXPOSURE_KIND_UKB_PPP


@needs_siblings
def test_pqtl_row_without_a_ukb_ppp_client_is_a_named_refusal(tmp_path: Path):
    """A `UKB_PPP_*` study id with no pQTL client raises Tier2NotAvailable (the
    documented fallback signal), not a TypeError from calling None."""
    mapping = composer.StudyIdMapping(
        ot_left_study_id="UKB_PPP_EUR_SORT1",
        ot_right_study_id="GCST90269602",
        gwas_catalog_accession="GCST90269602",
        ukb_ppp_protein_label="SORT1",
        ukb_ppp_ancestry="EUR",
    )
    with pytest.raises(composer.Tier2NotAvailable, match="ukb_ppp_client"):
        composer.render_tier2_for_lead(
            lead_variant_id="1_109274968_G_T", chromosome="1",
            lead_position_bp=109_274_968, window_bp=1_000_000,
            study_mapping=mapping, eqtl_client=object(), gwas_client=object(),
            ld_client=None, out_path=tmp_path / "x.png", ot_release="26.03",
            ukb_ppp_client=None,
        )


@needs_siblings
def test_cli_refuses_ukb_ppp_source_when_the_sibling_is_absent(tmp_path: Path, capsys):
    if composer.UKB_PPP_AVAILABLE:
        pytest.skip("ukb-ppp-region-fetch is installed; source=ukb_ppp is accepted here")
    config = {
        "lead": {"variant_id": "1_109274968_G_T", "chromosome": "1",
                 "position_bp": 109274968, "window_bp": 1000000},
        "exposure": {"trait_label": "SORT1 plasma sortilin",
                     "fetch": {"source": "ukb_ppp", "protein_label": "SORT1", "ancestry": "EUR"}},
        "outcome": {"trait_label": "y", "fetch": {"source": "gwas_catalog", "accession": "GCST90269602"}},
    }
    config_path = tmp_path / "c.json"
    config_path.write_text(json.dumps(config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(tmp_path / "out")])
    assert rc == 2
    err = capsys.readouterr().err
    assert "ukb-ppp-region-fetch" in err and "ukb_ppp" in err


# ----------------------- pre-fetched TSV input mode -----------------------


CANONICAL_TSV_HEADER = "variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\tp"


def _write_canonical_tsv(path: Path, rows: list[str]) -> None:
    path.write_text(CANONICAL_TSV_HEADER + "\n" + "\n".join(rows) + "\n")


@needs_siblings
def test_prefetched_load_sumstats_tsv_parses_canonical_schema(tmp_path: Path):
    from _prefetched import load_sumstats_tsv

    tsv = tmp_path / "exposure.tsv"
    _write_canonical_tsv(tsv, [
        "1_500000_A_T\t1\t500000\tA\tT\t0.5\t0.05\t1e-22",
        "1_500100_C_G\t1\t500100\tC\tG\t0.3\t0.05\t1e-9",
    ])
    variants = load_sumstats_tsv(tsv)
    assert len(variants) == 2
    assert variants[0].variant_id == "1_500000_A_T"
    assert variants[0].ref == "A" and variants[0].alt == "T"
    assert variants[0].beta == 0.5 and variants[0].se == 0.05
    assert variants[0].p_value == 1e-22


@needs_siblings
def test_prefetched_load_sumstats_tsv_rejects_missing_columns(tmp_path: Path):
    from _prefetched import PrefetchedSchemaError, load_sumstats_tsv

    tsv = tmp_path / "bad.tsv"
    tsv.write_text("variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\n")
    with pytest.raises(PrefetchedSchemaError, match="missing required columns"):
        load_sumstats_tsv(tsv)


@needs_siblings
def test_prefetched_load_sumstats_tsv_drops_rows_with_na(tmp_path: Path):
    from _prefetched import load_sumstats_tsv

    tsv = tmp_path / "exposure.tsv"
    _write_canonical_tsv(tsv, [
        "1_500000_A_T\t1\t500000\tA\tT\t0.5\t0.05\t1e-22",
        "1_500100_C_G\t1\t500100\tC\tG\tNA\tNA\tNA",
    ])
    variants = load_sumstats_tsv(tsv)
    assert [v.variant_id for v in variants] == ["1_500000_A_T"]


@needs_siblings
def test_prefetched_load_sumstats_tsv_does_not_range_check(tmp_path: Path):
    """The loader parses and drops NA rows; it does not validate ranges. A p of 2.0
    and a negative se are loaded as given (documented in references/input_schema.md)."""
    from _prefetched import load_sumstats_tsv

    tsv = tmp_path / "loose.tsv"
    _write_canonical_tsv(tsv, ["1_500000_A_G\t1\t500000\tA\tG\t0.5\t-0.05\t2.0"])
    variants = load_sumstats_tsv(tsv)
    assert len(variants) == 1 and variants[0].se == -0.05 and variants[0].p_value == 2.0


@needs_siblings
def test_prefetched_load_synthetic_ld_parses_two_columns(tmp_path: Path):
    from _prefetched import load_synthetic_ld

    tsv = tmp_path / "ld.tsv"
    tsv.write_text("partner_variant_id\tr2\n1_500100_C_G\t0.85\n1_499900_G_A\t0.42\n")
    assert load_synthetic_ld(tsv) == {"1_500100_C_G": 0.85, "1_499900_G_A": 0.42}


@needs_siblings
def test_prefetched_load_synthetic_ld_rejects_out_of_range(tmp_path: Path):
    from _prefetched import PrefetchedSchemaError, load_synthetic_ld

    tsv = tmp_path / "ld.tsv"
    tsv.write_text("partner_variant_id\tr2\n1_500100_C_G\t1.5\n")
    with pytest.raises(PrefetchedSchemaError, match="out of"):
        load_synthetic_ld(tsv)


@needs_siblings
def test_prefetched_load_synthetic_gene_track_parses(tmp_path: Path):
    from _prefetched import load_synthetic_gene_track

    tsv = tmp_path / "genes.tsv"
    tsv.write_text(
        "gene_symbol\tstart\tend\tstrand\tbiotype\n"
        "DEMOGENE_A\t100\t200\t+\tprotein_coding\n"
        "DEMOGENE_B\t450\t550\t-\tprotein_coding\n"
    )
    genes = load_synthetic_gene_track(tsv)
    assert [g.gene_symbol for g in genes] == ["DEMOGENE_A", "DEMOGENE_B"]
    assert genes[0].start == 100 and genes[0].end == 200 and genes[0].strand == "+"
    assert genes[1].strand == "-"


# ----------------------- offline synthetic demo, end to end -----------------------


@needs_siblings
def test_prefetched_synthetic_demo_runs_end_to_end_offline(tmp_path: Path):
    """01_synthetic_demo: TSV sumstats + synthetic LD matrix + synthetic gene track,
    no network. A PNG, a manifest and a report land in --output. The lead
    1_500000_A_T is an A/T SNP, so it is the one palindromic-excluded pair."""
    out = tmp_path / "out"
    rc = locuscompare.main(["--input", str(EXAMPLES / "01_synthetic_demo" / "config.json"),
                            "--output", str(out)])
    assert rc == 0, "01_synthetic_demo failed offline"
    plot = out / "1_500000_A_T_full_locuscompare.png"
    assert plot.is_file() and plot.stat().st_size > 10_000
    manifest = yaml.safe_load((out / "manifest.yaml").read_text())
    assert manifest["lead_variant_id"] == "1_500000_A_T"
    assert manifest["n_pairs"] == 200
    assert manifest["n_palindromic_excluded"] == 1
    block = manifest["render_block"]
    assert block["ld_panel"] == "synthetic" and block["plink_version"] == "prefetched"
    assert block["scatter_downsampled"] is False
    assert (out / "report.md").is_file()


@needs_siblings
def test_demo_flag_runs_the_synthetic_demo_by_name(tmp_path: Path):
    out = tmp_path / "out"
    rc = locuscompare.main(["--demo", "01_synthetic_demo", "--output", str(out)])
    assert rc == 0
    assert (out / "1_500000_A_T_full_locuscompare.png").is_file()


@needs_siblings
def test_lead_rs_id_propagates_into_manifest_and_report(tmp_path: Path):
    """Optional `lead.rs_id` flows config -> spec -> manifest + report. The join is
    on variant_id; rs_id is human-readable metadata only."""
    out = tmp_path / "out"
    base_dir = EXAMPLES / "01_synthetic_demo"
    base_config = json.loads((base_dir / "config.json").read_text())
    base_config["lead"]["rs_id"] = "rs99999999"
    base_config["exposure"]["sumstats_path"] = str(base_dir / "exposure.tsv")
    base_config["outcome"]["sumstats_path"] = str(base_dir / "outcome.tsv")
    base_config["ld"]["ld_matrix_path"] = str(base_dir / "ld_matrix.tsv")
    base_config["gene_track"]["genes_path"] = str(base_dir / "genes.tsv")
    override_path = tmp_path / "config_with_rsid.json"
    override_path.write_text(json.dumps(base_config))
    rc = locuscompare.main(["--input", str(override_path), "--output", str(out)])
    assert rc == 0
    manifest = yaml.safe_load((out / "manifest.yaml").read_text())
    assert manifest["lead_rs_id"] == "rs99999999"
    assert "rs99999999" in (out / "report.md").read_text()


def test_focal_gene_highlight_styles_matching_label():
    """render_gene_track bolds and tints the label of the gene matching
    `focal_gene_symbol` and leaves the others in the default style."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from regional_plot import FOCAL_GENE_COLOR, GeneTrackEntry, render_gene_track

    genes = [
        GeneTrackEntry(gene_symbol="GENE_A", start=100, end=200, strand="+"),
        GeneTrackEntry(gene_symbol="SORT1", start=450, end=550, strand="-"),
        GeneTrackEntry(gene_symbol="GENE_C", start=800, end=900, strand="+"),
    ]
    fig, ax = plt.subplots()
    render_gene_track(ax, genes, xlim_bp=(0, 1000), focal_gene_symbol="SORT1")
    texts = {t.get_text(): t for t in ax.texts}
    sort1_label = next(t for k, t in texts.items() if "SORT1" in k)
    other_label = next(t for k, t in texts.items() if "GENE_A" in k)
    assert sort1_label.get_fontweight() == "bold"
    assert sort1_label.get_color() == FOCAL_GENE_COLOR
    assert other_label.get_fontweight() == "normal"
    assert other_label.get_color() != FOCAL_GENE_COLOR
    plt.close(fig)
