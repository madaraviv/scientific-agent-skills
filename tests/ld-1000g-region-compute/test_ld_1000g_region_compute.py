"""Unit tests for the on-demand 1000G LD client (plink 1.9 subprocess) and its CLI.

The region fetch (requests / pysam) and the plink subprocess are mocked, so these
tests are offline and need neither plink, pysam, nor the network. The live smoke
test (real plink 1.9 against a 1000G window fetched from EBI) is in
test_live_ld_1000g_region_compute.py.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# Inject the skill's scripts dir so the bare-name modules resolve under the
# kebab-cased skill directory (`skills/ld-1000g-region-compute/`), which is not
# a valid Python identifier.
SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "ld-1000g-region-compute"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

import skill_contract  # noqa: E402

import ld_1000g_region_compute as cli_module  # noqa: E402
import ondemand_client as oc_module  # noqa: E402
from ondemand_client import (  # noqa: E402
    OnDemand1000GLDClient,
    OnDemandLDError,
    OnDemandLDPair,
    OnDemandLDResult,
    _parse_ld,
)

# Every argparse script under scripts/ answers --help in a subprocess.
CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)


# -----------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------


def _patch_plink_version(monkeypatch, version: str = "PLINK v1.90b6.27 64-bit (2023-05-09)") -> None:
    monkeypatch.setattr(oc_module, "_detect_plink_version", lambda _bin: version)


def _patch_samples(monkeypatch, samples: list[str] | None = None) -> None:
    samples = samples or ["NA12878", "NA12879", "NA12891"]
    monkeypatch.setattr(
        oc_module, "_resolve_super_pop_samples",
        lambda _super_pop, _cache_dir: list(samples),
    )


def _patch_region_fetch(monkeypatch, tmp_path: Path) -> Path:
    """Replace the tabix fetch with a no-op that returns an empty .vcf.gz path."""
    fake_vcf = tmp_path / "fake_region.vcf.gz"
    fake_vcf.write_bytes(b"")
    monkeypatch.setattr(
        oc_module, "_fetch_region_vcf",
        lambda chrom, start, end, cache_dir: fake_vcf,
    )
    return fake_vcf


def _write_ld(out_prefix: Path, lead: str, partner_pairs: list[tuple[str, float]]) -> None:
    """Emit a plink 1.9 .ld file at <out_prefix>.ld.

    plink 1.9 emits whitespace-separated columns:
    `CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2`.
    """

    def _to_panel(ot_id: str) -> str:
        return ot_id.replace("_", ":", 3)

    rows = [" CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2"]
    for partner, r2 in partner_pairs:
        rows.append(
            f"   2     36910110  {_to_panel(lead)}    "
            f"   2     36932656  {_to_panel(partner)}    {r2:.6f}"
        )
    Path(f"{out_prefix}.ld").write_text("\n".join(rows) + "\n")


# -----------------------------------------------------------------
# Construction / validation
# -----------------------------------------------------------------


def test_client_validates_plink_present(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    c = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    assert c.plink_version.startswith("PLINK v1.90")
    assert c.super_pop == "EUR"


def test_client_raises_when_plink_missing(monkeypatch, tmp_path):
    def boom(_bin):
        raise OnDemandLDError("plink not found")
    monkeypatch.setattr(oc_module, "_detect_plink_version", boom)
    with pytest.raises(OnDemandLDError, match="plink not found"):
        OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)


def test_detect_plink_version_names_the_missing_binary(monkeypatch):
    """A missing binary raises with the install hint, before any subprocess runs."""
    monkeypatch.setattr(oc_module.shutil, "which", lambda _bin: None)
    with pytest.raises(OnDemandLDError, match="plink binary not found"):
        oc_module._detect_plink_version("plink-not-here")


# -----------------------------------------------------------------
# r2_with_lead happy path
# -----------------------------------------------------------------


def test_r2_with_lead_parses_ld_output(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    captured_cmd: dict = {}

    def fake_run(cmd, capture_output, text, check, timeout):
        captured_cmd["cmd"] = cmd
        out_prefix = Path(cmd[cmd.index("--out") + 1])
        _write_ld(out_prefix, lead="2_36910110_C_T", partner_pairs=[
            ("2_36932656_A_G", 0.94),
            ("2_36905984_C_T", 0.81),
            ("2_36897612_T_C", 0.40),
        ])
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(
        lead="2_36910110_C_T",
        partners=["2_36932656_A_G", "2_36905984_C_T", "2_36897612_T_C"],
        chromosome="2",
        window_bp=1_000_000,
    )

    assert result.lead_variant_id == "2_36910110_C_T"
    assert result.super_pop == "EUR"
    assert result.n_partners_requested == 3
    assert result.n_partners_returned == 3
    r2_by_id = {p.partner_variant_id: p.r2 for p in result.pairs}
    assert r2_by_id["2_36932656_A_G"] == pytest.approx(0.94)
    assert r2_by_id["2_36905984_C_T"] == pytest.approx(0.81)
    assert r2_by_id["2_36897612_T_C"] == pytest.approx(0.40)

    # Verify cmd shape uses plink 1.9 flags.
    cmd = captured_cmd["cmd"]
    assert "--r2" in cmd
    assert "--ld-snp" in cmd
    assert "--ld-window-r2" in cmd
    assert "--ld-window-kb" in cmd
    assert "--set-missing-var-ids" in cmd
    # plink 1.9 should NOT see plink2's matrix-only flag.
    assert "--r2-unphased" not in cmd


def test_r2_with_lead_passes_ids_to_plink_verbatim(monkeypatch, tmp_path):
    """The client only swaps `_` for `:`; it does not reorder alleles.

    plink 1.9's `--set-missing-var-ids '@:#:$1:$2'` spells the two alleles in
    ASCII order, so a caller id whose alleles are not in that order (`G_A`,
    `T_C`, ...) never matches the panel id. The extract list and the --ld-snp
    argument must carry the caller's spelling unchanged for that to be visible.
    """
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch, ["HG00096", "HG00097"])
    _patch_region_fetch(monkeypatch, tmp_path)
    seen: dict = {}

    def fake_run(cmd, capture_output, text, check, timeout):
        seen["ld_snp"] = cmd[cmd.index("--ld-snp") + 1]
        seen["extract"] = Path(cmd[cmd.index("--extract") + 1]).read_text().split()
        seen["keep"] = Path(cmd[cmd.index("--keep") + 1]).read_text().splitlines()
        seen["window_kb"] = cmd[cmd.index("--ld-window-kb") + 1]
        _write_ld(Path(cmd[cmd.index("--out") + 1]), lead="1_100_T_G", partner_pairs=[])
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    client.r2_with_lead(lead="1_100_T_G", partners=["1_50_G_A"], chromosome="chr1", window_bp=250_000)

    assert seen["ld_snp"] == "1:100:T:G"
    assert seen["extract"] == ["1:100:T:G", "1:50:G:A"]
    # plink 1.9 + --vcf sets FID = IID = sample id, so --keep rows repeat the id.
    assert seen["keep"] == ["HG00096\tHG00096", "HG00097\tHG00097"]
    assert seen["window_kb"] == "250"


def test_r2_with_lead_strips_chr_prefix_and_centres_the_window(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    fetched: dict = {}

    def fake_fetch(chrom, start, end, cache_dir):
        fetched.update(chrom=chrom, start=start, end=end)
        vcf = tmp_path / "region.vcf.gz"
        vcf.write_bytes(b"")
        return vcf

    monkeypatch.setattr(oc_module, "_fetch_region_vcf", fake_fetch)

    def fake_run(cmd, capture_output, text, check, timeout):
        _write_ld(Path(cmd[cmd.index("--out") + 1]), lead="1_1000_A_C", partner_pairs=[])
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)
    client = OnDemand1000GLDClient(super_pop="AFR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(lead="1_1000_A_C", partners=["1_1200_A_G"], chromosome="chr1", window_bp=400)
    assert fetched == {"chrom": "1", "start": 800, "end": 1200}
    assert result.chromosome == "1"
    assert result.window_bp == 400
    assert result.n_partners_returned == 0


def test_r2_with_lead_handles_empty_partners(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(
        lead="2_1_C_T", partners=[], chromosome="2", window_bp=1_000,
    )
    assert result.n_partners_requested == 0
    assert result.pairs == []
    assert any("no partners requested" in n for n in result.notes)


def test_r2_with_lead_rejects_unparseable_lead(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    with pytest.raises(OnDemandLDError, match="cannot parse lead variant id"):
        client.r2_with_lead(lead="rs646776", partners=["1_2_A_G"], chromosome="1", window_bp=1_000)


def test_r2_with_lead_raises_on_plink_nonzero_exit(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    def boom(cmd, capture_output, text, check, timeout):
        return MagicMock(returncode=1, stdout="", stderr="missing variant id")
    monkeypatch.setattr("subprocess.run", boom)

    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    with pytest.raises(OnDemandLDError, match="exited with code 1"):
        client.r2_with_lead(
            lead="2_1_C_T", partners=["2_2_A_G"], chromosome="2", window_bp=1_000,
        )


def test_r2_with_lead_raises_when_ld_absent(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    def succeeds_but_writes_nothing(cmd, capture_output, text, check, timeout):
        return MagicMock(returncode=0, stdout="ok", stderr="")
    monkeypatch.setattr("subprocess.run", succeeds_but_writes_nothing)

    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    with pytest.raises(OnDemandLDError, match=r"no \.ld output"):
        client.r2_with_lead(
            lead="2_1_C_T", partners=["2_2_A_G"], chromosome="2", window_bp=1_000,
        )


# -----------------------------------------------------------------
# Super-population sample resolution (panel TSV parsing, offline)
# -----------------------------------------------------------------


def test_resolve_super_pop_samples_reads_cached_panel_without_network(tmp_path, monkeypatch):
    panel = tmp_path / "integrated_call_samples_v3.20130502.ALL.panel"
    panel.write_text(
        "sample\tpop\tsuper_pop\tgender\n"
        "HG00096\tGBR\tEUR\tmale\n"
        "HG00097\tGBR\tEUR\tfemale\n"
        "NA19017\tLWK\tAFR\tfemale\n"
    )

    def no_network(*_a, **_k):
        raise AssertionError("panel is cached; requests.get must not be called")

    monkeypatch.setattr(oc_module.requests, "get", no_network)
    assert oc_module._resolve_super_pop_samples("EUR", tmp_path) == ["HG00096", "HG00097"]
    assert oc_module._resolve_super_pop_samples("AFR", tmp_path) == ["NA19017"]
    with pytest.raises(OnDemandLDError, match="Valid super-pop codes"):
        oc_module._resolve_super_pop_samples("XYZ", tmp_path)


# -----------------------------------------------------------------
# Parser-only tests
# -----------------------------------------------------------------


def test_parse_ld_skips_rows_without_lead_match(tmp_path):
    """An .ld row whose SNP_A and SNP_B are both unrelated to the lead is ignored."""
    path = tmp_path / "ld_out.ld"
    rows = [
        " CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2",
        "   2            1   2:x:C:T      2            2   2:y:A:G        0.500000",
        "   2            1         lead    2            2   2:z:A:G        0.700000",
    ]
    path.write_text("\n".join(rows) + "\n")

    notes: list[str] = []
    pairs = _parse_ld(path, "lead", notes)
    assert len(pairs) == 1
    assert pairs[0].partner_variant_id == "2:z:A:G"
    assert pairs[0].r2 == pytest.approx(0.7)


def test_parse_ld_keeps_the_lead_self_row_and_notes_na(tmp_path):
    """plink 1.9 emits the lead-vs-itself row; the parser keeps it as a pair (r2 = 1),
    so `n_partners_returned` counts it. An `NA` r2 is dropped with a note."""
    path = tmp_path / "ld_out.ld"
    rows = [
        " CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2",
        "   1          100   1:100:G:T      1          100   1:100:G:T        1.000000",
        "   1          100   1:100:G:T      1          150   1:150:C:T        0.204900",
        "   1          100   1:100:G:T      1          200   1:200:C:T              NA",
    ]
    path.write_text("\n".join(rows) + "\n")

    notes: list[str] = []
    pairs = _parse_ld(path, "1:100:G:T", notes)
    assert [p.partner_variant_id for p in pairs] == ["1:100:G:T", "1:150:C:T"]
    assert pairs[0].r2 == 1.0
    assert notes == ["missing r² for partner 1:200:C:T"]


# -----------------------------------------------------------------
# CLI (ld_1000g_region_compute.py main)
# -----------------------------------------------------------------


def _fake_result(lead: str, partners: list[str], super_pop: str = "EUR") -> OnDemandLDResult:
    pairs = [OnDemandLDPair(partner_variant_id=lead, r2=1.0)]
    pairs += [OnDemandLDPair(partner_variant_id=p, r2=0.5 + 0.1 * i) for i, p in enumerate(partners)]
    return OnDemandLDResult(
        panel_id="1000g_phase3_v5b_grch38_basic", panel_version="5b_remote_2019_03_12",
        super_pop=super_pop, plink_version="PLINK v1.90b6.27 64-bit (2023-05-09)",
        chromosome="1", lead_variant_id=lead, window_bp=1_000_000,
        n_partners_requested=len(partners), n_partners_returned=len(pairs),
        pairs=pairs, fetched_at_utc="2026-01-01T00:00:00Z", notes=["fetched 1000G region VCF to /x"],
    )


@pytest.fixture
def fake_client(monkeypatch, tmp_path):
    """Swap the on-demand client for a stub and point the result cache at tmp_path."""
    constructed: list[dict] = []

    class FakeClient:
        def __init__(self, super_pop, plink_bin):
            constructed.append({"super_pop": super_pop, "plink_bin": plink_bin})

        def r2_with_lead(self, lead, partners, chromosome, window_bp):
            return _fake_result(lead, list(partners), super_pop=constructed[-1]["super_pop"])

    monkeypatch.setattr(oc_module, "OnDemand1000GLDClient", FakeClient)
    monkeypatch.setattr(cli_module, "DEFAULT_LD_RESULT_CACHE_DIR", tmp_path / "result_cache")
    return constructed


def _read_manifest(out: Path) -> dict:
    """manifest.yaml when PyYAML is installed, else the JSON fallback the CLI writes."""
    if (out / "manifest.yaml").is_file():
        yaml = pytest.importorskip("yaml")
        return yaml.safe_load((out / "manifest.yaml").read_text())
    return json.loads((out / "manifest.json").read_text())


def test_cli_list_demos_names_the_bundled_configs(capsys):
    assert cli_module.main(["--list-demos"]) == 0
    out = capsys.readouterr().out
    assert "default (default)" in out
    assert "sort1_locus_eur" in out
    assert "input" in out


def test_cli_requires_output_and_a_config_source(tmp_path):
    with pytest.raises(SystemExit) as exc:
        cli_module.main(["--output", str(tmp_path / "o")])
    assert exc.value.code == 2
    with pytest.raises(SystemExit) as exc:
        cli_module.main(["--demo"])
    assert exc.value.code == 2


def test_cli_writes_pairs_manifest_and_report(fake_client, tmp_path, capsys):
    cfg = {
        "lead": "1_109274968_G_T",
        "partners": ["1_109274968_G_T", "1_109272630_A_G", "1_109274623_C_T"],
        "chromosome": 1,
        "window_bp": 1000000,
        "super_pop": "AFR",
        "plink_bin": "/opt/plink19/plink",
    }
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(cfg))
    out = tmp_path / "out"

    assert cli_module.main(["--input", str(cfg_path), "--output", str(out)]) == 0

    # The lead is removed from its own partner list; config overrides reach the client.
    assert fake_client == [{"super_pop": "AFR", "plink_bin": "/opt/plink19/plink"}]

    lines = (out / "ld_pairs.tsv").read_text().splitlines()
    assert lines[0] == "# ld-1000g-region-compute v0.1.0"
    assert lines[4].split("\t") == [
        "lead_variant_id", "partner_variant_id", "r2", "dprime", "panel_id", "super_pop",
    ]
    rows = [line.split("\t") for line in lines[5:]]
    assert [r[1] for r in rows] == ["1_109274968_G_T", "1_109272630_A_G", "1_109274623_C_T"]
    assert rows[0][2] == "1.000000" and rows[0][3] == ""
    assert {r[0] for r in rows} == {"1_109274968_G_T"}
    assert {r[5] for r in rows} == {"AFR"}

    manifest = _read_manifest(out)
    assert manifest["skill"] == "ld-1000g-region-compute"
    assert manifest["lead"] == "1_109274968_G_T"
    assert manifest["chromosome"] == "1"
    assert manifest["super_pop"] == "AFR"
    assert manifest["n_partners_requested"] == 2
    assert manifest["n_partners_returned"] == 3
    assert manifest["panel_version"] == "5b_remote_2019_03_12"
    assert manifest["outputs"] == {"ld_pairs_tsv": "ld_pairs.tsv"}

    report = (out / "report.md").read_text()
    assert "# ld-1000g-region-compute report" in report
    assert "chr1 ±500 kb" in report
    assert "2 / 3" in report
    assert "3 pairs ->" in capsys.readouterr().out


def test_cli_result_cache_hit_skips_the_client(fake_client, tmp_path, capsys):
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps({
        "lead": "1_109274968_G_T", "partners": ["1_109272630_A_G"],
        "chromosome": "1", "window_bp": 1000000,
    }))
    first = tmp_path / "first"
    second = tmp_path / "second"
    assert cli_module.main(["--input", str(cfg_path), "--output", str(first)]) == 0
    assert cli_module.main(["--input", str(cfg_path), "--output", str(second)]) == 0
    assert len(fake_client) == 1, "second run must be served from the result cache"
    assert "(cache hit)" in capsys.readouterr().out
    assert (second / "ld_pairs.tsv").read_text() == (first / "ld_pairs.tsv").read_text()
    cached = list((tmp_path / "result_cache").glob("*.json"))
    assert len(cached) == 1
    assert cached[0].name.startswith("1_109274968_G_T__EUR__win1000000__")

    third = tmp_path / "third"
    assert cli_module.main(["--input", str(cfg_path), "--output", str(third), "--no-cache"]) == 0
    assert len(fake_client) == 2, "--no-cache must bypass the result cache"


def test_cli_demo_resolves_the_bundled_default(fake_client, tmp_path, capsys):
    out = tmp_path / "demo"
    assert cli_module.main(["--demo", "--output", str(out)]) == 0
    assert "using bundled demo default.json" in capsys.readouterr().err
    manifest = _read_manifest(out)
    assert manifest["lead"] == "1_109274968_G_T"
    assert manifest["n_partners_requested"] == 5


def test_cli_demo_unknown_name_lists_the_available_ones(fake_client, tmp_path):
    with pytest.raises(FileNotFoundError, match="sort1_locus_eur"):
        cli_module.main(["--demo", "nope", "--output", str(tmp_path / "o")])


def test_cli_rejects_unknown_config_extension(tmp_path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text("lead = 'x'")
    with pytest.raises(ValueError, match="unsupported config extension"):
        cli_module.main(["--input", str(cfg), "--output", str(tmp_path / "o")])
