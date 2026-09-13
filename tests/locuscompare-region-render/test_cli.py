"""Unit tests for the `cli.py` helpers: demo discovery and selection, relative-path
resolution against the config file's directory, and the provenance prefix. None of
these need the sibling skills; `cli.py` defers those imports to `_run`."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

import cli  # noqa: E402

NUMBERED_DEMOS = (
    "01_synthetic_demo",
    "02_eqtl_catalogue_x_gwas_catalog",
    "03_open_targets_followup",
    "04_gwas_lookup_followup",
    "05_sqtl_sort1_liver_txrev",
    "06_sceqtl_sort1_onek1k_cd14_mono",
    "07_pqtl_sort1_ukbppp_eur",
)


def test_list_demos_returns_the_seven_numbered_configs_and_skips_docs_dirs():
    paths = cli._list_demos()
    assert [p.parent.name for p in paths] == list(NUMBERED_DEMOS)
    names = {p.parent.name for p in paths}
    assert "recipes" not in names and "chains" not in names
    for p in paths:
        assert p.is_file() and p.name.startswith("config.")


def test_resolve_demo_path_default_picks_eqtl_x_gwas_demo():
    p = cli._resolve_demo_path("__default__")
    assert p.parent.name == "02_eqtl_catalogue_x_gwas_catalog"
    assert p.name.startswith("config.")


def test_resolve_demo_path_named_demo():
    assert cli._resolve_demo_path("01_synthetic_demo").parent.name == "01_synthetic_demo"


def test_resolve_demo_path_short_prefix_match():
    """`--demo 01` maps to `01_synthetic_demo/config.*`."""
    assert cli._resolve_demo_path("01").parent.name == "01_synthetic_demo"


def test_resolve_demo_path_unknown_demo_lists_available():
    with pytest.raises(FileNotFoundError) as excinfo:
        cli._resolve_demo_path("99_does_not_exist")
    msg = str(excinfo.value)
    assert "no bundled demo named" in msg
    assert "01_synthetic_demo" in msg


def test_print_available_demos_lists_each_numbered_demo(capsys):
    cli._print_available_demos()
    out = capsys.readouterr().out
    for name in NUMBERED_DEMOS:
        assert name in out
    assert "(default)" in out
    assert "[config." in out


def test_resolve_path_returns_absolute_unchanged(tmp_path):
    abs_path = tmp_path / "abs.tsv"
    abs_path.write_text("")
    assert cli._resolve_path(abs_path, tmp_path / "other_config_dir") == abs_path


def test_resolve_path_resolves_relative_against_config_dir(tmp_path):
    cfg_dir = tmp_path / "cfg"
    cfg_dir.mkdir()
    target = cfg_dir / "rel.tsv"
    target.write_text("")
    assert cli._resolve_path("rel.tsv", cfg_dir) == target.resolve()


def test_cli_list_demos_flag_returns_zero(capsys):
    assert cli.main(["--list-demos"]) == 0
    assert "01_synthetic_demo" in capsys.readouterr().out


def test_format_provenance_prefix_combines_ot_and_gwas_lookup():
    prefix = cli._format_provenance_prefix(
        {"provenance": {"ot_release": "26.03", "gwas_lookup_run_dir": "runs/sort1_vldl"}}
    )
    assert "OT release: 26.03" in prefix
    assert "gwas-lookup chain: runs/sort1_vldl" in prefix
    assert prefix.endswith(" | ")
