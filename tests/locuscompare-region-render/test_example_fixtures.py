"""The recorded runs of demos 05 (sqtl, txrev) and 07 (pQTL, UKB-PPP) under
fixtures/examples/<demo>/expected_output/ are checked against the current code
rather than kept on trust. Two things can drift: the manifest shape cli.py
writes (top-level key order, render_block key set) and the report's lead line,
whose window label is the half-window the fetch covered. The demo READMEs quote
the recorded counts, so those are checked against the fixtures too."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
import yaml

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

import cli as locuscompare  # noqa: E402

try:
    import locuscompare_region_render  # noqa: E402,F401  (resolves the siblings)
    SIBLING_ERROR = None
except ImportError as _e:
    SIBLING_ERROR = str(_e)

needs_siblings = pytest.mark.skipif(SIBLING_ERROR is not None, reason=SIBLING_ERROR or "")

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "examples"
EXAMPLES = SKILL_ROOT / "scripts" / "examples"
RECORDED_DEMOS = ("05_sqtl_sort1_liver_txrev", "07_pqtl_sort1_ukbppp_eur")


def _lead_line(manifest: dict) -> str:
    """The report's lead line as cli.py writes it: chromosome and position come
    from the variant id, the rs id is optional, and the window label is the
    half-window (window_bp // 2000 kb), never the full width."""
    lead = manifest["lead_variant_id"]
    chromosome, position_bp = lead.split("_")[:2]
    rs_id = manifest["lead_rs_id"]
    window_bp = manifest["render_block"]["window_bp"]
    return (
        f"- **Lead variant:** `{lead}`"
        + (f" ({rs_id}; " if rs_id else " (")
        + f"chr{chromosome}:{position_bp}, ±{window_bp//2000} kb)"
    )


def _load(demo: str) -> tuple[dict, list[str]]:
    out = FIXTURES / demo / "expected_output"
    manifest = yaml.safe_load((out / "manifest.yaml").read_text())
    report = (out / "report.md").read_text().splitlines()
    return manifest, report


@pytest.mark.parametrize("demo", RECORDED_DEMOS)
def test_recorded_report_lead_line_carries_the_half_window(demo: str):
    manifest, report = _load(demo)
    assert manifest["render_block"]["window_bp"] == 1_000_000
    assert report[0] == "# locuscompare report" and report[1] == ""
    assert report[2] == _lead_line(manifest)
    assert "±500 kb" in report[2] and "1000 kb" not in report[2]
    assert f"- **n_pairs:** {manifest['n_pairs']}" in report
    assert f"- **n_palindromic_excluded:** {manifest['n_palindromic_excluded']}" in report


@pytest.mark.parametrize("demo", RECORDED_DEMOS)
def test_demo_readme_quotes_the_recorded_counts(demo: str):
    manifest, _ = _load(demo)
    readme = (EXAMPLES / demo / "README.md").read_text()
    quoted = {k: int(v) for k, v in re.findall(r"`(n_pairs|n_palindromic_excluded): (\d+)`", readme)}
    assert quoted == {"n_pairs": manifest["n_pairs"],
                      "n_palindromic_excluded": manifest["n_palindromic_excluded"]}
    assert f"fixtures/examples/{demo}/expected_output/" in readme
    assert "test_example_fixtures.py" in readme


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[dict, list[str]]:
    """One offline emit of 01_synthetic_demo: the manifest shape and report the
    current cli.py writes, to compare the recorded runs against."""
    out = tmp_path_factory.mktemp("synthetic") / "out"
    rc = locuscompare.main(["--demo", "01_synthetic_demo", "--output", str(out)])
    assert rc == 0, "01_synthetic_demo failed offline"
    manifest = yaml.safe_load((out / "manifest.yaml").read_text())
    return manifest, (out / "report.md").read_text().splitlines()


@needs_siblings
def test_lead_line_template_matches_what_cli_writes(synthetic_run):
    """The template `_lead_line` reproduces is the one cli.py emits: a fresh run
    of the synthetic demo has to match it before it is trusted on the fixtures."""
    manifest, report = synthetic_run
    assert report[2] == _lead_line(manifest)


@needs_siblings
@pytest.mark.parametrize("demo", RECORDED_DEMOS)
def test_recorded_manifest_shape_matches_a_fresh_emit(demo: str, synthetic_run):
    fresh, _ = synthetic_run
    recorded, _ = _load(demo)
    assert list(recorded) == list(fresh), "top-level key order differs from what cli.py writes"
    assert set(recorded["render_block"]) == set(fresh["render_block"]), (
        "render_block key set differs from the composer's manifest_block")
    assert recorded["render_block"]["lead_variant_id"] == recorded["lead_variant_id"]
    assert recorded["render_block"]["n_pairs"] == recorded["n_pairs"]
