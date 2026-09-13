"""pytest configuration for the locuscompare-region-render suite.

Two jobs. First, the `live` marker (network / matplotlib smoke), gated on
RUN_LIVE_TESTS=1 or `pytest -m live`. Second, sibling-skill resolution: the
composer imports eqtl-catalogue-region-fetch, gwas-catalog-region-fetch and
ld-1000g-region-compute from `skills/<sibling>/scripts/`. When those skills are
installed in this checkout nothing more is needed. When they are not (this
branch on its own), set LOCUSCOMPARE_SIBLING_SKILLS_ROOT to a directory holding
the same `<sibling>/scripts/` layout and the suite runs against those copies;
without either, the tests that need the composer skip with the composer's own
ImportError message as the reason. No sibling is ever mocked.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"

SIBLING_SKILLS = (
    "eqtl-catalogue-region-fetch",
    "gwas-catalog-region-fetch",
    "ld-1000g-region-compute",
    "ukb-ppp-region-fetch",  # optional; the composer loads without it
)


def _sibling_roots() -> list[Path]:
    roots = [SKILL_ROOT.parent]
    override = os.environ.get("LOCUSCOMPARE_SIBLING_SKILLS_ROOT")
    if override:
        roots.append(Path(override))
    return roots


for _root in _sibling_roots():
    for _skill in SIBLING_SKILLS:
        _scripts = _root / _skill / "scripts"
        if _scripts.is_dir() and str(_scripts) not in sys.path:
            sys.path.append(str(_scripts))


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "live: smoke test that runs the full matplotlib render pipeline end to end",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    run_live = os.environ.get("RUN_LIVE_TESTS") == "1" or "live" in (config.getoption("-m") or "")
    if run_live:
        return
    skip_live = pytest.mark.skip(reason="live test (set RUN_LIVE_TESTS=1 or pytest -m live)")
    for item in items:
        if "live" in item.keywords:
            item.add_marker(skip_live)
