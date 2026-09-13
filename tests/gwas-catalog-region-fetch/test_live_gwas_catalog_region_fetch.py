"""Live smoke test for gwas-catalog-region-fetch.

Reads a small slice of the real harmonised file for GCST90269602 (cholesterol in medium
VLDL, 1 Mb around SORT1) from the EBI GWAS Catalog FTP via tabix. Gated by
@pytest.mark.live and RUN_LIVE_TESTS=1 so the offline suite stays network-free.

Run with:
    RUN_LIVE_TESTS=1 python -m pytest tests/gwas-catalog-region-fetch -m live
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "gwas-catalog-region-fetch"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

from gwas_catalog_region_fetch import (  # noqa: E402
    GWASCatalogClient,
    RegionResult,
)

pytestmark = pytest.mark.live


@pytest.mark.skipif(
    os.getenv("RUN_LIVE_TESTS") != "1",
    reason="live tests gated on RUN_LIVE_TESTS=1",
)
def test_live_fetch_region_sort1_cholesterol_smoke(tmp_path, monkeypatch) -> None:
    """A real tabix slice of the SORT1 window returns variants with sane fields.

    pysam writes the remote .tbi into the working directory, so the test runs from
    tmp_path to keep the repository clean.
    """
    monkeypatch.chdir(tmp_path)
    client = GWASCatalogClient()
    result = client.fetch_region(
        accession="GCST90269602",
        chromosome="1",
        start_bp=109_774_000,
        end_bp=109_775_000,
    )
    assert isinstance(result, RegionResult)
    assert result.variants, "expected a non-empty tabix slice for the SORT1 window"
    assert result.release.harmonised_url.endswith(
        "/GCST90269001-GCST90270000/GCST90269602/harmonised/GCST90269602.h.tsv.gz"
    )
    for v in result.variants:
        assert v.chromosome == "1"
        assert 109_774_000 <= v.position <= 109_775_000
        assert v.variant_id == f"1_{v.position}_{v.ref}_{v.alt}"
        assert v.se is None or v.se > 0
        assert v.p_value is None or 0 <= v.p_value <= 1
    assert (tmp_path / "GCST90269602.h.tsv.gz.tbi").is_file(), "htslib's index download lands in cwd"
