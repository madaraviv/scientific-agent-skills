"""Live smoke test: real plink 1.9 against a 1000G region tabix-fetched from EBI.

Gated by RUN_LIVE_TESTS=1 (or `pytest -m live`) and skipped when plink is not on
PATH (the binary is GPL-3 and not packaged via pip; install via brew / apt / conda,
or set PLINK_BIN). Network-dependent by design; not part of every-PR CI.

The unit suite mocks both the fetch and plink, so it cannot notice the panel
moving or plink changing its output layout; this test can.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "ld-1000g-region-compute"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

from ondemand_client import (  # noqa: E402
    OnDemand1000GLDClient,
    OnDemandLDResult,
)

pytestmark = pytest.mark.live


def _plink_bin() -> str | None:
    explicit = os.getenv("PLINK_BIN")
    if explicit and os.path.isfile(explicit):
        return explicit
    return shutil.which("plink")


@pytest.mark.skipif(
    _plink_bin() is None,
    reason="plink not on PATH and PLINK_BIN not set",
)
def test_live_r2_with_lead_sort1_smoke() -> None:
    """The bundled SORT1 demo lead against its three ASCII-ordered partners, EUR.

    The other two demo partners (`1_109270398_G_A`, `1_109274857_G_C`) are spelled
    with alleles out of ASCII order and are dropped by plink's id rewrite, which is
    the documented gotcha rather than a panel gap; they are left out here so the
    assertion is about the fetch and the compute.
    """
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin=_plink_bin())
    lead = "1_109274968_G_T"
    partners = ["1_109272630_A_G", "1_109274570_A_G", "1_109274623_C_T"]
    result = client.r2_with_lead(
        lead=lead, partners=partners, chromosome="1", window_bp=1_000_000,
    )
    assert isinstance(result, OnDemandLDResult)
    assert result.plink_version.startswith("PLINK v1.9")
    r2_by_id = {p.partner_variant_id: p.r2 for p in result.pairs}
    # plink emits the lead-vs-itself row and the parser keeps it.
    assert r2_by_id[lead] == pytest.approx(1.0)
    for partner in partners:
        assert partner in r2_by_id, f"{partner} missing from the 1000G EUR panel result"
        assert 0.0 <= r2_by_id[partner] <= 1.0
    assert result.n_partners_returned == len(partners) + 1
