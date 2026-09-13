"""Unit tests for the `_prefetched.py` stand-in clients (PrefetchedEQTLClient /
PrefetchedGWASClient / PrefetchedLDClient), which satisfy the composer's fetcher
interface from TSVs, plus the `gwas_variants_from_eqtl` reshape. These construct
the sibling skills' RegionResult dataclasses, so they need the siblings."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

try:
    import locuscompare_region_render  # noqa: E402,F401  (resolves the siblings)
except ImportError as _e:
    pytest.skip(str(_e), allow_module_level=True)

from _prefetched import (  # noqa: E402
    PrefetchedEQTLClient,
    PrefetchedGWASClient,
    PrefetchedLDClient,
    gwas_variants_from_eqtl,
)
from eqtl_catalogue_region_fetch import RegionVariant as EQTLRegionVariant  # noqa: E402
from gwas_catalog_region_fetch import RegionVariant as GWASRegionVariant  # noqa: E402


def _eqtl_variant(vid="1_100_A_G", ref="A", alt="G", beta=0.1, se=0.01, p=1e-6):
    return EQTLRegionVariant(
        variant_id=vid, chromosome="1", position=int(vid.split("_")[1]),
        ref=ref, alt=alt, beta=beta, se=se, p_value=p,
        maf=0.2, effect_allele_frequency=0.2, raw={},
    )


def test_prefetched_eqtl_client_fetch_region_returns_loaded_variants():
    variants = [_eqtl_variant("1_100_A_G"), _eqtl_variant("1_200_C_T", ref="C", alt="T")]
    client = PrefetchedEQTLClient(
        variants=variants, dataset_id="QTD000276", dataset_release="prefetched",
        study_label="SORT1 expression - minor salivary gland",
        tissue_label="minor salivary gland", quant_method="ge",
    )
    result = client.fetch_region(dataset_id="ignored-passthrough", chromosome="1",
                                 start_bp=1, end_bp=1_000_000)
    assert result.dataset_id == "QTD000276"
    assert result.n_variants == 2
    assert result.region_start_bp == 1 and result.region_end_bp == 1_000_000
    assert result.variants[0].variant_id == "1_100_A_G"
    assert any("pre-fetched" in n for n in result.notes)
    assert result.release.dataset_release == "prefetched"
    assert result.release.tissue_label == "minor salivary gland"


def test_prefetched_eqtl_client_falls_back_to_caller_dataset_id_when_empty():
    client = PrefetchedEQTLClient(variants=[_eqtl_variant()], dataset_id="")
    result = client.fetch_region(dataset_id="caller-supplied", chromosome="1",
                                 start_bp=1, end_bp=1_000)
    assert result.dataset_id == "caller-supplied"


def test_prefetched_gwas_client_fetch_region_returns_loaded_variants():
    gwas_variants = gwas_variants_from_eqtl([_eqtl_variant("1_100_A_G")])
    client = PrefetchedGWASClient(variants=gwas_variants, accession="GCST90269602")
    result = client.fetch_region(accession="ignored-passthrough", chromosome="1",
                                 start_bp=1, end_bp=1_000_000)
    assert result.accession == "GCST90269602"
    assert result.n_variants == 1
    assert result.variants[0].variant_id == "1_100_A_G"
    assert any("pre-fetched" in n for n in result.notes)


def test_gwas_variants_from_eqtl_preserves_core_fields():
    eqtl = [
        _eqtl_variant("1_100_A_G", ref="A", alt="G", beta=0.5, se=0.05, p=1e-10),
        _eqtl_variant("1_200_C_T", ref="C", alt="T", beta=-0.3, se=0.04, p=1e-8),
    ]
    gwas = gwas_variants_from_eqtl(eqtl)
    assert len(gwas) == 2
    assert all(isinstance(v, GWASRegionVariant) for v in gwas)
    assert [v.variant_id for v in gwas] == ["1_100_A_G", "1_200_C_T"]
    assert gwas[0].beta == pytest.approx(0.5)
    assert gwas[1].beta == pytest.approx(-0.3)
    assert gwas[0].odds_ratio is None and gwas[1].odds_ratio is None


def test_prefetched_ld_client_returns_only_known_partners():
    """Partner ids missing from the matrix are dropped and the counts say so."""
    client = PrefetchedLDClient(
        r2_by_partner={"1_100_A_G": 0.95, "1_200_C_T": 0.40},
        panel_id="synthetic-test", panel_version="1.0",
        plink_version="prefetched", super_pop="synthetic",
    )
    result = client.r2_with_lead(lead="1_500_A_T",
                                 partners=["1_100_A_G", "1_999_T_C", "1_200_C_T"],
                                 chromosome="1", window_bp=1_000_000)
    assert result.lead_variant_id == "1_500_A_T"
    assert result.window_bp == 1_000_000
    assert result.panel_id == "synthetic-test"
    assert result.plink_version == "prefetched"
    assert result.n_partners_requested == 3
    assert result.n_partners_returned == 2
    r2_by_id = {p.partner_variant_id: p.r2 for p in result.pairs}
    assert r2_by_id == {"1_100_A_G": pytest.approx(0.95), "1_200_C_T": pytest.approx(0.40)}
    assert any("pre-fetched synthetic matrix" in n for n in result.notes)


def test_prefetched_ld_client_handles_empty_partner_list():
    client = PrefetchedLDClient(r2_by_partner={"1_100_A_G": 0.9})
    result = client.r2_with_lead(lead="1_500_A_T", partners=[], chromosome="1", window_bp=1_000_000)
    assert result.n_partners_requested == 0
    assert result.n_partners_returned == 0
    assert result.pairs == []
