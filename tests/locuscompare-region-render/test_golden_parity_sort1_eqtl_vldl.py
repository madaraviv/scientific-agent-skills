"""Golden-parity test: SORT1 eQTL x cholesterol-VLDL, replayed offline.

Two cassettes captured 2026-05-15 (eQTL Catalogue QTD000276 for ENSG00000134243
and GWAS Catalog harmonised GCST90269602, both +/-500 kb of 1_109274968_G_T)
are served by stand-in clients, so the composer's join, allele harmonisation,
palindromic flagging, caveat and manifest logic run on real-shaped data with no
network. The manifest block is compared with `fixtures/golden/sort1_eqtl_vldl/
expected.yaml`. A failing assertion means the composer changed behaviour.

Needs the sibling skills (the cassettes are rebuilt into their RegionResult
dataclasses); skips with the composer's ImportError otherwise.
"""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

try:
    from locuscompare_region_render import (  # noqa: E402
        EXPOSURE_KIND_EQTL_CATALOGUE,
        LocusCompareSpec,
        Tier2NotAvailable,
        render_locuscompare_for_lead,
    )
except ImportError as _e:
    pytest.skip(str(_e), allow_module_level=True)

from eqtl_catalogue_region_fetch import EQTLCatalogueRelease  # noqa: E402
from eqtl_catalogue_region_fetch import RegionResult as EQTLRegionResult  # noqa: E402
from eqtl_catalogue_region_fetch import RegionVariant as EQTLRegionVariant  # noqa: E402
from gwas_catalog_region_fetch import GWASCatalogRelease  # noqa: E402
from gwas_catalog_region_fetch import RegionResult as GWASRegionResult  # noqa: E402
from gwas_catalog_region_fetch import RegionVariant as GWASRegionVariant  # noqa: E402
from ld_1000g_region_compute import SuperPop  # noqa: E402

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "golden" / "sort1_eqtl_vldl"
INPUTS_DIR = FIXTURE_DIR / "inputs"


def _load_eqtl_cassette() -> EQTLRegionResult:
    with gzip.open(INPUTS_DIR / "eqtl_exposure_QTD000276_chr1_109_274_968_window.json.gz", "rt") as f:
        data = json.load(f)
    variants = [
        EQTLRegionVariant(**{k: v for k, v in row.items() if k != "raw"} | {"raw": row.get("raw", {})})
        for row in data["variants"]
    ]
    return EQTLRegionResult(
        dataset_id=data["dataset_id"], chromosome=data["chromosome"],
        region_start_bp=data["region_start_bp"], region_end_bp=data["region_end_bp"],
        n_variants=data["n_variants"], variants=variants,
        release=EQTLCatalogueRelease(**data["release"]), notes=list(data.get("notes", [])),
    )


def _load_gwas_cassette() -> GWASRegionResult:
    with gzip.open(INPUTS_DIR / "gwas_outcome_GCST90269602_chr1_109_274_968_window.json.gz", "rt") as f:
        data = json.load(f)
    variants = [
        GWASRegionVariant(**{k: v for k, v in row.items() if k != "raw"} | {"raw": row.get("raw", {})})
        for row in data["variants"]
    ]
    return GWASRegionResult(
        accession=data["accession"], chromosome=data["chromosome"],
        region_start_bp=data["region_start_bp"], region_end_bp=data["region_end_bp"],
        n_variants=data["n_variants"], variants=variants,
        release=GWASCatalogRelease(**data["release"]), notes=list(data.get("notes", [])),
    )


class _CassetteClient:
    """Returns the cassette for any fetch_region call and records the call."""

    def __init__(self, cassette) -> None:
        self._cassette = cassette
        self.calls: list[dict[str, Any]] = []

    def fetch_region(self, **kwargs: Any):
        self.calls.append(kwargs)
        return self._cassette


def _eqtl_spec() -> LocusCompareSpec:
    return LocusCompareSpec(
        lead_variant_id="1_109274968_G_T", chromosome="1", lead_position_bp=109_274_968,
        window_bp=1_000_000, eqtl_dataset_id="QTD000276", molecular_trait_id="ENSG00000134243",
        gwas_accession="GCST90269602", exposure_kind=EXPOSURE_KIND_EQTL_CATALOGUE,
        exposure_gene_symbol="SORT1", outcome_trait_label="cholesterol in medium VLDL",
        release_tag="26.03", super_pop=SuperPop.EUR,
        prefetched_gene_track=[],  # no Ensembl call; the caveat is fine for the assertion
    )


@pytest.fixture(scope="module")
def expected_block() -> dict[str, Any]:
    return yaml.safe_load((FIXTURE_DIR / "expected.yaml").read_text())["eqtl_render"]["manifest_block"]


def _render(tmp_path, *, eqtl_cassette=None, gwas_cassette=None):
    eqtl_client = _CassetteClient(eqtl_cassette or _load_eqtl_cassette())
    gwas_client = _CassetteClient(gwas_cassette or _load_gwas_cassette())
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(), eqtl_client=eqtl_client, gwas_client=gwas_client,
        ld_client=None, out_path=tmp_path / "sort1_eqtl.png", ukb_ppp_client=None,
    )
    return result, eqtl_client, gwas_client


def test_sort1_eqtl_vldl_replay_matches_the_locked_manifest_block(tmp_path, expected_block):
    result, eqtl_client, gwas_client = _render(tmp_path)
    diffs = [f"  {k}: got {result.manifest_block.get(k)!r}, expected {v!r}"
             for k, v in expected_block.items() if result.manifest_block.get(k) != v]
    assert not diffs, "manifest_block drifted from golden expected.yaml:\n" + "\n".join(diffs)
    assert result.n_pairs == expected_block["n_pairs"]
    assert result.n_palindromic_excluded == expected_block["n_palindromic_excluded"]
    assert result.plot_path.is_file() and result.plot_path.stat().st_size > 10_000
    # The composer asks the exposure fetcher for the parent gene via `gene_id`
    # and both fetchers for the same +/-500 kb window.
    assert len(eqtl_client.calls) == 1 and len(gwas_client.calls) == 1
    assert eqtl_client.calls[0]["gene_id"] == "ENSG00000134243"
    assert gwas_client.calls[0]["accession"] == "GCST90269602"
    for call in (eqtl_client.calls[0], gwas_client.calls[0]):
        assert (call["start_bp"], call["end_bp"]) == (108_774_968, 109_774_968)


def test_golden_block_key_set_is_the_manifest_key_set_minus_the_two_volatile_keys(tmp_path, expected_block):
    """The fixture must name every stable key, so a new manifest field cannot ship
    unlocked and a renamed one cannot vanish from the comparison."""
    result, _, _ = _render(tmp_path)
    assert set(result.manifest_block) - set(expected_block) == {"fetched_at", "plot_artifact"}


def _caveat_marker(qm: str) -> str:
    return f"credible-set-filtered (eQTL Catalogue .cc.tsv.gz; quant_method={qm})"


@pytest.mark.parametrize("quant_method", ["txrev", "leafcutter", "exon", "tx"])
def test_non_ge_quant_methods_add_the_credible_set_caveat(tmp_path, quant_method):
    cassette = _load_eqtl_cassette()
    cassette.release.quant_method = quant_method
    result, _, _ = _render(tmp_path, eqtl_cassette=cassette)
    caveats = result.manifest_block["ancestry_caveats"]
    assert any(_caveat_marker(quant_method) in c for c in caveats), caveats


@pytest.mark.parametrize("quant_method", ["ge", "microarray", ""])
def test_full_nominal_pass_quant_methods_omit_the_credible_set_caveat(tmp_path, quant_method):
    cassette = _load_eqtl_cassette()
    cassette.release.quant_method = quant_method or None
    result, _, _ = _render(tmp_path, eqtl_cassette=cassette)
    assert not any("credible-set-filtered" in c for c in result.manifest_block["ancestry_caveats"])


def test_fetcher_notes_propagate_to_data_source_warnings_not_caveats(tmp_path):
    """Each fetcher's RegionResult.notes lands in `data_source_warnings` with the
    fetcher name as prefix, and never in the curated `ancestry_caveats`."""
    eqtl_cassette = _load_eqtl_cassette()
    eqtl_cassette.notes = ["row column count 18 != header 19; skipping"]
    gwas_cassette = _load_gwas_cassette()
    gwas_cassette.notes = ["row column count 14 != header 15; skipping"]
    result, _, _ = _render(tmp_path, eqtl_cassette=eqtl_cassette, gwas_cassette=gwas_cassette)
    warnings = result.manifest_block["data_source_warnings"]
    assert "eqtl_catalogue: row column count 18 != header 19; skipping" in warnings
    assert "gwas_catalog: row column count 14 != header 15; skipping" in warnings
    assert not any(c.startswith(("eqtl_catalogue:", "gwas_catalog:", "ukb_ppp:"))
                   for c in result.manifest_block["ancestry_caveats"])


def test_zero_exposure_variants_is_a_tier2_not_available(tmp_path):
    cassette = _load_eqtl_cassette()
    cassette.variants = []
    cassette.n_variants = 0
    with pytest.raises(Tier2NotAvailable, match="zero variants"):
        _render(tmp_path, eqtl_cassette=cassette)


def test_no_joinable_variants_is_a_tier2_not_available(tmp_path):
    """Exposure and outcome sharing no variant_id cannot be plotted; the composer
    refuses by name instead of rendering empty panels."""
    gwas_cassette = _load_gwas_cassette()
    for v in gwas_cassette.variants:
        v.variant_id = "X_" + v.variant_id
    with pytest.raises(Tier2NotAvailable, match="no joinable variants"):
        _render(tmp_path, gwas_cassette=gwas_cassette)
