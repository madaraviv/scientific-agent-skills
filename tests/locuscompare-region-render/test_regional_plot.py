"""Harmoniser and renderer tests for `regional_plot.py`: allele flip / palindromic
flag / drop rules of `harmonise_regions_for_locuscompare`, and
`render_full_locuscompare` end to end with and without a gene track. The renderer
is pure (fully resolved input in, PNG out), so no sibling skill is needed."""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

from regional_plot import (  # noqa: E402
    LD_R2_BINS,
    GeneTrackEntry,
    HarmonisedRegionPair,
    RegionalLocusCompareInput,
    _r2_color,
    _safe_neg_log10,
    _stratified_downsample,
    harmonise_regions_for_locuscompare,
    render_full_locuscompare,
)
from _wald_ratio_types import LocusVariant  # noqa: E402


def _lv(vid, ref, alt, beta, se=0.05, p=1e-8, pip=0.5):
    return LocusVariant(
        variant_id=vid, chromosome="1", position=int(vid.split("_")[1]),
        ref=ref, alt=alt, pip=pip, beta=beta, se=se, p_value=p,
        is95=True, is99=True, r2_lead=1.0,
    )


def _png_pixel_dims(path: Path) -> tuple[int, int]:
    """(width, height) from the PNG IHDR chunk, without Pillow."""
    raw = path.read_bytes()
    assert raw.startswith(b"\x89PNG\r\n\x1a\n"), "not a PNG file"
    return struct.unpack(">I", raw[16:20])[0], struct.unpack(">I", raw[20:24])[0]


def _pair(vid, pos, ref, alt, be, bo, pe, po, r2):
    return HarmonisedRegionPair(
        variant_id=vid, chromosome="1", position=pos, ref=ref, alt=alt,
        beta_exposure=be, se_exposure=0.04, p_exposure=pe,
        beta_outcome=bo, se_outcome=0.05, p_outcome=po, r2_with_lead=r2,
    )


def _stub_input(*, lead_id="1_500_A_G", with_gene_track=True):
    pairs = [
        _pair(lead_id, 500, "A", "G", 0.5, 0.4, 1e-12, 1e-10, 1.0),
        _pair("1_400_C_T", 400, "C", "T", 0.30, 0.20, 1e-6, 1e-5, 0.65),
        _pair("1_600_C_T", 600, "C", "T", -0.10, -0.05, 1e-2, 1e-1, 0.25),
    ]
    exposure_track = [_lv(p.variant_id, p.ref, p.alt, p.beta_exposure) for p in pairs]
    outcome_track = [_lv(p.variant_id, p.ref, p.alt, p.beta_outcome) for p in pairs]
    gene_track = []
    if with_gene_track:
        gene_track = [
            GeneTrackEntry(gene_symbol="DEMOGENE_A", start=100, end=300, strand="+"),
            GeneTrackEntry(gene_symbol="SORT1", start=450, end=550, strand="-"),
            GeneTrackEntry(gene_symbol="DEMOGENE_C", start=700, end=900, strand="+"),
        ]
    return RegionalLocusCompareInput(
        pairs=pairs, lead_variant_id=lead_id, chromosome="1", window_bp=1000,
        ld_panel_label="synthetic LD (test)", window_label=f"+/-500 bp of lead {lead_id}",
        exposure_label="synthetic exposure", outcome_label="synthetic outcome",
        provenance_label="test provenance", caveats=["test caveat"],
        title="regional_plot smoke",
        exposure_track_variants=exposure_track, outcome_track_variants=outcome_track,
        r2_by_variant={p.variant_id: p.r2_with_lead or 0.0 for p in pairs},
        exposure_short_label="synthetic eQTL", outcome_short_label="synthetic GWAS",
        gene_track=gene_track, focal_gene_symbol="SORT1" if with_gene_track else None,
    )


# ----- harmonise_regions_for_locuscompare


def test_harmonise_regions_keeps_matching_alleles_unchanged():
    pairs = harmonise_regions_for_locuscompare(
        [_lv("1_100_A_G", "A", "G", 0.5)], [_lv("1_100_A_G", "A", "G", 0.3)],
        {"1_100_A_G": 0.9}, "1_100_A_G",
    )
    assert len(pairs) == 1
    p = pairs[0]
    assert p.flip_outcome_beta is False
    assert p.palindromic_excluded is False
    assert p.beta_outcome == pytest.approx(0.3)
    assert p.r2_with_lead == 1.0  # the lead is always seeded to 1.0


def test_harmonise_regions_flips_swapped_outcome_beta():
    pairs = harmonise_regions_for_locuscompare(
        [_lv("1_100_A_G", "A", "G", 0.5)], [_lv("1_100_A_G", "G", "A", 0.3)],
        {"1_100_A_G": 0.9}, "1_999_T_C",
    )
    assert len(pairs) == 1
    assert pairs[0].flip_outcome_beta is True
    assert pairs[0].beta_outcome == pytest.approx(-0.3)


def test_harmonise_regions_flags_palindromic_without_dropping():
    pairs = harmonise_regions_for_locuscompare(
        [_lv("1_100_A_T", "A", "T", 0.5)], [_lv("1_100_A_T", "A", "T", 0.3)], {}, "1_999_T_C",
    )
    assert len(pairs) == 1
    assert pairs[0].palindromic_excluded is True


def test_harmonise_regions_drops_irreconcilable_alleles():
    pairs = harmonise_regions_for_locuscompare(
        [_lv("1_100_A_G", "A", "G", 0.5)], [_lv("1_100_A_G", "C", "T", 0.3)], {}, "1_999_T_C",
    )
    assert pairs == []


def test_harmonise_regions_drops_variants_only_on_one_side():
    pairs = harmonise_regions_for_locuscompare(
        [_lv("1_100_A_G", "A", "G", 0.5), _lv("1_200_C_T", "C", "T", 0.1)],
        [_lv("1_200_C_T", "C", "T", 0.05)], {}, "1_999_T_C",
    )
    assert [p.variant_id for p in pairs] == ["1_200_C_T"]


def test_harmonise_regions_attaches_r2_for_non_lead_variants():
    pairs = harmonise_regions_for_locuscompare(
        [_lv("1_400_C_T", "C", "T", 0.2)], [_lv("1_400_C_T", "C", "T", 0.1)],
        r2_by_variant={"1_400_C_T": 0.42}, lead_variant_id="1_500_A_G",
    )
    assert pairs[0].r2_with_lead == pytest.approx(0.42)


# ----- colour bins, p-value guard, downsampling


def test_r2_bins_are_the_locuszoom_five_and_none_is_grey():
    labels = [label for *_, label in LD_R2_BINS]
    assert labels == ["0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "0.0 - 0.2"]
    assert _r2_color(None)[1] == "0.0 - 0.2"
    assert _r2_color(0.85)[1] == "0.8 - 1.0"
    assert _r2_color(1.0)[1] == "0.8 - 1.0"
    assert _r2_color(0.6)[1] == "0.6 - 0.8"


def test_p_of_zero_yields_no_point_rather_than_a_floor():
    """p <= 0, None or non-finite gives None: the variant is skipped in every panel,
    it is not plotted at a floored value."""
    assert _safe_neg_log10(0.0) is None
    assert _safe_neg_log10(None) is None
    assert _safe_neg_log10(float("nan")) is None
    assert _safe_neg_log10(1e-10) == pytest.approx(10.0)


def test_stratified_downsample_keeps_the_top_r2_bin_and_caps_the_rest():
    pairs = [_pair(f"1_{i}_A_G", i, "A", "G", 0.1, 0.1, 1e-3, 1e-3, 0.9) for i in range(1, 101)]
    pairs += [_pair(f"1_{i}_A_G", i, "A", "G", 0.1, 0.1, 1e-3, 1e-3, 0.1) for i in range(101, 1101)]
    kept = _stratified_downsample(pairs, target_max=200)
    assert len(kept) <= 200
    assert sum(1 for p in kept if p.r2_with_lead == 0.9) == 100  # top bin kept whole
    assert len(_stratified_downsample(pairs[:150], target_max=200)) == 150  # under the cap: untouched


# ----- render_full_locuscompare


def test_render_full_locuscompare_with_gene_track(tmp_path):
    out_path = tmp_path / "smoke_with_gene_track.png"
    assert render_full_locuscompare(_stub_input(with_gene_track=True), out_path) == out_path
    assert out_path.stat().st_size > 10_000
    width, height = _png_pixel_dims(out_path)
    assert width > 0 and height > 0


def test_render_full_locuscompare_without_gene_track_is_shorter(tmp_path):
    """An empty gene track collapses that row; the other panels still render."""
    with_track = tmp_path / "with.png"
    without = tmp_path / "without.png"
    render_full_locuscompare(_stub_input(with_gene_track=True), with_track)
    render_full_locuscompare(_stub_input(with_gene_track=False), without)
    assert without.stat().st_size > 10_000
    assert _png_pixel_dims(without)[1] < _png_pixel_dims(with_track)[1]


def test_render_full_locuscompare_creates_parent_directory(tmp_path):
    out_path = tmp_path / "nested" / "sub" / "smoke.png"
    render_full_locuscompare(_stub_input(), out_path)
    assert out_path.is_file()


def test_render_full_locuscompare_falls_back_to_joined_pairs_for_manhattan(tmp_path):
    """With no per-side track lists the Manhattans render from the joined pairs."""
    inp = _stub_input(with_gene_track=True)
    inp.exposure_track_variants = []
    inp.outcome_track_variants = []
    out_path = tmp_path / "smoke_joined_only.png"
    render_full_locuscompare(inp, out_path)
    assert out_path.stat().st_size > 10_000


def test_render_default_title_reports_the_half_window(tmp_path, monkeypatch):
    """window_bp is the full width of the fetched region. With no title override
    the suptitle reports half of it, the distance covered on each side of the
    lead: 1,000,000 bp renders as "±500 kb", never "±1000 kb". Captured from the
    Figure at savefig time, since the renderer closes the figure before returning."""
    from matplotlib.figure import Figure

    captured: list[str] = []
    real_savefig = Figure.savefig

    def spy(self, *args, **kwargs):
        captured.append(self.get_suptitle())
        return real_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy)
    inp = _stub_input(lead_id="1_500000_A_G", with_gene_track=False)
    inp.title = None
    inp.window_bp = 1_000_000
    render_full_locuscompare(inp, tmp_path / "half_window.png")
    assert captured == ["Regional LocusCompare: 3 variants joined (±500 kb of 1_500000_A_G)"]
