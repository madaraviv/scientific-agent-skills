"""Unit tests for `_wald_ratio_types.py`, the shared dataclasses the renderer and
harmoniser use: `LocusVariant.from_ot_row`, `harmonise_locus_intersection` and
`wald_ratio_at_lead` (the slope drawn on the effect-size panel)."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

from _wald_ratio_types import (  # noqa: E402
    HarmonisedPair,
    LocusVariant,
    PALINDROMIC_PAIRS,
    WrCall,
    harmonise_locus_intersection,
    wald_ratio_at_lead,
)


def test_palindromic_pairs_cover_at_and_gc():
    assert {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")} <= set(PALINDROMIC_PAIRS)
    assert ("A", "G") not in PALINDROMIC_PAIRS
    assert ("C", "T") not in PALINDROMIC_PAIRS


def test_locus_variant_from_ot_row_combines_mantissa_and_exponent():
    row = {
        "variant": {"id": "1_109274968_G_T", "chromosome": "1", "position": 109_274_968,
                    "referenceAllele": "G", "alternateAllele": "T"},
        "posteriorProbability": 0.42, "beta": 0.31, "standardError": 0.05,
        "pValueMantissa": 3.0, "pValueExponent": -8,
        "is95CredibleSet": True, "is99CredibleSet": True, "r2Overall": 0.92,
    }
    lv = LocusVariant.from_ot_row(row)
    assert (lv.variant_id, lv.chromosome, lv.position) == ("1_109274968_G_T", "1", 109_274_968)
    assert (lv.ref, lv.alt) == ("G", "T")
    assert lv.pip == 0.42 and lv.beta == 0.31 and lv.se == 0.05
    assert lv.p_value == pytest.approx(3.0e-8)
    assert lv.is95 is True and lv.is99 is True and lv.r2_lead == 0.92


def test_locus_variant_from_ot_row_handles_missing_pvalue_components():
    row = {"variant": {"id": "1_1_A_G", "referenceAllele": "A", "alternateAllele": "G"},
           "posteriorProbability": 0.1, "beta": 0.0, "standardError": 0.1}
    assert LocusVariant.from_ot_row(row).p_value is None


def _lv(vid, ref, alt, beta, *, pip=0.5, se=0.05, p_value=1e-8):
    return LocusVariant(variant_id=vid, chromosome="1", position=int(vid.split("_")[1]),
                        ref=ref, alt=alt, pip=pip, beta=beta, se=se, p_value=p_value,
                        is95=True, is99=True, r2_lead=1.0)


def test_harmonise_intersection_drops_variants_only_on_one_side():
    pairs = harmonise_locus_intersection(
        [_lv("1_100_A_G", "A", "G", 0.4), _lv("1_200_C_T", "C", "T", 0.2)],
        [_lv("1_100_A_G", "A", "G", 0.5)],
    )
    assert [p.variant_id for p in pairs] == ["1_100_A_G"]
    assert pairs[0].sign_flipped is False and pairs[0].palindromic is False
    assert pairs[0].pip_product == pytest.approx(0.25)


def test_harmonise_intersection_flips_swapped_alleles():
    pairs = harmonise_locus_intersection([_lv("1_100_A_G", "A", "G", 0.4)],
                                         [_lv("1_100_A_G", "G", "A", 0.6)])
    assert len(pairs) == 1
    assert pairs[0].sign_flipped is True
    assert pairs[0].beta_left == pytest.approx(0.4)
    assert pairs[0].beta_right == pytest.approx(-0.6)


def test_harmonise_intersection_keeps_palindromic_but_flags_them():
    pairs = harmonise_locus_intersection([_lv("1_100_A_T", "A", "T", 0.5)],
                                         [_lv("1_100_A_T", "A", "T", 0.6)])
    assert len(pairs) == 1 and pairs[0].palindromic is True


def test_harmonise_intersection_drops_irreconcilable_alleles():
    assert harmonise_locus_intersection([_lv("1_100_A_G", "A", "G", 0.4)],
                                        [_lv("1_100_A_G", "C", "T", 0.6)]) == []


def test_harmonise_intersection_drops_missing_beta():
    right_missing = LocusVariant(variant_id="1_100_A_G", chromosome="1", position=100,
                                 ref="A", alt="G", pip=0.5, beta=None, se=None, p_value=None,
                                 is95=None, is99=None, r2_lead=None)
    assert harmonise_locus_intersection([_lv("1_100_A_G", "A", "G", 0.4)], [right_missing]) == []


def _pair(vid, *, beta_left, beta_right, pip_left=0.5, pip_right=0.5,
          se_left=0.05, se_right=0.05, palindromic=False):
    return HarmonisedPair(variant_id=vid, pip_left=pip_left, pip_right=pip_right,
                          beta_left=beta_left, beta_right=beta_right,
                          se_left=se_left, se_right=se_right, p_left=1e-8, p_right=1e-8,
                          sign_flipped=False, palindromic=palindromic)


def test_wald_ratio_picks_pip_product_argmax():
    res = wald_ratio_at_lead([
        _pair("v1", beta_left=0.5, beta_right=0.2, pip_left=0.10, pip_right=0.10),
        _pair("v2", beta_left=0.5, beta_right=0.2, pip_left=0.80, pip_right=0.80),
        _pair("v3", beta_left=0.5, beta_right=0.2, pip_left=0.40, pip_right=0.40),
    ])
    assert res.lead is not None and res.lead.variant_id == "v2"
    assert res.wr == pytest.approx(0.4)
    assert res.call is WrCall.ALIGNED
    assert (res.n_intersected, res.n_eligible, res.n_palindromic_excluded) == (3, 3, 0)


def test_wald_ratio_opposite_sign_call_when_betas_disagree():
    res = wald_ratio_at_lead([_pair("v1", beta_left=0.5, beta_right=-0.3, pip_left=0.9, pip_right=0.9)])
    assert res.wr == pytest.approx(-0.6)
    assert res.call is WrCall.OPPOSITE
    assert res.se_wr is not None and res.se_wr > 0
    assert res.ci_lo < res.wr < res.ci_hi


def test_wald_ratio_excludes_palindromic_from_lead_selection():
    res = wald_ratio_at_lead([
        _pair("vpal", beta_left=10.0, beta_right=10.0, pip_left=0.99, pip_right=0.99, palindromic=True),
        _pair("vok", beta_left=0.5, beta_right=0.2, pip_left=0.4, pip_right=0.4),
    ])
    assert res.lead is not None and res.lead.variant_id == "vok"
    assert res.n_palindromic_excluded == 1
    assert any("palindromic" in n for n in res.notes)


def test_wald_ratio_undetermined_when_no_eligible_pair():
    res = wald_ratio_at_lead([_pair("vpal", beta_left=0.5, beta_right=0.2,
                                    pip_left=0.9, pip_right=0.9, palindromic=True)])
    assert res.lead is None and res.call is WrCall.UNDETERMINED and res.wr is None


def test_wald_ratio_undefined_when_exposure_beta_zero():
    res = wald_ratio_at_lead([_pair("v1", beta_left=0.0, beta_right=0.5, pip_left=0.9, pip_right=0.9)])
    assert res.call is WrCall.UNDETERMINED and res.wr is None
    assert any("zero" in n.lower() for n in res.notes)


def test_wald_ratio_handles_missing_right_se():
    res = wald_ratio_at_lead([_pair("v1", beta_left=0.5, beta_right=0.2, se_right=None,
                                    pip_left=0.9, pip_right=0.9)])
    assert res.wr == pytest.approx(0.4)
    assert res.se_wr is None or math.isfinite(res.se_wr)


def test_wald_ratio_empty_input_returns_undetermined():
    res = wald_ratio_at_lead([])
    assert res.lead is None and res.call is WrCall.UNDETERMINED
    assert (res.n_intersected, res.n_eligible) == (0, 0)
