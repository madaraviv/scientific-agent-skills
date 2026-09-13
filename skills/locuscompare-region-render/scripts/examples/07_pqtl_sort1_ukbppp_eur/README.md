# 07: plasma pQTL x GWAS (SORT1 plasma sortilin, UKB-PPP European x cholesterol-VLDL)

The protein-side companion of demo 02: the same lead and outcome, with plasma sortilin
abundance (UKB-PPP, Olink) as the exposure. `exposure.fetch.source: ukb_ppp` dispatches to the
`ukb-ppp-region-fetch` skill.

## Requirement

`ukb-ppp-region-fetch` is optional. When it is not installed beside this skill under
`skills/`, this demo is refused with exit code 2 and the message
`exposure.fetch.source=ukb_ppp needs the sibling skill ukb-ppp-region-fetch ...`, and a pQTL
exposure through the Python API raises `Tier2NotAvailable`. The config and this note document
the pQTL entry vector either way.

## Run (with the sibling installed)

```bash
python cli.py --demo 07_pqtl_sort1_ukbppp_eur --output runs/pqtl_sort1_ukbppp/
```

## What the recorded run showed

The `manifest.yaml` and `report.md` of a run on 2026-05-24 with the sibling installed
(`n_pairs: 4606`, `n_palindromic_excluded: 635`) live in this repository's test tree, not
in the skill, at
`tests/locuscompare-region-render/fixtures/examples/07_pqtl_sort1_ukbppp_eur/expected_output/`;
the PNG is not shipped. `tests/locuscompare-region-render/test_example_fixtures.py` checks
them against the current code: the manifest's key layout against a fresh offline run, the
report's lead line against the line `cli.py` writes, and the counts above against the
manifest. The window label on line 3 of that `report.md` (`±500 kb`) was re-derived from
the current code, which reports the half-window the fetch covered, rather than copied
from the 2026-05-24 run, whose report printed the full width.
The exposure is UKB-PPP release 1 (Sun 2023, Nature 2023;622:329-338,
doi:10.1038/s41586-023-06592-6), SORT1 (UniProt Q99523, Olink OID20213), European discovery
cohort. At rs12740374 the T allele raises plasma sortilin (beta +0.12 in that release) while
the same allele lowers hepatic SORT1 mRNA; the mRNA-versus-protein discordance is a property of
the locus, not an artefact of the figure.

## Alternative without the sibling

Harmonise a UKB-PPP per-protein file with `recipes/ukb_ppp_pqtl/harmonise.sh` (registration
required) and supply it through `exposure.sumstats_path`; that path needs only the three
required siblings.
