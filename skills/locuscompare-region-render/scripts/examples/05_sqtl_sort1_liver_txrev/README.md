# 05: splicing QTL x GWAS (SORT1 transcript usage in GTEx liver x cholesterol-VLDL)

A transcript-usage (`txrev`) exposure: the same SORT1 locus, with the QTL side measuring the
fraction of SORT1 expression carried by one transcript rather than total expression.

## Dataset choice

eQTL Catalogue dataset `QTD000269` (GTEx liver, txrev). Counted on 2026-05-15 over
chr1:108.77-109.77 Mb (the +/-500 kb window), the four GTEx liver splicing-type datasets held
these SORT1 rows in their credible-set files:

| Quantification | Dataset | SORT1 rows in the window |
|---|---|---|
| `exon` | QTD000267 | 2,877 (one exon) |
| `tx` | QTD000268 | 5,754 (two transcripts) |
| `txrev` (this demo) | QTD000269 | 2,877 (one transcript-usage trait, ENST00000483508) |
| `leafcutter` | QTD000270 | 0 (no SORT1 intron cluster in the credible-set file) |

`leafcutter` cannot render SORT1 here because the credible-set file keeps no SORT1 cluster; the
other three work by changing `dataset_id`.

## Credible-set file

For splicing, exon, transcript and transcript-usage datasets the eQTL Catalogue fetcher reads
`<dataset>.cc.tsv.gz`, the credible-set-filtered rows (the strongest molecular trait per
fine-mapped signal), so the exposure track is sparser than a gene-expression run. The
composer writes that as a caveat in the caption and in `render_block.ancestry_caveats`.

## Run

```bash
python cli.py --demo 05_sqtl_sort1_liver_txrev --output runs/sqtl_sort1_liver_txrev/
```

The `manifest.yaml` and `report.md` of a run on 2026-05-24 (`n_pairs: 2648`,
`n_palindromic_excluded: 349`, LD from 1000 Genomes EUR via plink 1.9) live in this
repository's test tree, not in the skill, at
`tests/locuscompare-region-render/fixtures/examples/05_sqtl_sort1_liver_txrev/expected_output/`;
the PNG is not shipped. `tests/locuscompare-region-render/test_example_fixtures.py` checks
them against the current code: the manifest's key layout against a fresh offline run (both of these checks run only when the sibling fetch skills are installed; without them the test skips them), the
report's lead line against the line `cli.py` writes, and the counts above against the
manifest. The window label on line 3 of that `report.md` (`±500 kb`) was re-derived from
the current code, which reports the half-window the fetch covered, rather than copied
from the 2026-05-24 run, whose report printed the full width.

## Reading it

A `txrev` beta is the per-allele change in the fraction of expression carried by the
transcript, not in abundance; it is not comparable to the gene-expression beta of demo 02. A
diagonal here says the same lead alters transcript usage as well as total expression.
