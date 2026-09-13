# 01: synthetic offline demo

A 200-variant synthetic locus where the exposure and the outcome share one causal variant at
chr1:500000 (A>T). No network, no external data; the run took 0.51 s wall clock on one laptop
run.

## Run

```bash
python cli.py --demo 01_synthetic_demo --output /tmp/locuscompare_demo
# equivalently
python cli.py --input examples/01_synthetic_demo/config.json --output /tmp/locuscompare_demo
```

## Files

- `config.json`: the config, pointing at the four fixtures by relative path.
- `exposure.tsv`, `outcome.tsv`: the two summary-statistics slices in the harmonised TSV
  format (`references/input_schema.md`), 200 rows each.
- `ld_matrix.tsv`: r² between the lead and each of the 199 other variants.
- `genes.tsv`: three synthetic protein-coding genes; `DEMOGENE_B` straddles the lead and is
  the focal gene.
- `generate_synthetic_fixtures.py`: the deterministic generator (seed 20260524). Re-run it
  after changing it; the test suite checks that the shipped TSVs match its output byte for
  byte, so do not hand-edit them.

## What you should see

Both Manhattan tracks peak at chr1:500000, colours fade from orange at the lead to grey at the
window edges, and the two scatters show a diagonal (`-log10 p`) and a positive slope (`β`).
The manifest reports `n_pairs: 200` and `n_palindromic_excluded: 1`: the lead itself is an
A/T SNP, so it is flagged palindromic and does not appear as a diamond in the two scatters
(it is still annotated on both Manhattan tracks). `ld_panel` is `synthetic` and
`plink_version` is `prefetched`.
