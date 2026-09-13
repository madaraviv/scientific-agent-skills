# 06: single-cell eQTL x GWAS (SORT1 in OneK1K CD14+ monocytes x cholesterol-VLDL)

The eQTL Catalogue harmonises single-cell cohorts (OneK1K, Perez 2022, Nathan 2022, Randolph
2021) into the same gene-expression schema as bulk tissue, so a single-cell dataset runs
through the same path as demo 02 with a different `dataset_id`.

## Dataset

eQTL Catalogue `QTD000609` (OneK1K, CD14+ monocytes, gene expression), one of the OneK1K
PBMC cell-type datasets; monocytes were chosen because SORT1 is expressed there. Hepatocytes,
the textbook SORT1 tissue, are not a OneK1K cell type.

## Run

```bash
python cli.py --demo 06_sceqtl_sort1_onek1k_cd14_mono --output runs/sceqtl_sort1_onek1k/
```

A run on 2026-05-24 with plink 1.9 available reported `n_pairs: 2648` and
`n_palindromic_excluded: 349` in its manifest.

## Caveat

Per-cell-type sample sizes are smaller than bulk-tissue ones, so the exposure signal can be
weak or the number of joined variants small; read `n_pairs` in the report before reading the
scatters. The config's `caveats:` entry is printed in the caption for that reason.
