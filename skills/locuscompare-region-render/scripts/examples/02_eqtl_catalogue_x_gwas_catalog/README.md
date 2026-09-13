# 02: eQTL Catalogue x GWAS Catalog (the default live demo)

Both sides come from the bundled fetchers; you supply the dataset id, the accession and the
lead.

## The biology

SORT1 x cholesterol in medium VLDL is the 1p13.3 LDL locus of Musunuru 2010 (Nature
2010;466:714-719, doi:10.1038/nature09266): the minor allele at rs12740374 creates a C/EBP
binding site, raises hepatic SORT1 expression and lowers plasma LDL. The exposure here is the
GTEx minor salivary gland gene-expression dataset (eQTL Catalogue QTD000276), an Open Targets
colocalisation row for this pair; liver would be the textbook tissue.

## Run

```bash
python cli.py --demo --output runs/sort1_vldl/
```

Network, all anonymous: the eQTL Catalogue FTP (QTD000276 slice), the GWAS Catalog FTP
(GCST90269602 slice), the 1000 Genomes FTP (chr1 region VCF, needs plink 1.9 locally) and
Ensembl REST (gene track).

## What you should see

- Top: the cholesterol-VLDL Manhattan with a sharp peak at chr1:109274968.
- Middle: the SORT1 eQTL Manhattan peaking at the same position.
- Gene track: CELSR2, PSRC1, SORT1 (bold red), SARS1 and neighbours.
- LocusCompare panel: a diagonal with the lead at the top right.
- Effect-size panel: a negative slope (the allele that raises SORT1 expression lowers
  VLDL cholesterol).

Replaying the recorded 2026-05-15 slices offline (the golden fixture in the test suite) gives
`n_pairs: 2547` and `n_palindromic_excluded: 333`; a live run at a later catalogue release may
differ.

## Attribution

eQTL Catalogue (CC-BY 4.0): Kerimov 2021. GWAS Catalog: Sollis 2023. 1000 Genomes: Auton
2015. GENCODE: Frankish 2021. Full citations are in the skill's SKILL.md.
