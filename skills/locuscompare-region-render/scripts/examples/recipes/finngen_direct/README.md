# Recipe: FinnGen direct download

Convert a FinnGen R12+ phenotype TSV into the canonical locuscompare schema.

## Source

FinnGen direct download: <https://finngen.gitbook.io/documentation/data-download>

- **Format**: tab-separated values, bgzip + tabix
- **Coordinate system**: GRCh38 (since R6 onwards)
- **License**: open-access for sumstats since FinnGen R12; check the current
  release's data-use agreement at <https://www.finngen.fi/en/access_results>
- **Auth**: none for sumstats (was registration-gated until R11; relaxed in R12)
- **Why use direct over GWAS Catalog**: more recent releases, bigger N, and
  FinnGen can carry Finnish-enriched signals. In the FinnGen release 5 analysis (data freeze 5, 224,737 participants)
  (Kurki 2023, cited below), fine-mapping implicated 148 coding variants; 91 of
  those had an allele frequency below 5% in non-Finnish Europeans, and 62 of the
  91 were enriched more than twofold in Finland.

## Quick start

```bash
# 1. Download a phenotype (e.g. heart failure I9_HEARTFAIL)
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_I9_HEARTFAIL.gz

# 2. Harmonise to canonical schema
bash examples/recipes/finngen_direct/harmonise.sh \
    finngen_R12_I9_HEARTFAIL.gz \
    finngen_R12_I9_HEARTFAIL.canonical.tsv

# 3. Drop in to a locuscompare config
# outcome:
#   trait_label: "heart failure (FinnGen R12 I9_HEARTFAIL)"
#   sumstats_path: "finngen_R12_I9_HEARTFAIL.canonical.tsv.gz"
```

## Caveats

- **Ancestry**: FinnGen is a Finnish cohort. The 1000 Genomes EUR super-population
  is the closest bundled LD reference but not a Finnish panel, so r² colours are an
  approximation. Say so in the config's `caveats:` block; it is printed in the
  figure's caption.
- **`af_alt_cases` vs `af_alt_controls`**: this recipe uses `af_alt`
  (population pooled). For balanced trait phenotypes that's fine; for highly
  imbalanced (case fraction < 5%), prefer `af_alt_controls` as the maf
  estimate. Edit the awk mapping accordingly.
- **`rsids`**: FinnGen ships comma-separated rsids when a position is multi-
  mapped in dbSNP. This recipe takes the first; the join key is `variant_id`
  (chr_pos_ref_alt), so this is informational only.

## Sign-check status

FinnGen R12 sumstats are sign-checked against the alt allele as the effect
allele (matches our canonical convention). No flip needed.

## Citation

> Kurki et al. (2023) *FinnGen provides genetic insights from a well-phenotyped
> isolated population.* Nature 613, 508-518. doi:10.1038/s41586-022-05473-8,
> PMID 36653562
