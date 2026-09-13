# Harmonised summary-statistics TSV for `sumstats_path` inputs

The TSV contract for a side supplied by `sumstats_path:` instead of `fetch:`. The loader is
`load_sumstats_tsv` in `scripts/_prefetched.py`; the four recipes under
`scripts/examples/recipes/` write this format from FinnGen, Pan-UKBB, UKB-PPP and GTEx v10
native files.

## File format

- UTF-8, Unix line endings, literal TAB separator.
- Plain `.tsv`, or gzip / bgzip compressed (`.tsv.gz`, `.tsv.bgz`, `.gz`, `.bgz`); the loader
  opens by extension.
- Lines starting with `#` are skipped (the recipes write `# locuscompare-schema-version: 1.0`
  and `# source: ...` there), so are blank lines. The first remaining line is the header.
- Sort order is not required by the loader; sort by chromosome and position if you also want
  to tabix-index the file for other tools.

## Required columns

| Column | Type | Meaning |
|---|---|---|
| `variant_id` | string | `chrom_pos_ref_alt`, GRCh38, no `chr` prefix (`1_109274968_G_T`). The join key between the two sides and the id the LD client is asked about. |
| `chromosome` | string | Without `chr`: `1` .. `22`, `X`, `Y`, `MT`. |
| `position_bp` | integer | 1-based GRCh38 position. |
| `allele_a` | string | The non-effect allele. |
| `allele_b` | string | The effect allele; `beta` is per copy of `allele_b`. |
| `beta` | float | Effect estimate. |
| `se` | float | Standard error of `beta`. |
| `p` | float | Two-sided p-value. |

## Optional columns

`maf` (minor-allele frequency) and `eaf` (frequency of `allele_b`) are read into the variant
record when present. Any other column (`n`, `rsid`, `info`, `molecular_trait_id`, `study_id`,
GRCh37 coordinates) is ignored by the loader; keep them for your own provenance.

## What the loader checks, and what it does not

Checks (each failure raises `PrefetchedSchemaError` with the file and line):

1. All eight required columns are in the header.
2. Every data row has at least as many fields as the header.
3. `position_bp` parses as an integer; `beta`, `se`, `p`, `maf`, `eaf` parse as floats, with
   an empty field or `NA` meaning missing.

Behaviour that is not a check:

- A row with a missing `beta`, `se` or `p` is dropped silently before harmonisation.
- No range check: `p = 2.0`, `se = -0.05` or `maf = 0.9` are loaded as given. `p <= 0` later
  yields no point in any panel (see the layout reference).
- No consistency check between `variant_id` and the `chromosome` / `position_bp` /
  `allele_a` / `allele_b` columns. The join uses `variant_id`; the Manhattan x-position uses
  `position_bp`; allele reconciliation uses `allele_a` / `allele_b`. If they disagree the
  figure is quietly wrong, so build the id from the four columns as the recipes do.
- No duplicate-id check. With duplicate `variant_id` rows on the outcome side the last row wins
  in the join; on the exposure side each duplicate row joins separately.

## Conventions the renderer assumes

- **Coordinates:** GRCh38 throughout. Lift over GRCh37 sources first and drop variants that do
  not map uniquely.
- **Effect allele:** `allele_b`. This matches the eQTL Catalogue (`alt` is the effect allele),
  the GWAS Catalog harmonised files (`hm_effect_allele`) and the Open Targets locus tables. A
  source whose effect allele is the reference must have `beta` negated and the two allele
  columns swapped before it is written in this format.
- **Palindromic SNPs (A/T, G/C):** write them; the composer flags them, counts them in
  `n_palindromic_excluded`, keeps them in the Manhattan tracks and leaves them out of the
  two scatters. Do not pre-filter.
- **Multi-allelic sites:** one row per (chrom, pos, ref, alt); never several alts in one row.
- **Indels:** left-aligned, `bcftools norm` style, with the longest common prefix removed. The
  reconciliation compares allele strings exactly, so both sides must use the same
  representation.

## Companion files

The two synthetic inputs of the offline demo use the same conventions:

- `ld.ld_matrix_path`: TSV with `partner_variant_id` and `r2` (0 to 1; the lead's own row is
  not needed, it is seeded to 1.0). Partners absent from the file render grey.
- `gene_track.genes_path`: TSV with `gene_symbol`, `start`, `end`, `strand` and an optional
  `biotype` (default `protein_coding`). No exons; genes draw as bars.

## Recipes

Each recipe is a bash + awk pipeline that writes the TSV, sorts it, bgzips it and builds a
tabix index (`bgzip` and `tabix` from htslib on PATH):

| Recipe | Source | Output columns beyond the required eight |
|---|---|---|
| `scripts/examples/recipes/finngen_direct/harmonise.sh` | FinnGen phenotype TSV (GRCh38) | `eaf`, `rsid` |
| `scripts/examples/recipes/pan_ukbb_direct/harmonise.sh` | Pan-UKBB phenotype TSV, one ancestry's columns | `eaf` |
| `scripts/examples/recipes/ukb_ppp_pqtl/harmonise.sh` | UKB-PPP per-protein TSV (registration-gated) | `n` |
| `scripts/examples/recipes/gtex_v10_direct/harmonise.sh` | GTEx v10 all-pairs file, filtered to one gene | `molecular_trait_id` |

Each recipe's README documents the source URL, the column map and the effect-allele
orientation of that source. Point `sumstats_path` at the `.gz` the recipe writes.
