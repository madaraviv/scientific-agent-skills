# Examples

Seven numbered demos (runnable with `--demo <NAME>` or `--demo <NN>`), one chain pattern, and
four harmonisation recipes for sources without a bundled fetcher. Run from the skill's
`scripts/` directory.

## Numbered demos

| Demo | Entry vector | Network | Use when |
|---|---|---|---|
| `01_synthetic_demo/` | pre-fetched TSVs + synthetic LD and gene track | no | smoke-testing the install; CI; offline use |
| `02_eqtl_catalogue_x_gwas_catalog/` | bundled fetchers | yes | you have an eQTL Catalogue dataset id and a GWAS Catalog accession (the default for bare `--demo`) |
| `03_open_targets_followup/` | Open Targets colocalisation row | yes | you found a row in an Open Targets colocalisation query |
| `04_gwas_lookup_followup/` | an rsID lookup | yes | you found a variant through a per-variant lookup and want the regional view |
| `05_sqtl_sort1_liver_txrev/` | bundled fetchers (credible-set file) | yes | a splicing / transcript-usage exposure; SORT1 in GTEx liver |
| `06_sceqtl_sort1_onek1k_cd14_mono/` | bundled fetchers | yes | a single-cell eQTL exposure; SORT1 in OneK1K CD14+ monocytes |
| `07_pqtl_sort1_ukbppp_eur/` | bundled fetchers (`source: ukb_ppp`) | yes | a plasma pQTL exposure; needs the optional `ukb-ppp-region-fetch` skill, which is not in this collection |

Demos 02 and 04 to 07 share the outcome GCST90269602 (cholesterol in medium VLDL) and the
SORT1 lead `1_109274968_G_T` (rs12740374); 03 is the LDLR locus. Live demos need plink 1.9
for LD colours; without it the points are grey and the manifest says so.

## Chain pattern (not listed by `--list-demos`)

`chains/finemapping_chain/` shows how to point `sumstats_path:` at the harmonised slices a
fine-mapping run already produced for both sides, so nothing is re-fetched. It has no
self-contained data.

## Recipes for sources without a bundled fetcher

`recipes/` ships bash + awk + bgzip + tabix scripts that convert a source's native file into
the TSV described in the skill's `references/input_schema.md`. Point `sumstats_path:` at the
output.

| Recipe | Source | Access |
|---|---|---|
| `finngen_direct/` | FinnGen phenotype summary statistics | open (check the release's data-use terms) |
| `pan_ukbb_direct/` | Pan-UKBB per-ancestry summary statistics | open |
| `ukb_ppp_pqtl/` | UKB-PPP plasma pQTL | registration-gated (UK Biobank + Synapse) |
| `gtex_v10_direct/` | GTEx v10 cis-eQTL all-pairs files | open |
