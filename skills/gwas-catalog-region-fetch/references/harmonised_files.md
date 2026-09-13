# GWAS Catalog harmonised files as this skill reads them

## Where a file lives

`https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<bucket>/<GCST>/harmonised/`

The bucket is a range of 1,000 accessions: `GCST90269602` sits in
`GCST90269001-GCST90270000/`. `gcst_url_base()` in the script computes it; the accession must
match `^GCST\d+$` or the function raises `ValueError` before any request is made.

A harmonised directory listed on 2026-09-13 (`GCST90269602`) held:

| File | Role |
|---|---|
| `GCST90269602.h.tsv.gz` | The harmonised summary statistics, bgzip-compressed, sorted |
| `GCST90269602.h.tsv.gz.tbi` | Tabix index; what makes a range read possible |
| `GCST90269602.h.tsv.gz-meta.yaml` | Sidecar with trait, sample, build and harmonisation status |
| `GCST90269602.running.log` | The harmoniser's run log |
| `md5sum.txt` | Checksums |

The raw deposited file sits one directory up. The script never opens it.

## File naming

EBI's root-level `harmonised_list.txt` listed 156,116 harmonised files on 2026-09-13
(155,551 distinct accessions):

| Pattern | Count | Example |
|---|---|---|
| `<GCST>.h.tsv.gz` | 136,441 | `GCST90269602.h.tsv.gz` |
| `<PMID>-<GCST>-<EFO>.h.tsv.gz` | 19,635 | `34927100-GCST90019016-EFO_0000676.h.tsv.gz` |
| other | 40 | `GCST90267286_buildGRCh37.h.tsv.gz` (13 of the 15 inspected), one dot-prefixed name, one `.harmonised/` directory |

`harmonised_file_url()` HEADs the simple name first; on any status other than 200, or a request
exception, it GETs the directory listing and returns the first `href` ending in `.h.tsv.gz`. If
neither yields a name it raises `GWASCatalogFetchError`. A `requests.Session` can be injected
(`session=`) to make this step deterministic in tests.

## Two column layouts

Sample of 17 files read on 2026-09-13: the three demo studies, GCST90019016, and every
13,000th line of `harmonised_list.txt`. The column names were read the way the script reads
them (a 32 KB HTTP `Range` read of the file head, gunzipped to the first line).

| Accession | File name | Layout | `.tbi` |
|---|---|---|---|
| GCST90269602 | simple | GWAS-SSF | 200 |
| GCST90691573 | simple | GWAS-SSF | 200 |
| GCST90691576 | simple | GWAS-SSF | 200 |
| GCST90019016 | prefixed | `hm_` columns | 404 |
| GCST90615056 | simple | GWAS-SSF | 200 |
| GCST90242247 | simple | `hm_` columns | 404 |
| GCST90288376 | simple | GWAS-SSF | 200 |
| GCST90243213 | simple | `hm_` columns | 404 |
| GCST90002381 | simple | GWAS-SSF | 200 |
| GCST90399910 | simple | GWAS-SSF (no `beta`; has `z_score`) | 200 |
| GCST90232839 | simple | GWAS-SSF | 200 |
| GCST90801823 | simple | GWAS-SSF | 200 |
| GCST90571827 | simple | GWAS-SSF | 200 |
| GCST90803464 | simple | GWAS-SSF | 200 |
| GCST90385960 | simple | GWAS-SSF | 200 |
| GCST90266254 | simple | GWAS-SSF | 200 |
| GCST90077735 | prefixed | `hm_` columns | 404 |

**`hm_` layout** (4 of 17): `hm_variant_id`, `hm_rsid`, `hm_chrom`, `hm_pos`,
`hm_other_allele`, `hm_effect_allele`, `hm_beta`, `hm_odds_ratio`, `hm_ci_lower`,
`hm_ci_upper`, `hm_effect_allele_frequency`, `hm_code`, followed by the depositor's original
columns (`chromosome`, `base_pair_location`, `effect_allele`, `standard_error`, `p_value`, ...).

**GWAS-SSF layout** (13 of 17): `chromosome`, `base_pair_location`, `effect_allele`,
`other_allele`, `beta`, `standard_error`, `effect_allele_frequency`, `p_value`, `rsid`, then
study-specific extras, then `hm_coordinate_conversion` and `hm_code`. The harmonised values
are written into the standard columns; the sidecar records `genome_assembly: GRCh38`,
`coordinate_system: 1-based`, `is_harmonised: true`.

`_normalise_row()` maps a row as follows, first name wins:

| `RegionVariant` field | Columns tried |
|---|---|
| `chromosome` | `hm_chrom`, `chromosome` (a `chr` prefix is stripped) |
| `position` | `hm_pos`, `base_pair_location` |
| `ref` | `hm_other_allele`, `other_allele` (upper-cased) |
| `alt` | `hm_effect_allele`, `effect_allele` (upper-cased) |
| `variant_id` | `hm_variant_id`, else `<chromosome>_<position>_<ref>_<alt>` |
| `beta` | `hm_beta`, `beta` |
| `se` | `standard_error`, `se` |
| `p_value` | `p_value`, `pvalue` |
| `odds_ratio` | `hm_odds_ratio`, `odds_ratio` |
| `effect_allele_frequency` | `hm_effect_allele_frequency`, `effect_allele_frequency` |

A row missing any of chromosome, position, ref or alt is dropped without a note. `NA`, empty
and non-numeric values become `None`; `NaN` also becomes `None`. A row whose field count differs
from the header's is skipped with a note in `RegionResult.notes`. The whole raw row is kept in
`RegionVariant.raw`.

In the 17-file sample every file without an index was in the `hm_` layout and every file in
the GWAS-SSF layout had one. That is an observation about 17 files, not a rule to code against:
the script does not look at the layout to decide anything.

## The index gap

GCST90019016 on 2026-09-12: `34927100-GCST90019016-EFO_0000676.h.tsv.gz` served (HTTP 200,
261 MB), no `.tbi` in the directory listing, and 2 of 4 sampled accessions with a harmonised
file had no index. Its sidecar reads `is_sorted: false`; tabix requires a sorted bgzip file.
On such a file `pysam.TabixFile(url)` fails, htslib prints
`[E::hts_idx_load3] Could not load local index file '<url>.tbi' : No such file or directory`,
and the script raises `GWASCatalogFetchError("could not open tabix index for <url>: ...")`.
There is no fallback in this script.

## The `-meta.yaml` sidecar

Fields seen in the sidecars of GCST90269602 and GCST90019016 (2026-09-13): `gwas_id`,
`gwas_catalog_api`, `date_metadata_last_modified`, `trait_description` (list),
`ontology_mapping` (list of EFO ids, may be empty), `genome_assembly`, `coordinate_system`,
`genotyping_technology` (list), `samples` (list of `sample_ancestry_category` list plus
`sample_size`), `sex`, `data_file_name`, `file_type`, `data_file_md5sum`,
`minor_allele_freq_lower_limit`, `is_harmonised`, `is_sorted`, `harmonisation_reference`, and
sometimes a free-text `author_notes`. The script does not read it; the fields above are what a
caller can get with one small GET when the REST API is not wanted.

## What the harmoniser does to a row

From the pipeline's documentation
([Introduction](https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/),
read 2026-09-13):

1. **Genome build mapping.** Position updated to GRCh38 by rsID lookup against Ensembl
   (release 95, dbSNP 151); `hm_coordinate_conversion = rs`. When no rsID maps, UCSC liftOver;
   `hm_coordinate_conversion = lo`. When neither works the row is removed. Liftover without an
   rsID can map two distinct input variants onto one GRCh38 position; the documentation's
   Limitations page shows such a duplicate.
2. **Palindromic variants (A/T, G/C).** Strand is inferred from a strand-consensus rate
   computed on a 10% sample of non-palindromic sites (forward / (forward + reverse)): at or
   above 0.995 the inferred strand is used; between 0.9 and 0.995 the rate is recomputed on all
   non-palindromic sites and used if above 0.99, else palindromes are dropped; at or below 0.9
   they are dropped.
3. **Harmonising to the reference.** Each row is matched to the Ensembl VCF by chromosome,
   position and alleles. Reverse-strand alleles are reverse-complemented. When the effect and
   other alleles are the reverse of the reference pair they are swapped, with
   `beta = -beta`, `odds_ratio = 1 / odds_ratio`, the confidence interval bounds exchanged and
   inverted, and `effect_allele_frequency = 1 - effect_allele_frequency`.
4. **Quality control.** Rows lacking a variant id, chromosome, position or p-value are removed.
5. **Output.** Bgzip-compressed, sorted, tabix-indexed `.h.tsv.gz` with sidecar and log.

## `hm_code`

From the pipeline's reference guide (read 2026-09-13). Rows with codes 9 and 14 to 18 are not in
the harmonised file.

| Code | Meaning |
|---|---|
| 1 | Palindromic; strand inferred; forward strand; alleles correct |
| 2 | Palindromic; strand inferred; forward strand; alleles flipped |
| 3 | Palindromic; strand inferred; reverse strand; alleles correct |
| 4 | Palindromic; strand inferred; reverse strand; alleles flipped |
| 5 | Palindromic; forward strand assumed; alleles correct |
| 6 | Palindromic; forward strand assumed; alleles flipped |
| 7 | Palindromic; reverse strand assumed; alleles correct |
| 8 | Palindromic; reverse strand assumed; alleles flipped |
| 9 | Palindromic; dropped; not harmonised |
| 10 | Forward strand; alleles correct |
| 11 | Forward strand; alleles flipped |
| 12 | Reverse strand; alleles correct |
| 13 | Reverse strand; alleles flipped |
| 14 | Required fields not known; not harmonised |
| 15 | No matching variant in the reference VCF; not harmonised |
| 16 | Multiple matching variants in the reference VCF; not harmonised |
| 17 | Palindromic; strand inferred; EAF or reference AF not known; not harmonised |
| 18 | Palindromic; strand inferred; EAF below the minor allele frequency threshold; not harmonised |

`hm_coordinate_conversion`: `rs` (position from rsID lookup) or `lo` (position from liftOver).
