---
name: gwas-catalog-region-fetch
description: Fetch a region of GWAS summary statistics for one study (GCST accession) from the NHGRI-EBI GWAS Catalog harmonised collection with a tabix range read on the EBI FTP. Returns per-variant rows (variant_id, chromosome, position, ref, alt, beta, SE, p-value, effect-allele frequency; GRCh38, effect allele = alt) for colocalisation, fine-mapping, regional plotting or Mendelian randomisation. Use when an agent needs every variant's summary statistics inside one chromosomal window for one study, not a single-variant lookup or genome-wide top hits.
license: MIT
compatibility: Python 3.10+; pysam (htslib) for tabix range reads over HTTPS; requests; network access to the EBI GWAS Catalog FTP (harmonised summary statistics). PyYAML only for YAML configs and manifest.yaml (manifest.json is written without it).
metadata:
  version: "1.0"
  skill-author: Aviv Madar
---

# GWAS Catalog Region Fetch

## Overview

The [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/) (Sollis 2023) redistributes the
summary statistics deposited with it and runs each file through its harmonisation pipeline
([gwas-sumstats-harmoniser](https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/)):
positions on GRCh38, alleles on the forward strand, and the effect allele with its beta, odds
ratio and frequency aligned to the Ensembl reference. The harmonised file for a study lives at
`.../summary_statistics/<bucket>/<GCST>/harmonised/<name>.h.tsv.gz`, bgzip-compressed, with a
sibling `.tbi` tabix index (13 of 17 sampled files had one; Gotcha 2), so one chromosomal window
can be read without downloading the file. EBI's `harmonised_list.txt` named 156,116 harmonised
files on 2026-09-13.

This skill ships one script, `scripts/gwas_catalog_region_fetch.py`.
`GWASCatalogClient.fetch_region(accession, chromosome, start_bp, end_bp)` resolves the harmonised
file name, opens it with pysam's tabix reader over HTTPS, reads the rows in the window and maps
each to a `RegionVariant` (variant id as `chr_pos_ref_alt`, ref = other allele, alt = effect
allele, beta, SE, p-value, odds ratio, effect-allele frequency, and the raw row). The CLI wraps
the client and writes `variants.tsv`, `manifest.yaml` and `report.md`. The demo study's file is
600 MB; its 1 Mb window read returned 5,117 variants in 4.1 s (one measurement, 2026-09-13).

## Trigger

**Fire when** the user (or an upstream agent step) wants:

- The summary statistics (beta, SE, p-value, effect-allele frequency) of every variant in a
  chromosomal window, for one GWAS Catalog study.
- The outcome side of a colocalisation against an eQTL or pQTL signal, a regional association
  plot, a fine-mapping input, or the outcome for Mendelian randomisation at one locus.
- Provenance with the rows: the exact file URL and the fetch time, in a manifest.

Trigger phrases include "GWAS region fetch", "GWAS Catalog region", "GWAS sumstats slice",
"GCST harmonised tabix", "summary statistics around a locus".

**Do NOT fire when** the user wants:

- **One variant in one study**: a point lookup; use a database-query skill such as
  `database-lookup` against the GWAS Catalog REST API.
- **Genome-wide lead associations for a trait**: the GWAS Catalog REST API `/associations`
  endpoint serves curated top hits; this skill reads full per-variant files.
- **A phenome-wide scan of one variant** across many studies: that is many studies, not one
  window of one study.
- **FinnGen, Pan-UKBB, Biobank Japan or UKB-PPP files from their own hosts**: only what the
  GWAS Catalog has harmonised is reachable here.
- **Fine-mapping credible sets or posterior probabilities**: the harmonised files carry none.
- **LD between the variants**: pair this skill with an LD skill such as
  `ld-1000g-region-compute`; this one returns association statistics only.
- **The raw deposited file** (the depositor's build and allele conventions): this skill reads
  `harmonised/` only.

## Quick start

```bash
# Bundled demo: cholesterol in medium VLDL (GCST90269602), 1 Mb around SORT1, chr1
python scripts/gwas_catalog_region_fetch.py --demo --output ./out/sort1_vldl

# Your own config
python scripts/gwas_catalog_region_fetch.py --input my_region.json --output ./out/my_region

# List the bundled configs
python scripts/gwas_catalog_region_fetch.py --list-demos
```

The demo, run on 2026-09-13:

```
info: using bundled demo default.json
gwas-catalog-region-fetch: 5117 variants -> out/sort1_vldl/variants.tsv
  source: GWAS Catalog harmonised | accession GCST90269602
```

From Python (run from `scripts/`, or put it on `sys.path`):

```python
from gwas_catalog_region_fetch import GWASCatalogClient

result = GWASCatalogClient().fetch_region(
    accession="GCST90269602", chromosome="1",
    start_bp=108_774_968, end_bp=109_774_968,
)
print(result.n_variants, result.release.harmonised_url)
for v in result.variants[:3]:
    print(v.variant_id, v.ref, v.alt, v.beta, v.se, v.p_value)
```

`scripts/examples/run_example.sh` does the same for `scripts/examples/input.json`.

`scripts/examples/expected_output.md` predates the current file and is illustrative only; the demo fetch on 2026-09-13 (n = 5,117 rows) returned 5,117 variants, not the ~3000 it states, the lead row `1_109274968_G_T` is rs12740374 with p = 1.5e-78 (rs646776 is `1_109275908_C_T`), and rows carry `effect_allele_frequency`, not `maf`.

Bundled configs (`scripts/examples/`), traits read from each file's `-meta.yaml` sidecar on
2026-09-13:

| Config | Study | Trait | N | Window (GRCh38) |
|---|---|---|---|---|
| `default.json`, `sort1_cholesterol_vldl.json`, `input.json` | GCST90269602 | Cholesterol in medium VLDL | 88,329 | chr1:108,774,968-109,774,968 (SORT1; Musunuru 2010) |
| `il6r_crp.yaml` | GCST90691573 | C-reactive protein levels | 400,094 | chr1:153,925,508-154,925,508 (IL6R) |
| `tcf7l2_hba1c.json` | GCST90691576 | Glycated haemoglobin (HbA1c) levels | 400,825 | chr10:112,531,629-113,531,629 (TCF7L2) |

## Config schema

JSON or YAML (YAML needs PyYAML):

```json
{
  "accession": "GCST90269602",
  "chromosome": "1",
  "start_bp": 108774968,
  "end_bp": 109774968
}
```

| Field | Type | Description | Required |
|---|---|---|---|
| `accession` | string | GWAS Catalog study accession, `GCST` followed by digits only; anything else is rejected before a URL is built. | Yes |
| `chromosome` | string | Chromosome name; a `chr` prefix is stripped, and the other spelling is retried if the index uses it. | Yes |
| `start_bp` | integer | Window start, 1-based GRCh38 inclusive. | Yes |
| `end_bp` | integer | Window end, 1-based GRCh38 inclusive. | Yes |

Other keys (the bundled configs carry a `_description`) are ignored.

## Output files

`<output>/variants.tsv`: three `#` lines (schema version, source, accession), a header, then one
row per variant. `allele_a` is the other allele, `allele_b` the effect allele, and `beta` is per
copy of `allele_b`. Numbers are written with six significant digits; a missing value is an empty
cell. First rows of the demo:

```
# locuscompare-schema-version: 1.0
# source: gwas_catalog
# accession: GCST90269602
variant_id	chromosome	position_bp	allele_a	allele_b	beta	se	p	eaf	study_id
1_108775337_C_T	1	108775337	C	T	-0.00913978	0.00497078	0.0659628	0.339325	GCST90269602
1_108775456_T_A	1	108775456	T	A	0.0580302	0.105754	0.583193	0.00114873	GCST90269602
```

`<output>/manifest.yaml` (or `manifest.json` without PyYAML): `skill`, `version` (the script's
own, `0.1.0`), `accession`, `region` (`chromosome`, `start_bp`, `end_bp`), `n_variants`,
`release` (`accession`, `harmonised_url`, `harmoniser_version`, `fetched_at_utc`) and `outputs`.
`harmoniser_version` is always `null`: the script does not read the `-meta.yaml` sidecar.

`<output>/report.md`: the same facts as a five-line Markdown list.

The Python result carries more than the files: `RegionVariant.odds_ratio`, the full raw row
(`RegionVariant.raw`, every column of the source file including `hm_code`), and
`RegionResult.notes` (one note per row skipped for a column-count mismatch). None of these reach
`variants.tsv`, the manifest or the report; the cache JSON (see Operational notes) has all of them.

## Gotchas

1. **The harmonised files come in two column layouts, and the script reads both.** One layout
   carries the harmonised values in `hm_`-prefixed columns (`hm_chrom`, `hm_pos`,
   `hm_other_allele`, `hm_effect_allele`, `hm_beta`, `hm_odds_ratio`,
   `hm_effect_allele_frequency`, with `standard_error` and `p_value` unprefixed); the other
   writes them into the standard GWAS-SSF columns (`chromosome`, `base_pair_location`,
   `other_allele`, `effect_allele`, `beta`, `standard_error`, `effect_allele_frequency`,
   `p_value`) and prefixes only `hm_coordinate_conversion` and `hm_code`. In a sample of 17 files
   read on 2026-09-13 (the three demo studies, GCST90019016, and every 13,000th line of
   `harmonised_list.txt`), 13 had the second layout and 4 the first. The script tries the
   `hm_` column first and falls back to the standard one, so both layouts yield the same
   `RegionVariant`; in the second layout the standard columns are already GRCh38 and
   forward-strand (the sidecar says `genome_assembly: GRCh38`, `is_harmonised: true`), and the
   variant id is built as `chr_pos_ref_alt` because no `hm_variant_id` exists. Details in
   [references/harmonised_files.md](references/harmonised_files.md).

2. **EBI publishes some harmonised files without a tabix index, and this script then raises its
   generic fetch error.** Measured 2026-09-12 on GCST90019016: the `.h.tsv.gz` is served (HTTP
   200) and its `.tbi` is not (404); 2 of 4 sampled accessions with a harmonised file had no
   index. Re-measured 2026-09-13 over the 17-file sample above: 4 had no index, and those 4 were
   exactly the 4 in the `hm_` layout (a correlation in a sample of 17, not a rule).
   GCST90019016's sidecar says `is_sorted: false`, which is a file tabix cannot index. The call
   fails with `GWASCatalogFetchError: could not open tabix index for <url>: ...` after htslib
   prints `[E::hts_idx_load3] Could not load local index file ...`. The client's docstring
   mentions a streaming fallback; none is implemented: the call raises, and the caller decides
   whether to download the whole file and index it locally.

3. **Two file-naming conventions, resolved by a HEAD then a directory listing.** Of the 156,116
   files in `harmonised_list.txt` on 2026-09-13, 136,441 are named `<GCST>.h.tsv.gz`, 19,635
   `<PMID>-<GCST>-<EFO>.h.tsv.gz`, and 40 something else (of the 15 inspected: 13
   `<GCST>_buildGRCh37.h.tsv.gz`, one dot-prefixed name, one `.harmonised/` directory). The
   script HEADs the simple name and, on anything but 200, lists `harmonised/` and takes the
   first `*.h.tsv.gz` link. That is one extra request for a prefixed file, and in a directory
   with more than one `.h.tsv.gz` the first listed wins silently, so read `harmonised_url` in
   the manifest.

4. **pysam writes the remote index into the current working directory.** After the demo the
   working directory held `GCST90269602.h.tsv.gz.tbi` (1.6 MB), which htslib downloaded to open
   the file; the script does not clean it up, and it is not the skill's cache. Run from a
   scratch directory, or expect one `.tbi` per study.

5. **A chromosome the index does not know raises `ValueError`, not `GWASCatalogFetchError`.**
   The script queries the name without `chr`, retries with `chr` when tabix rejects it, and
   lets the second rejection propagate as pysam's `ValueError`. Catch both when looping over
   chromosomes, and do not read an empty result as "no variants": an absent contig errors and an
   empty window returns `n_variants: 0`.

6. **Beta is per copy of the effect allele (`allele_b`), and not every file has one.** The
   harmoniser flips the sign, inverts the odds ratio and takes `1 - EAF` when it swaps alleles
   to match the reference, but the TSV does not record which allele of the pair is the
   reference's alternate; join other sources on `variant_id` and re-check alleles before
   comparing effects
   ([references/effect_allele_harmonisation.md](references/effect_allele_harmonisation.md)).
   One of the 17 sampled files (GCST90399910) has `z_score` and no `beta` column: its rows land
   with an empty `beta` cell. A binary-trait file's `hm_odds_ratio` is kept on the Python object
   and in the cache JSON only; the TSV has no odds-ratio column.

7. **A variant missing from the slice may have been dropped by the harmoniser, not by the
   study.** Rows the pipeline could not place (`hm_code` 9 and 14 to 18: dropped palindromes, no
   or several reference matches, missing fields) are absent from the file. Absence is not a null
   effect. The `hm_code` of each returned row is in `RegionVariant.raw`, not in the TSV.

8. **The result cache is unversioned.** Every CLI run first looks in
   `~/.clawbio/gwas_catalog_region_fetch_cache/` (override with `GWAS_CATALOG_CACHE_DIR`) for
   `<GCST>__chr<C>_<start>_<end>.json` and returns it without touching the network; a re-harmonised
   file at EBI is never noticed. `--no-cache` skips both the read and the write. The Python
   client does not cache.

9. **The manifest names the study, not the trait.** Trait, sample size and ancestry are not
   fetched. Read them from the `-meta.yaml` beside the harmonised file or from the REST record
   ([references/study_metadata.md](references/study_metadata.md)); the REST field is
   `fullPvalueSet`, and it is the depositor's declaration, not proof that a file with an index
   is on the FTP.

## Operational notes

- **Coordinates are GRCh38, 1-based inclusive on input**; the tabix query is issued as
  `start_bp - 1` to `end_bp`, which is the 0-based half-open form of the same window.
- **Requests per fetch, besides the tabix reads themselves.** One HEAD (plus a directory
  listing when the simple name is absent) and a 32 KB `Range` read of the file's head to recover
  the column names, because GWAS-SSF files do not `#`-prefix their header line and tabix
  reports no header for them.
- **Timing, one measurement each, 2026-09-13, on a laptop.** Demo fetch 4.1 s; the same
  window served from the cache 0.0 s.
- **Attribution.** Cite the GWAS Catalog (Sollis 2023) and the depositing study's publication;
  the PubMed id is in the REST record and, for prefixed file names, in the file name itself.

## Safety

**Not for clinical decisions.** The output is research-grade summary statistics from a public
database; it is not a basis for diagnosis, treatment or risk communication to a person.

**Winner's curse.** Effect sizes of variants discovered in the same study are inflated; an MR
estimate built on them needs independent instruments or a correction.

**Ancestry.** Effects from a single-ancestry cohort do not transfer to other ancestries without
validation; the study's ancestry is in its metadata, not in this skill's output.

## Agent boundary

The skill returns the harmonised summary statistics of every variant in one window of one
study. The agent should:

- **Feed the rows to colocalisation, fine-mapping or MR tooling**, and not infer causation from
  a p-value.
- **Not select variants by p-value alone**; the window is returned whole so that downstream
  methods can use it whole.
- **Harmonise alleles before comparing effects across sources** (Gotcha 6).
- **Name the study in full in the reply**: accession, trait, sample size and ancestry from the
  study's metadata, and the `harmonised_url` from the manifest, not the bare `GCST` id.
- **Say what was not fetched**: a raised index error (Gotcha 2) means the study's data exist and
  could not be range-read, which is a different statement from "no summary statistics".

## References

Read only when the question needs them:

- [references/harmonised_files.md](references/harmonised_files.md): the FTP layout, the naming
  and column-layout conventions with the 17-file sample, the `-meta.yaml` sidecar, what the
  harmoniser does to each row, and the `hm_code` table.
- [references/effect_allele_harmonisation.md](references/effect_allele_harmonisation.md):
  why cross-source effect comparison needs an allele check, the failure modes, and what this
  skill's columns give you for it.
- [references/study_metadata.md](references/study_metadata.md): finding a GCST, the REST
  record's fields, the sidecar as a network-light alternative, and mapping a study's ancestry
  to a 1000 Genomes super-population.

## Sources & licensing

- NHGRI-EBI GWAS Catalog summary statistics, EBI FTP: available under
  [EMBL-EBI's terms of use](https://www.ebi.ac.uk/about/terms-of-use); cite the Catalog and
  the depositing study. Sollis E, Mosaku A, Abid A, et al. The NHGRI-EBI GWAS Catalog:
  knowledgebase and deposition resource. Nucleic Acids Res. 2023;51(D1):D977-D985.
  doi:10.1093/nar/gkac1010. PMID: 36350656.
- Harmonisation pipeline: EBISPOT
  [gwas-sumstats-harmoniser](https://github.com/EBISPOT/gwas-sumstats-harmoniser) and its
  [documentation](https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/). File
  format: Hayhurst J, Buniello A, Harris L, et al. A community driven GWAS summary statistics
  standard. bioRxiv. 2022. doi:10.1101/2022.07.15.500230 (preprint).
- Demo locus: Musunuru K, Strong A, Frank-Kamenetsky M, et al. From noncoding variant to
  phenotype via SORT1 at the 1p13 cholesterol locus. Nature. 2010;466(7307):714-719.
  doi:10.1038/nature09266. PMID: 20686566.
- pysam (MIT) and htslib (MIT): the tabix reader.
- Skill code: **MIT**.
