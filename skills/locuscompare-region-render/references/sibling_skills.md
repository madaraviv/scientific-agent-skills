# Sibling skills, stand-in clients and the manifest block

## How the siblings are found

`scripts/locuscompare_region_render.py` computes `skills/` as the grandparent of its own
directory and, for each sibling it needs, inserts `skills/<sibling>/scripts/` at the front of
`sys.path` when that directory exists. It then checks that each required module can be found
and raises

```
ImportError: locuscompare-region-render needs the sibling skill 'eqtl-catalogue-region-fetch'
(Python module 'eqtl_catalogue_region_fetch'), expected at <skills>/eqtl-catalogue-region-fetch/scripts.
Install that skill beside this one under <skills>, or add its scripts/ directory to sys.path or PYTHONPATH.
```

for the first one that cannot. `ukb-ppp-region-fetch` is optional: when it is absent
`UKB_PPP_AVAILABLE` is `False`, `UKBPPPClient` and `UKBPPPRegionResult` are `None`, and a
pQTL exposure raises `Tier2NotAvailable` (API) or is refused with exit code 2 naming the skill
(CLI). The test suite honours the same rule and additionally reads
`LOCUSCOMPARE_SIBLING_SKILLS_ROOT`, a directory with the same `<sibling>/scripts/` layout, for
running against siblings that are not installed in the checkout.

## What the composer calls

| Sibling | Call | Returns | Failure mapped to `Tier2NotAvailable` |
|---|---|---|---|
| `eqtl-catalogue-region-fetch` | `EQTLCatalogueClient().fetch_region(dataset_id=, chromosome=, start_bp=, end_bp=, gene_id=<molecular_trait_id>)` | `RegionResult` with `.variants` (`variant_id`, `chromosome`, `position`, `ref`, `alt`, `beta`, `se`, `p_value`, `maf`, `effect_allele_frequency`), `.release` (`dataset_release`, `study_label`, `tissue_label`, `condition_label`, `sample_group`, `quant_method`) and `.notes` | `EQTLCatalogueAPIError`, or zero variants |
| `gwas-catalog-region-fetch` | `GWASCatalogClient().fetch_region(accession=, chromosome=, start_bp=, end_bp=)` | `RegionResult` with the same variant fields plus `odds_ratio`, `.release.fetched_at_utc`, `.notes` | `GWASCatalogFetchError`, or zero variants |
| `ld-1000g-region-compute` | `OnDemand1000GLDClient(super_pop=).r2_with_lead(lead=, partners=[every exposure variant id except the lead], chromosome=, window_bp=)` | `OnDemandLDResult` with `.pairs` (`partner_variant_id`, `r2`), `.plink_version`, `.panel_id`, `.panel_version` | not raised: `LDComputeError` becomes grey points plus a caveat; `ld_client=None` does the same |
| `ukb-ppp-region-fetch` (optional) | `UKBPPPClient().fetch_region(protein_label=, ancestry=, chromosome=, start_bp=, end_bp=)` | `RegionResult` with the same variant fields and a release carrying `release_label`, `synapse_id`, `protein_label`, `protein_hgnc`, `olink_panel`, `olink_reagent_id`, `ancestry`, `ancestry_label`, `n_samples` | `UKBPPPAccessError`, or zero variants |

The gene track is not a sibling: `scripts/_fetchers/gencode_ondemand.py` calls Ensembl REST
`overlap/region` (GRCh38) and caches the JSON per region; a `requests` error or a bad payload
becomes a figure without a gene track plus the caveat `gene track unavailable (Ensembl REST
error)`.

The composer joins the exposure and outcome variants with
`harmonise_regions_for_locuscompare` (same alleles: keep; swapped: negate the outcome beta;
A/T or G/C: flag palindromic; anything else: drop), then builds the render input and the
manifest block. Order of operations: exposure fetch, outcome fetch, LD, join, gene track,
render, manifest.

## Stand-in clients for pre-fetched inputs

`scripts/_prefetched.py` gives the CLI's `sumstats_path`, `ld.source: synthetic` and
`gene_track.source: synthetic` paths objects with the same method signatures, built from TSVs
and returning the siblings' own result dataclasses (which is why the required siblings must be
installed even for an offline run):

- `PrefetchedEQTLClient(variants, dataset_id, study_label, tissue_label, quant_method)`: the
  release is tagged `prefetched`; `quant_method` defaults to `ge`, so no credible-set caveat
  fires unless you set it.
- `PrefetchedGWASClient(variants, accession)`; `gwas_variants_from_eqtl` reshapes loaded rows
  into the outcome dataclass.
- `PrefetchedLDClient(r2_by_partner)`: reports `panel_id="synthetic"`,
  `plink_version="prefetched"`, and only the partners present in the matrix.
- `load_sumstats_tsv`, `load_synthetic_ld`, `load_synthetic_gene_track`: the three loaders,
  each raising `PrefetchedSchemaError` with file and line on a malformed input.

Any object with the same `fetch_region` / `r2_with_lead` signature works as a client; the test
suite replays recorded JSON cassettes this way.

## Entering from an Open Targets colocalisation row

`render_tier2_for_lead(lead_variant_id=, chromosome=, lead_position_bp=, window_bp=,
study_mapping=StudyIdMapping(...), eqtl_client=, gwas_client=, ld_client=, out_path=,
ot_release=, super_pop=, intersected_pip_product=, extra_caveats=, ukb_ppp_client=)` builds a
`LocusCompareSpec` from a row and calls the same pipeline. `StudyIdMapping` carries
`ot_left_study_id` (the QTL side; a `UKB_PPP_<ancestry>_<protein>` id dispatches to the pQTL
fetcher, anything else to the eQTL Catalogue), `ot_right_study_id`, `gwas_catalog_accession`,
and either `eqtl_catalogue_dataset_id` or `ukb_ppp_protein_label` + `ukb_ppp_ancestry`, plus
`outcome_trait_label`, `exposure_gene_symbol` and `notes`. The Ensembl gene id is taken from the
trailing `_ensg<digits>` of the QTL study id. `load_study_id_mappings(yaml_path)` reads a
`mappings:` list of such rows keyed by (left, right) study id. `intersected_pip_product`, when
given, is printed in the window caption as `PIP_L x PIP_R = <value>`.

## `render_block` in `manifest.yaml`

| Key | Value |
|---|---|
| `ot_release` | `provenance.ot_release` from the config (`release_tag` on the spec), else empty |
| `exposure_source` | `eqtl_catalogue` (also for `sumstats_path` inputs) or `ukb_ppp` |
| `exposure_source_release` | the eQTL Catalogue dataset release (`v7+` when the fetcher reports none), `prefetched` for a TSV, or the UKB-PPP release label |
| `exposure_study_id` | dataset id, the config's `exposure.study_id`, or the UKB-PPP Synapse id |
| `exposure_protein_label`, `exposure_ancestry`, `exposure_ancestry_label` | pQTL renders only; empty strings otherwise |
| `exposure_harmonisation_version`, `outcome_harmonisation_version` | always empty; neither source reports one |
| `outcome_source` | always `gwas_catalog_harmonised` (also for `sumstats_path` inputs) |
| `outcome_source_release` | the date part of the outcome fetch timestamp |
| `outcome_study_id` | the accession, or the config's `outcome.study_id` |
| `ld_panel`, `ld_panel_super_pop`, `ld_panel_version`, `plink_version` | from the LD result; `none` / empty when LD was unavailable; `synthetic` / `prefetched` for a matrix file. `ld_panel_super_pop` is the requested super-population even when no LD ran |
| `window_bp`, `lead_variant_id` | as configured |
| `n_pairs` | joined pairs after allele reconciliation, palindromic pairs included |
| `n_palindromic_excluded` | joined A/T or G/C pairs |
| `scatter_downsampled`, `scatter_downsample_target` | `n_pairs > 5000`, and 5000 |
| `ancestry_caveats` | the curated caveat list: LD fallback, gene-track fallback, the credible-set-filter note, and the config's `caveats:` |
| `data_source_warnings` | each fetcher's `notes`, prefixed `eqtl_catalogue:`, `gwas_catalog:` or `ukb_ppp:`; for TSV inputs this holds `sumstats loaded from pre-fetched TSV (no live fetch)` |
| `plot_artifact` | the PNG basename |
| `fetched_at` | render time, UTC |
