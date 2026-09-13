---
name: locuscompare-region-render
description: Render the four-panel regional LocusCompare diagnostic (Liu 2019 convention) for one lead variant. Stacks a GWAS Manhattan track over a QTL Manhattan track and a gene track, then draws a cross-trait -log10(p) scatter and an effect-size scatter, every point coloured by LD r² to the lead from 1000 Genomes. Joins an exposure slice (eQTL Catalogue fetch or a harmonised TSV) to an outcome slice (GWAS Catalog fetch or a harmonised TSV) on variant id, flips swapped alleles, flags palindromic SNPs, and writes a PNG plus a provenance manifest and a short report. Use when an agent needs a visual, auditable check that two association signals at a locus share one causal variant, as a companion to a colocalisation posterior it already has.
license: MIT
compatibility: Python 3.10+ with matplotlib, numpy, PyYAML and requests. Needs the sibling skills eqtl-catalogue-region-fetch, gwas-catalog-region-fetch and ld-1000g-region-compute installed beside it under skills/ (ukb-ppp-region-fetch is optional). Network access to the EBI FTPs and Ensembl REST for live fetches; none for pre-fetched TSV inputs. plink 1.9 on PATH for LD colouring, otherwise points render grey.
metadata:
  version: "1.0"
  skill-author: Aviv Madar
---

# LocusCompare Region Render

## Overview

A colocalisation posterior says how likely two association signals share one causal variant;
it does not show you the locus. This skill draws it: the four-panel regional LocusCompare figure
of Liu 2019 for one lead variant, an exposure study and an outcome study. From top to bottom the
figure stacks the outcome (GWAS) Manhattan track, the exposure (QTL) Manhattan track and a gene
track on one shared genomic axis, then puts two scatters side by side underneath: the
LocusCompare panel (`-log10 p` exposure against `-log10 p` outcome, with the identity line) and
the effect-size panel (`β` exposure against `β` outcome, with the slope through the lead). Every
variant is coloured by its r² to the lead in the LocusZoom five-bin palette; the lead is a
purple diamond. A six-line caption under the figure records the LD panel, the window, both
sources and the caveats, so the PNG carries its own provenance.

The skill is a composer. It calls three sibling skills for the data (the eQTL Catalogue and
GWAS Catalog region fetchers and the 1000 Genomes LD client), harmonises the two slices on
`variant_id`, renders, and writes `manifest.yaml` and `report.md` next to the PNG. Each side can
instead be a harmonised TSV on disk, in which case nothing is fetched and a run takes well under
a second (the bundled synthetic demo rendered in 0.51 s wall clock on one laptop run, 200
variants). It does not compute a colocalisation posterior, fine-map, or estimate a causal effect;
those numbers come from upstream tools and this figure is read next to them.

## Trigger

**Fire when** the user (or an upstream agent step) wants:

- A regional four-panel view of a colocalisation candidate, given a lead variant, an eQTL
  Catalogue dataset (or a harmonised exposure TSV) and a GWAS Catalog accession (or a harmonised
  outcome TSV).
- A visual check on an Open Targets colocalisation row (a `left_studyId` / `right_studyId` pair
  with its lead) before quoting its H4 in a write-up.
- The picture to go with an eQTL or sQTL, single-cell eQTL, or (with the optional pQTL sibling)
  plasma-protein QTL joined to a GWAS at one locus.
- A LocusZoom-style regional plot of two traits with LD colouring and a gene track, from
  summary statistics already on disk.

Trigger phrases include "regional locuscompare", "locuscompare plot", "regional coloc plot",
"GWAS eQTL overlay", "colocalization plot for a lead variant", "do these two signals share a
causal variant".

**Do NOT fire when** the user wants:

- **A colocalisation posterior** (H3 / H4). Take it from Open Targets, or run coloc, coloc-SuSiE
  or another colocalisation tool; this skill only draws the data and can display the number
  you pass it in a caption.
- **A causal-effect estimate from the two slices.** The effect-size panel draws the ratio of
  the lead's two betas as a slope for orientation only; no standard error, instrument check or
  sensitivity analysis is computed here.
- **A single-trait regional plot.** The figure is defined by the join of two slices; with one
  slice there is no LocusCompare panel.
- **Fine-mapping or credible sets.** Nothing here computes posterior inclusion probabilities.
- **A lead variant.** The lead is an input. Without one (from a colocalisation row, a fine-mapping
  result or a lookup), there is nothing to centre the window on.
- **A trans-QTL.** The exposure window is centred on the lead; a trans signal sits outside the
  cis window and the exposure track would be empty.

## Quick start

Install the three required sibling skills beside this one (`skills/eqtl-catalogue-region-fetch`,
`skills/gwas-catalog-region-fetch`, `skills/ld-1000g-region-compute`); the composer resolves
each at `skills/<sibling>/scripts/` and raises an `ImportError` naming the missing skill
otherwise. Then:

```bash
# Offline demo: pre-fetched TSVs, a synthetic LD matrix, a synthetic gene track. No network.
python scripts/cli.py --demo 01_synthetic_demo --output ./out/synthetic

# Live demo: SORT1 eQTL (GTEx minor salivary gland) x cholesterol in medium VLDL.
# Fetches from the eQTL Catalogue FTP, the GWAS Catalog FTP, the 1000 Genomes FTP (plink 1.9
# needed for the LD step) and Ensembl REST (gene track). Bare --demo picks this one.
python scripts/cli.py --demo --output ./out/sort1_vldl

# Your own config (JSON or YAML)
python scripts/cli.py --input config.yaml --output ./out/my_locus

python scripts/cli.py --list-demos
```

Seven demos ship under `scripts/examples/` (details in `scripts/examples/README.md`):

| Demo | Exposure | Network | Needs |
|---|---|---|---|
| `01_synthetic_demo` | synthetic TSVs | no | nothing beyond the required siblings |
| `02_eqtl_catalogue_x_gwas_catalog` | SORT1 gene expression, GTEx minor salivary gland (QTD000276) | yes | plink 1.9 for LD colours |
| `03_open_targets_followup` | LDLR transcript expression, GTEx skin (QTD000318), from an Open Targets row | yes | same |
| `04_gwas_lookup_followup` | SORT1 gene expression (QTD000276), entered from an rsID lookup | yes | same |
| `05_sqtl_sort1_liver_txrev` | SORT1 transcript usage, GTEx liver (QTD000269) | yes | same |
| `06_sceqtl_sort1_onek1k_cd14_mono` | SORT1 gene expression, OneK1K CD14+ monocytes (QTD000609) | yes | same |
| `07_pqtl_sort1_ukbppp_eur` | SORT1 plasma protein, UKB-PPP European | yes | the optional `ukb-ppp-region-fetch` skill, which is not part of this collection; the run is refused by name without it |

Every live demo shares the outcome GCST90269602 (cholesterol in medium VLDL) except 03.

## Config schema

Two blocks are required, `lead` and both of `exposure` / `outcome`. Each side is supplied
either by a `fetch:` block (live tabix fetch through a sibling skill) or by `sumstats_path:`
(a harmonised TSV in the schema of [references/input_schema.md](references/input_schema.md));
the two sides may mix.

```yaml
lead:
  variant_id: "1_109274968_G_T"   # chr_pos_ref_alt, GRCh38, no "chr" prefix; the join key
  rs_id: "rs12740374"             # optional; manifest and report only, never used to join
  chromosome: "1"
  position_bp: 109274968
  window_bp: 1000000              # full width; the window is +/- window_bp/2 around the lead
exposure:
  trait_label: "SORT1 expression, minor salivary gland"
  gene_symbol: "SORT1"            # optional; this gene is drawn bold and red on the gene track
  fetch:
    source: eqtl_catalogue        # or ukb_ppp (needs the optional sibling)
    dataset_id: QTD000276
    molecular_trait_id: ENSG00000134243   # passed to the fetcher as the gene filter
outcome:
  trait_label: "cholesterol in medium VLDL"
  fetch:
    source: gwas_catalog
    accession: GCST90269602
ld:                               # optional
  source: 1000g_phase3_grch38     # default; or synthetic (+ ld_matrix_path)
  super_pop: EUR                  # EUR / AFR / AMR / EAS / SAS
  plink_bin: plink                # optional; PLINK_BIN env var is read by the LD sibling too
gene_track:                       # optional
  source: gencode_v39             # default: Ensembl REST on demand; or synthetic (+ genes_path)
  gtf_path: null                  # accepted, not read: the track is REST or synthetic (see below)
  biotypes: [protein_coding]      # default
provenance:                       # optional; folded into the caption and the manifest
  ot_release: "26.03"
  gwas_lookup_run_dir: "runs/gwas_lookup/"
caveats:                          # optional; printed in the caption and stored in the manifest
  - "FinnGen cohort; 1000G EUR used as the LD reference"
```

| Key | Type | Notes |
|---|---|---|
| `lead.variant_id` | string | Required. Must appear in both slices to be drawn as the lead; the join does not depend on it. |
| `lead.chromosome`, `lead.position_bp`, `lead.window_bp` | string, int, int | Required. Fetch window is `position_bp +/- window_bp // 2`. |
| `exposure.fetch.source` | `eqtl_catalogue` or `ukb_ppp` | `eqtl_catalogue` needs `dataset_id`; `ukb_ppp` needs `protein_label` (and `ancestry`, default EUR) and the optional sibling. |
| `exposure.sumstats_path` / `outcome.sumstats_path` | path | Relative paths resolve against the config file's directory. `exposure.study_id` / `outcome.study_id` label the manifest (default `prefetched`). |
| `outcome.fetch.source` | `gwas_catalog` only | With `accession` (`GCST...`). |
| `ld.source: synthetic` | | `ld_matrix_path` is a two-column TSV (`partner_variant_id`, `r2`). |
| `gene_track.source: synthetic` | | `genes_path` is a TSV (`gene_symbol`, `start`, `end`, `strand`, optional `biotype`). |

Keys the CLI accepts silently and does not act on: `schema_version`, `harmonisation`
(palindromic exclusion is always on; see Gotchas), `ld.region_vcf_path`, the chain demo's
`fine_mapping_outputs`, and `gene_track.gtf_path` / the `GENCODE_GTF` environment variable (a
local GTF is never parsed; setting it only makes the CLI skip its own Ensembl call, after which
the composer performs the same call). They are in the bundled configs for documentation only.

## Output files

```
<output>/
├── <lead_variant_id>_full_locuscompare.png   # the figure, 160 dpi, 13 in wide
├── manifest.yaml                             # provenance, counts, render_block
└── report.md                                 # lead, both labels, n_pairs, n_palindromic_excluded, notes
```

`manifest.yaml` top level: `lead_variant_id`, `lead_rs_id`, `n_pairs` (variants joined across
the two slices after allele reconciliation), `n_palindromic_excluded` (joined A/T or G/C SNPs,
kept in the Manhattan tracks and out of the two scatters), `plot_path`, `notes` (LD and
gene-track fallbacks), and `render_block`. The `render_block` keys are `ot_release`,
`exposure_source`, `exposure_source_release`, `exposure_study_id`, `exposure_protein_label`,
`exposure_ancestry`, `exposure_ancestry_label`, `exposure_harmonisation_version`,
`outcome_source`, `outcome_source_release`, `outcome_study_id`,
`outcome_harmonisation_version`, `ld_panel`, `ld_panel_super_pop`, `ld_panel_version`,
`plink_version`, `window_bp`, `lead_variant_id`, `n_pairs`, `n_palindromic_excluded`,
`scatter_downsampled`, `scatter_downsample_target`, `ancestry_caveats` (the curated caveat list,
including the credible-set-filter note for non-gene-expression datasets and the grey-points
note when LD is unavailable), `data_source_warnings` (each fetcher's own notes, prefixed
`eqtl_catalogue:` / `gwas_catalog:` / `ukb_ppp:`), `plot_artifact` and `fetched_at`.
Field-by-field detail is in [references/sibling_skills.md](references/sibling_skills.md).

## Python API

`scripts/locuscompare_region_render.py` exposes the composer directly:

```python
# from the skill's scripts/ directory, or with it on sys.path
from locuscompare_region_render import LocusCompareSpec, Tier2NotAvailable, render_locuscompare_for_lead
from eqtl_catalogue_region_fetch import EQTLCatalogueClient
from gwas_catalog_region_fetch import GWASCatalogClient
from ondemand_client import OnDemand1000GLDClient   # from ld-1000g-region-compute

spec = LocusCompareSpec(
    lead_variant_id="1_109274968_G_T", chromosome="1", lead_position_bp=109_274_968,
    window_bp=1_000_000, eqtl_dataset_id="QTD000276", molecular_trait_id="ENSG00000134243",
    gwas_accession="GCST90269602", exposure_gene_symbol="SORT1",
    outcome_trait_label="cholesterol in medium VLDL",
)
result = render_locuscompare_for_lead(
    spec, eqtl_client=EQTLCatalogueClient(), gwas_client=GWASCatalogClient(),
    ld_client=OnDemand1000GLDClient(super_pop="EUR"),   # or None: grey points plus a caveat
    out_path="out/sort1.png",
)
result.n_pairs, result.n_palindromic_excluded, result.manifest_block, result.notes
```

`Tier2NotAvailable` is raised, with the reason in the message, when a fetcher cannot resolve
the study, a fetch returns zero variants, the two slices share no variant, or a pQTL exposure is
requested without a pQTL client. Catch it and fall back to a credible-set-only display rather
than retrying. `render_tier2_for_lead` is the same pipeline keyed by an Open Targets row through
a `StudyIdMapping` (`ot_left_study_id`, `ot_right_study_id`, `gwas_catalog_accession`,
`eqtl_catalogue_dataset_id` or the UKB-PPP protein and ancestry); it dispatches on the study id
prefix (`UKB_PPP_*` to the pQTL fetcher, anything else to the eQTL Catalogue) and
`load_study_id_mappings(yaml_path)` reads a table of such rows.

## Sibling skills

| Skill | Module | Used for | Required |
|---|---|---|---|
| `eqtl-catalogue-region-fetch` | `eqtl_catalogue_region_fetch` | exposure slice (gene expression, exon, transcript, transcript usage, splicing, microarray, single-cell datasets) | yes |
| `gwas-catalog-region-fetch` | `gwas_catalog_region_fetch` | outcome slice (harmonised summary statistics by `GCST` accession) | yes |
| `ld-1000g-region-compute` | `ld_1000g_region_compute`, `ondemand_client` | r² of every exposure variant to the lead in one 1000 Genomes super-population, via plink 1.9 | yes (the client is optional at run time) |
| `ukb-ppp-region-fetch` | `ukb_ppp_region_fetch` | plasma pQTL exposure (`source: ukb_ppp`, `UKB_PPP_*` study ids) | no; not in this collection |

The composer inserts `skills/<sibling>/scripts/` on `sys.path` for each of these when the
directory exists. A required sibling that cannot be imported raises
`ImportError: locuscompare-region-render needs the sibling skill '<name>' ... expected at
<path>`. Without `ukb-ppp-region-fetch` the module still loads with
`UKB_PPP_AVAILABLE = False`; a pQTL exposure then raises `Tier2NotAvailable` from the API and
the CLI refuses `source: ukb_ppp` with a message naming the skill. pQTL exposures are unavailable
in this collection until that skill exists here. The stand-in clients in
`scripts/_prefetched.py` still build the siblings' result types, so `sumstats_path:` inputs
need the required siblings installed too.

## Gotchas

1. **Palindromic SNPs are always excluded from the two scatters and always kept in the
   Manhattan tracks.** An A/T or G/C SNP cannot be strand-aligned across two sources without
   allele frequencies, so its joined pair is flagged, counted in `n_palindromic_excluded`, and
   left out of the LocusCompare and effect-size panels (Hemani 2018 convention). In the bundled
   SORT1 fixture that is 333 of 2,547 joined pairs (13%); in the synthetic demo it is the lead
   itself (`1_500000_A_T` is A/T), so that figure has no lead diamond in the bottom row. There is
   no switch: the `harmonisation.exclude_palindromic` key in the example configs is not read.

2. **Allele reconciliation is exact.** Same `ref`/`alt` on both sides: kept as is. Swapped
   (`ref`/`alt` on one side is `alt`/`ref` on the other): the outcome beta is negated and
   `flip_outcome_beta` is recorded. Any other allele pair at the same `variant_id` is dropped
   from the join. Indels are compared as strings, so normalise with `bcftools norm` upstream if
   the two sources represent them differently.

3. **The lead has to be in both slices to be drawn.** There is no proxy substitution: a lead
   absent from one side simply has no diamond and no annotation, while the rest of the figure
   renders and `n_pairs` counts the variants that did join. Check the figure for the diamond
   before describing the lead's position; the colours and the identity line stay meaningful.

4. **LD colouring is a fallback, not a failure.** With `ld_client=None`, without plink 1.9, or
   when the LD computation raises, every point is grey (the `r² 0.0 - 0.2` bin), `ld_panel` is
   `none`, and a caveat is written into the caption and `ancestry_caveats`. The LD sibling only
   supports plink 1.9; the super-population is the whole ancestry statement, so a FinnGen or
   Biobank Japan outcome coloured with EUR or EAS 1000 Genomes r² is an approximation you
   should name in `caveats:`.

5. **Non-gene-expression eQTL Catalogue datasets are credible-set-filtered.** For splicing,
   exon, transcript and transcript-usage datasets the eQTL Catalogue fetcher reads the
   `.cc.tsv.gz` file (credible-set rows only), so the exposure track is sparser than a
   gene-expression run. The composer writes the caveat `sumstats are credible-set-filtered
   (eQTL Catalogue .cc.tsv.gz; quant_method=<code>)` whenever the dataset's quantification
   method is not `ge` or `microarray`; quote it with the figure.

6. **Windows are centred on the lead, not on the gene.** The exposure fetch covers
   `position_bp +/- window_bp // 2`; the eQTL Catalogue tests variants within 1 Mb of the
   gene's transcription start, so an off-centre lead leaves one side of the exposure track
   empty. The renderer shades any part of the window where a source has no tested variants and
   labels it `shaded = no source data tested`; that shading is coverage, not absence of signal.

7. **`p = 0` and `p = NA` make no point.** A p-value of zero, `NA` or a non-finite value gives
   no `-log10 p`, so the variant is skipped in every panel (it is not floored to a tiny value).
   The `sumstats_path` loader drops rows missing any of `beta`, `se`, `p` before harmonisation
   and does not range-check the values it keeps.

8. **Dense windows are downsampled in the scatters only.** Past 5,000 joined non-palindromic
   pairs the two bottom panels keep every pair in the `r² 0.8 - 1.0` bin and sample the other
   bins proportionally (seeded, so repeatable); `scatter_downsampled` records it. The Manhattan
   tracks are never downsampled.

9. **`molecular_trait_id` is passed to the exposure fetcher as a gene id.** The composer calls
   the eQTL Catalogue fetcher with `gene_id=<molecular_trait_id>`, which is the right filter for
   gene-expression datasets and for the parent gene of splicing / transcript datasets (the
   fetcher matches its `gene_id` column). Passing a transcript or intron id here returns zero
   rows and a `Tier2NotAvailable`.

## Operational notes

- **Caches.** The gene-track fetcher caches each Ensembl REST region under
  `~/.clawbio/locuscompare_cache/gencode/` (the path is inherited from the code's origin);
  set `LOCUSCOMPARE_CACHE_DIR` to move it, read per call. The eQTL Catalogue and GWAS Catalog
  fetchers and the LD client keep their own caches, documented in their skills.
- **Network.** Live runs read `ftp.ebi.ac.uk` (eQTL Catalogue and GWAS Catalog harmonised
  summary statistics), `ftp.1000genomes.ebi.ac.uk` (the region VCF and the sample panel) and
  `rest.ensembl.org` (`overlap/region`, GRCh38). A `sumstats_path` run with `ld.source:
  synthetic` and `gene_track.source: synthetic` touches nothing.
- **Gene track sources.** `gene_track.source: synthetic` reads `genes_path` (no exons, so
  genes draw as lines with strand arrows); otherwise the track comes from Ensembl REST
  `overlap/region` with exons drawn as boxes. If the REST call fails the figure renders without
  the track and the caveat `gene track unavailable (Ensembl REST error)` is recorded. Only
  `protein_coding` genes are drawn unless `gene_track.biotypes` says otherwise.
- **Timing.** The offline synthetic demo (200 variants, no fetch) took 0.51 s wall clock on
  one run on a laptop; a live run is bounded by the three fetches and the plink step, which the
  sibling skills document. Time your own first live run before promising a duration.
- **Exit codes.** `cli.py` returns 0 on success, 2 for a config error (missing block,
  unsupported source, pQTL without the sibling), 1 when the render raises
  `Tier2NotAvailable` (the reason is printed as `render failed: ...`).

## Safety

**Not for clinical decisions.** The figure is an interpretation aid built from public research
summary statistics. Do not use it for diagnosis, treatment selection or any clinical decision.

**Local-first.** Inputs stay on the user's machine; the only outbound traffic is to the public
EBI, 1000 Genomes and Ensembl endpoints named above, and only for the sides configured with
`fetch:`. Nothing is uploaded.

**Provenance travels with the figure.** The caption and `manifest.yaml` record the LD panel and
super-population, the window, both study identifiers and releases, the palindromic count, the
fetchers' own warnings and every fallback taken. Quote them; a plot without its caveats is a
different claim.

## Agent boundary

The skill returns a figure and a manifest for one (lead, exposure, outcome) tuple. The agent
should:

- **Read the figure next to a colocalisation posterior, never instead of one.** A clean
  diagonal in the LocusCompare panel is consistent with one shared causal variant; two clusters
  are consistent with distinct causal variants in LD. Neither pattern is proof, and two causal
  variants in tight LD can look diagonal. See
  [references/visual_interpretation.md](references/visual_interpretation.md).
- **Say what was drawn.** In the reply, name the lead (with rs id when known), the window,
  both study identifiers with their labels, the LD panel and super-population, and
  `n_pairs` / `n_palindromic_excluded`. Expand codes: `quant_method=txrev` is "transcript
  usage", `EUR` is "European super-population (1000 Genomes)".
- **Surface every caveat verbatim**: grey points from a missing LD panel, the
  credible-set-filter note, a missing gene track, an ancestry mismatch, and the fetchers'
  `data_source_warnings`.
- **Say when the lead is not drawn** (absent from one slice, or palindromic and therefore out
  of the scatters) rather than describing a diamond that is not there.
- **Not read a slope as an effect estimate.** The dashed line in the effect-size panel is
  `β_outcome / β_exposure` at the lead, drawn for orientation; causal inference needs its own
  method and its own assumptions.
- **Not change the super-population silently.** If the cohort's ancestry does not match
  `ld.super_pop`, say so and ask before proceeding.

## References

Read only when the question needs them:

- [references/input_schema.md](references/input_schema.md): the harmonised TSV contract for
  `sumstats_path` inputs, what the loader checks and what it does not, and the four
  harmonisation recipes (FinnGen, Pan-UKBB, UKB-PPP, GTEx v10) under `scripts/examples/recipes/`.
- [references/four_panel_layout.md](references/four_panel_layout.md): panel by panel, what is
  drawn from which data, the r² palette, the downsampling rule and the caption fields.
- [references/visual_interpretation.md](references/visual_interpretation.md): reading the
  LocusCompare panel (diagonal versus two clusters), what the pattern does not prove, and the
  numbers to cite alongside it.
- [references/sibling_skills.md](references/sibling_skills.md): the calls made into each
  sibling skill, the stand-in clients for pre-fetched inputs, the Open Targets row entry point,
  and the manifest `render_block` field by field.

## Sources & licensing

- Skill code: MIT.
- eQTL Catalogue summary statistics: CC-BY 4.0 (per-study attribution in the catalogue's study
  table). Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human
  gene expression and splicing quantitative trait loci. Nat Genet. 2021;53(9):1290-1299.
  doi:10.1038/s41588-021-00924-w. PMID: 34493866.
- GWAS Catalog harmonised summary statistics: open access. Sollis E, Mosaku A, Abid A, et al.
  The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. Nucleic Acids Res.
  2023;51(D1):D977-D985. doi:10.1093/nar/gkac1010. PMID: 36350656.
- 1000 Genomes Phase 3: open access. 1000 Genomes Project Consortium, Auton A, Brooks LD, et
  al. A global reference for human genetic variation. Nature. 2015;526(7571):68-74.
  doi:10.1038/nature15393. PMID: 26432245.
- Gene track: Ensembl REST (GRCh38) / GENCODE. Frankish A, Diekhans M, Jungreis I, et al.
  GENCODE 2021. Nucleic Acids Res. 2021;49(D1):D916-D923. doi:10.1093/nar/gkaa1087.
  PMID: 33270111.
- plink 1.9 (GPL-3, run as a subprocess by the LD sibling, never bundled): Chang CC, Chow CC,
  Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge
  of larger and richer datasets. Gigascience. 2015;4:7. doi:10.1186/s13742-015-0047-8.
  PMID: 25722852.
- Figure convention: Liu B, Gloudemans MJ, Rao AS, Ingelsson E, Montgomery SB. Abundant
  associations with gene expression complicate GWAS follow-up. Nat Genet. 2019;51(5):768-769.
  doi:10.1038/s41588-019-0404-0. PMID: 31043754. Regional-plot and r² palette convention:
  Pruim RJ, Welch RP, Sanna S, et al. LocusZoom: regional visualization of genome-wide
  association scan results. Bioinformatics. 2010;26(18):2336-2337.
  doi:10.1093/bioinformatics/btq419. PMID: 20634204.
- Palindromic-SNP exclusion: Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform
  supports systematic causal inference across the human phenome. Elife. 2018;7:e34408.
  doi:10.7554/eLife.34408. PMID: 29846171.
- The SORT1 demo locus: Musunuru K, Strong A, Frank-Kamenetsky M, et al. From noncoding
  variant to phenotype via SORT1 at the 1p13 cholesterol locus. Nature. 2010;466(7307):714-719.
  doi:10.1038/nature09266. PMID: 20686566.
