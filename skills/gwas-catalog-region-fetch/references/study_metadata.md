# Study metadata: finding a GCST and describing it

The skill takes a `GCST` accession and returns rows; it fetches no trait, sample-size or
ancestry information. This page is how to get the accession in, and the description out.

## Finding a GCST

1. **GWAS Catalog search.** [www.ebi.ac.uk/gwas](https://www.ebi.ac.uk/gwas/) searches by
   trait, gene, variant or author; each study page shows its accession and whether full summary
   statistics are available.
2. **REST API by trait.** `https://www.ebi.ac.uk/gwas/rest/api/efoTraits/<EFO id>/studies`
   lists the studies annotated to one EFO term. The Catalog's API documentation is at
   [www.ebi.ac.uk/gwas/docs/api](https://www.ebi.ac.uk/gwas/docs/api).
3. **The FTP listing.** `harmonised_list.txt` at the root of
   `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/` names every harmonised file
   (156,116 on 2026-09-13). It is the one list that says which studies this skill can reach at
   all; the REST flag below does not.

## The REST study record

`https://www.ebi.ac.uk/gwas/rest/api/studies/<GCST>` returns JSON. Top-level keys observed for
GCST90269602 on 2026-09-13: `accessionId`, `ancestries`, `cohort`, `diseaseTrait`,
`fullPvalueSet`, `genotypingTechnologies`, `gxe`, `gxg`, `imputed`, `initialSampleSize`,
`platforms`, `pooled`, `publicationInfo`, `qualifier`, `replicationSampleSize`, `snpCount`,
`studyDesignComment`, `userRequested`.

| Field | Value for GCST90269602 | Use |
|---|---|---|
| `diseaseTrait.trait` | `Cholesterol in medium VLDL (UKB data field 23505)` | Trait label |
| `initialSampleSize` | `88,329 European ancestry individuals` | Free text; parse with care |
| `ancestries[].type` / `numberOfIndividuals` / `ancestralGroups[].ancestralGroup` | `initial` / `88329` / `European` | Structured N and ancestry |
| `publicationInfo.pubmedId` | `36764567` | Citation of the depositing study |
| `fullPvalueSet` | `true` | The depositor declared full summary statistics |

`fullPvalueSet` is a declaration. Whether a harmonised file with a tabix index is on the FTP is
answered by the directory listing (`harmonised_files.md`), and the two disagree: in the 17-file
sample read on 2026-09-13, 4 harmonised files had no index. There is no `hasSummaryStats` key in
this record.

## The sidecar instead of the API

Beside every harmonised file sits `<name>.h.tsv.gz-meta.yaml` (3 KB for GCST90019016, one GET,
no API). It carries `trait_description`, `ontology_mapping` (EFO ids, sometimes empty),
`samples[].sample_ancestry_category` and `sample_size`, `sex`, `genome_assembly`,
`coordinate_system`, `is_harmonised`, `is_sorted` and `data_file_md5sum`. For the three bundled
demo studies it gave, on 2026-09-13: GCST90269602 cholesterol in medium VLDL, N = 88,329,
European; GCST90691573 C-reactive protein levels, N = 400,094, European; GCST90691576 glycated
haemoglobin (HbA1c) levels, N = 400,825, European. `is_sorted: false` on GCST90019016 is
consistent with its missing index: tabix requires a sorted file.

## What to say in a reply

Expand the study every time a number from it is quoted: accession, trait, N (cases and
controls for a binary trait), ancestry, and the `harmonised_url` from the manifest. A bare
`GCST90269602` tells the reader nothing; `GCST90269602, cholesterol in medium VLDL, N = 88,329,
European ancestry (UK Biobank)` does.

## Mapping ancestry to an LD reference

When the rows feed a regional plot or an instrument-pruning step, the LD panel must match the
study's ancestry. A companion skill, `ld-1000g-region-compute`, computes r² from 1000 Genomes
Phase 3 by super-population. The mapping from the Catalog's `ancestralGroup` labels:

| GWAS Catalog ancestry | 1000G super-population | Note |
|---|---|---|
| European | EUR | Finnish-only cohorts are an isolate; say so |
| African | AFR | AFR includes admixed samples (ASW, ACB); African American cohorts are admixed |
| East Asian | EAS | |
| South Asian | SAS | |
| Hispanic or Latin American | AMR | 1000G AMR is itself admixed; r² is approximate |
| Several groups (a multi-ancestry meta-analysis) | one run per group | do not pick one and present it as the study's LD |

The `ancestries[]` block can list several `initial` entries; use all of them, not the first.
