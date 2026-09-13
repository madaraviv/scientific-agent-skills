# The 1000 Genomes Phase 3 GRCh38 panel as this skill uses it

## What is fetched

Two files from the EBI 1000 Genomes FTP, both read by `scripts/ondemand_client.py`:

| File | URL | Role |
|---|---|---|
| Per-chromosome phased VCF, GRCh38, biallelic SNVs and indels (NYGC re-called, release 2019-03-12) | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr{N}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz` | Tabix-indexed; only the requested window is read (pysam `VariantFile.fetch` over HTTPS) |
| Sample panel | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel` | Maps each sample to a population and super-population; drives `--keep` |

Contig names are bare (`1`, not `chr1`); the client strips a `chr` prefix from the input.
Positions are GRCh38, 1-based. There is no liftover: a GRCh37 coordinate is silently absent.

The manifest names the panel as `panel_id: 1000g_phase3_v5b_grch38_basic` and
`panel_version: 5b_remote_2019_03_12`. These are labels the skill assigns to the URLs above,
not strings published by the consortium; quote them together with the release date.

## Super-populations

Counted from the sample panel file fetched on 2026-09-13 (2,504 samples in total):

| Code | Label | n | Populations |
|---|---|---|---|
| EUR | European | 503 | CEU, TSI, FIN, GBR, IBS |
| AFR | African | 661 | YRI, LWK, GWD, MSL, ESN, ASW, ACB |
| AMR | Admixed American | 347 | MXL, PUR, CLM, PEL |
| EAS | East Asian | 504 | CHB, JPT, CHS, CDX, KHV |
| SAS | South Asian | 489 | GIH, PJL, BEB, STU, ITU |

`super_pop` is the only ancestry control the skill has; it selects samples, nothing finer. A
single population (say FIN alone) needs a different keep file than this skill writes.

## Panel file versus region VCF

The two files come from different releases and do not list the same samples. Measured on the
cached demo window, 2026-09-13: the region VCF's `#CHROM` header carries 2,548 sample columns
(the 2019-03-12 release; plink's log reads `2548 people ... loaded from .fam`), while the
2013 panel file above lists 2,504. The keep file is written from the panel file, so 45 VCF
samples appear in no keep file and never enter a computation, and one panel sample, NA18498
(YRI, AFR), is absent from the VCF. Effective n per super-population is therefore the table
above except AFR, which computes on 660 (plink `--keep` for EUR reports 503 people
remaining, so EUR is unaffected).

## Shape of a fetched region

The bundled demo window (chromosome 1, 108,774,968 to 109,774,968, 1 Mb around the SORT1
lead rs12740374), fetched 2026-09-13: 26,405 variant rows, 1,779 of them indels, none
multi-allelic (the release is split to biallelic rows), every ID column `.`, 3.15 MB gzipped
with 2,548 sample columns. The ID column matters: plink's `--set-missing-var-ids` only names
variants whose ID is missing, and in this release that is all of them.

## Variants the panel cannot answer for

- Variants absent from 1000G Phase 3: rare alleles seen only in larger sequencing cohorts,
  array-specific probes, many structural variants. The panel returns nothing for them and the
  client does not note the gap (see `references/r2_computation.md`).
- Variants monomorphic in the chosen super-population: no r² is defined and plink emits no
  row.
- Rare variants that are present: r² is computed from very few carriers and is unstable.

## Licence and attribution

The 1000 Genomes data are open access; the IGSR data-use statement asks for attribution.
Cite both the Phase 3 paper and the IGSR resource paper whenever an r² from this skill is
reported:

- 1000 Genomes Project Consortium, Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, et
  al. A global reference for human genetic variation. Nature. 2015;526(7571):68-74.
  doi:10.1038/nature15393. PMID: 26432245.
- Clarke L, Fairley S, Zheng-Bradley X, Streeter I, Perry E, Lowy E, et al. The international
  Genome sample resource (IGSR): a worldwide collection of genome variation incorporating the
  1000 Genomes Project data. Nucleic Acids Res. 2017;45(D1):D854-D859. doi:10.1093/nar/gkw829.
  PMID: 27638885.

A methods sentence that covers both the data and the tool:

```
LD r² was computed with plink 1.9 (Chang 2015) on the 1000 Genomes Phase 3 GRCh38
phased release of 2019-03-12 (Auton 2015; Clarke 2017), restricted to the <CODE>
super-population (n = <N>).
```
