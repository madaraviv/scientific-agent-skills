# Choosing `super_pop` for a study

LD is ancestry-specific: the same pair of variants can be tightly linked in one population
and independent in another. The r² this skill returns describes the 1000 Genomes samples of
the chosen super-population, and it is only informative about a study to the extent that the
study's participants share that ancestry. Pick `super_pop` from the study's ancestry, never
from convenience.

## Procedure

1. Read the study's ancestry from its own metadata: the eQTL Catalogue exposes a population
   field per dataset; the GWAS Catalog lists ancestral groups per study with participant
   counts; a consortium GWAS states it in the paper.
2. Map it to one of the five 1000G super-populations with the table below.
3. Pass that code as `super_pop`, and carry the caveat column into whatever you render.

When a study mixes ancestries and reports the weights, call the skill once per
super-population that carries a substantial share and report each panel separately with its
weight; a fixed cut such as 25 percent is a workable rule, provided the reply says which
groups were left out and why. One panel cannot stand in for a mixture.

## Mapping table

| Ancestry as stated by the source | `super_pop` | Caveat to carry |
|---|---|---|
| European, EUR, British, German, Scandinavian, Italian, Iberian | EUR | none |
| Finnish (FinnGen and other Finnish cohorts) | EUR | Finnish-specific LD; see below |
| African (continental), Yoruba, Luhya, Gambian, Esan, Mende | AFR | none |
| African American, African Caribbean | AFR | the AFR panel includes ASW and ACB, both admixed; the study's admixture may differ from the panel's |
| East Asian, Han Chinese, Japanese, Korean, Vietnamese | EAS | none when the population is named; flag a bare "East Asian" |
| South Asian, Indian, Pakistani, Bangladeshi, Sri Lankan | SAS | none |
| Hispanic, Latino, Mexican, Puerto Rican, Colombian, Peruvian, AMR | AMR | the AMR samples are admixed with varying European, Native American and African ancestry; r² reflects the panel's mixture, not the study's |
| Multi-ancestry with weights | one call per substantial group | label each panel with its weight |
| Unspecified, mixed without weights, blank | EUR as an explicit fallback | say so: "study ancestry not stated; 1000G EUR used; LD is approximate" |

## Finnish cohorts

The 1000G EUR panel includes 99 FIN samples among its 503, so Finnish LD is represented but
diluted. Finland's population history (a founder bottleneck followed by expansion) gives it a
rare-variant spectrum distinct from the rest of Europe (Locke 2019), and that is where the
EUR proxy is weakest; for common lead variants use the EUR panel and state the caveat. For rare-variant fine-mapping in a Finnish cohort use a Finnish reference,
which this skill does not provide.

## Admixed cohorts

1000G's AMR super-population, and the ASW and ACB populations inside AFR, are themselves
admixed, so the panel does carry some admixed LD. It carries the panel's admixture, not the
study's: at loci where local ancestry varies between individuals (the HLA region, for
example) the study cohort's LD can differ from the panel's. Say this in the reply when
`super_pop` is AMR or when the study is African American; for a cohort with a specific,
known admixture that matches none of the five groups, say that no 1000G super-population is a
match rather than silently picking one.

## When 1000G is the wrong reference

- Populations not in 1000G (Native American, Pacific Islander, Roma, Indigenous Australian,
  among others): state that no panel is available rather than substituting one.
- Rare-variant work: with 347 to 661 samples per super-population, r² at MAF below 0.01
  rests on a handful of carriers. Use a larger sequenced reference (TOPMed, gnomAD), which is
  outside this skill.
- Instrument selection for Mendelian randomisation where the GWAS's own LD reference is
  available: prefer the in-sample reference; use the 1000G proxy for visualisation.

## Citations

- 1000 Genomes Project Consortium, Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, et
  al. A global reference for human genetic variation. Nature. 2015;526(7571):68-74.
  doi:10.1038/nature15393. PMID: 26432245.
- Locke AE, Steinberg KM, Chiang CWK, Service SK, Havulinna AS, Stell L, et al. Exome
  sequencing of Finnish isolates enhances rare-variant association power. Nature.
  2019;572(7769):323-328. doi:10.1038/s41586-019-1457-z. PMID: 31367044.
