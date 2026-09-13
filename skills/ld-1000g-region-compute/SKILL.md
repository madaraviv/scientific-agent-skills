---
name: ld-1000g-region-compute
description: Compute pairwise LD r² between a lead variant and a set of partner variants in a window from the 1000 Genomes Phase 3 GRCh38 reference panel, stratified by super-population (EUR / AFR / AMR / EAS / SAS). Tabix-fetches only the region from the EBI 1000 Genomes FTP and runs plink 1.9 locally, so no multi-GB panel download. Use when an agent needs r² for LD colouring of a regional plot, LD pruning of MR instruments, or a check that two nearby GWAS hits tag one signal.
license: MIT
compatibility: Requires the plink 1.9 binary on PATH or in PLINK_BIN (GPL-3, run as a subprocess, not bundled; plink 2.x is not supported) and network access to the EBI 1000 Genomes FTP for the region VCF and the sample panel; Python packages pysam and requests (PyYAML only for YAML configs).
metadata:
  version: "1.0"
  skill-author: Aviv Madar
---

# LD 1000G Region Compute

## Overview

Pairwise r² between a lead variant and a set of partners is the input to LD colouring of a
regional association plot, to LD pruning of a candidate instrument set, and to the question
"are these two nearby hits one signal or two". This skill computes it from the [1000 Genomes
Phase 3](https://www.internationalgenome.org/) GRCh38 phased release of 2019-03-12 (Auton
2015; Clarke 2017), for one super-population at a time.

It ships one client, `OnDemand1000GLDClient` (`scripts/ondemand_client.py`): tabix-fetch the
window of the per-chromosome VCF from the EBI FTP (the 1 Mb demo window is 3.15 MB), select
the super-population's samples with the consortium's panel file, and run `plink --r2` on the
result. Nothing is downloaded beyond the window: one small fetch, then a plink run that took
0.77 s on the cached demo window (one measurement, see Operational notes). The
CLI (`scripts/ld_1000g_region_compute.py`) wraps the client and writes a pairs table, a
provenance manifest, and a short report.

The supported binary is plink 1.9 (Chang 2015), invoked as a subprocess; it is GPL-3 and is
never bundled. Install it with `brew install brewsci/bio/plink`, `apt-get install plink1.9`,
or `conda install -c bioconda plink`, or download it from
[cog-genomics.org/plink/1.9](https://www.cog-genomics.org/plink/1.9/) and point `PLINK_BIN`
at it. plink 2.x names its flags differently and is not supported.

## Trigger

**Fire when** the user (or an upstream agent step) wants:

- r² between a lead variant and a list of partner variants in one chromosomal window, in a
  named 1000G super-population.
- LD colouring input for a regional plot (LocusCompare or LocusZoom style).
- LD-pruning input for Mendelian randomisation instrument selection.
- A check that two GWAS hits at nearby positions tag the same signal (high r²) or separate
  signals (low r²).
- An ancestry-matched LD reference for colocalisation or fine-mapping inputs.

Trigger phrases include "LD around lead", "r² 1000G", "1000 Genomes LD panel",
"ancestry-stratified LD", "LocusZoom LD colouring".

**Do NOT fire when** the user wants:

- **r² across several populations at once**: call this skill once per super-population and
  report each; it does not merge them.
- **LD on UK Biobank, gnomAD, TOPMed, HRC or any other panel**: 1000G Phase 3 only.
- **A precomputed genome-wide LD matrix**: this is an on-demand window compute.
- **Haplotype blocks or D' matrices**: the client reports r² only (`dprime` is always empty).
- **Cross-population LD comparisons**: use a dedicated tool such as LDlink.
- **LD at a single population finer than a super-population** (FIN alone, say): the keep
  file this skill writes is per super-population.

## Quick start

```bash
# Bundled demo: SORT1 lead rs12740374, five partners within 5 kb, 1 Mb window, EUR
python scripts/ld_1000g_region_compute.py --demo --output ./out/sort1_demo

# Your own config
python scripts/ld_1000g_region_compute.py --input my_locus.json --output ./out/my_locus

# List the bundled configs
python scripts/ld_1000g_region_compute.py --list-demos
```

The demo, run on 2026-09-13 with plink v1.90p:

```
info: using bundled demo default.json
ld-1000g-region-compute: 4 pairs -> out/sort1_demo/ld_pairs.tsv
  panel: 1000g_phase3_v5b_grch38_basic (EUR) | plink PLINK v1.90p 64-bit (6 Sep 2023)
```

Five partners were requested and four pairs came back: three partners plus the lead's own
row. The two missing partners are spelled with alleles out of ASCII order; see Gotcha 2. That
is the demo as shipped, kept so the behaviour is visible on first run.

From Python (run from `scripts/`, or put it on `sys.path`):

```python
from ondemand_client import OnDemand1000GLDClient

client = OnDemand1000GLDClient(super_pop="EUR")          # checks that plink is reachable
result = client.r2_with_lead(
    lead="1_109274968_G_T",
    partners=["1_109272630_A_G", "1_109274570_A_G", "1_109274623_C_T"],
    chromosome="1",
    window_bp=1_000_000,
)
for pair in result.pairs:
    print(pair.partner_variant_id, pair.r2)
```

`scripts/examples/run_example.sh` does the same for `scripts/examples/input.json`, and
`scripts/examples/expected_output.md` shows the measured result.

## Config schema

JSON or YAML (YAML needs PyYAML):

```json
{
  "lead": "1_109274968_G_T",
  "partners": ["1_109272630_A_G", "1_109274570_A_G", "1_109274623_C_T"],
  "chromosome": "1",
  "window_bp": 1000000,
  "super_pop": "EUR"
}
```

| Field | Type | Description | Required |
|---|---|---|---|
| `lead` | string | Lead variant as `chr_pos_ref_alt`, GRCh38 (the Open Targets variant id form). Its position centres the window. | Yes |
| `partners` | list of strings | Partner variants in the same form, same chromosome. The lead is removed from the list if present. An empty list returns no pairs and a note. | Yes |
| `chromosome` | string | Chromosome, with or without `chr`. | Yes |
| `window_bp` | integer | Total window width; the fetch covers `lead ± window_bp / 2`. Partners outside it get no row. | Yes |
| `super_pop` | string | `EUR`, `AFR`, `AMR`, `EAS` or `SAS`. Default `EUR`. | No |
| `plink_bin` | string | Path to the plink 1.9 binary; overrides `PLINK_BIN` and PATH. | No |

There is no `partners: null` mode: the client only computes r² for the ids you list.

## Output files

`<output>/ld_pairs.tsv`: four `#` header lines (skill version, panel, super-population, plink
version), then one row per pair:

```
lead_variant_id	partner_variant_id	r2	dprime	panel_id	super_pop
1_109274968_G_T	1_109272630_A_G	0.478979		1000g_phase3_v5b_grch38_basic	EUR
1_109274968_G_T	1_109274570_A_G	1.000000		1000g_phase3_v5b_grch38_basic	EUR
1_109274968_G_T	1_109274623_C_T	0.552982		1000g_phase3_v5b_grch38_basic	EUR
1_109274968_G_T	1_109274968_G_T	1.000000		1000g_phase3_v5b_grch38_basic	EUR
```

`dprime` is always empty (the plink call requests r² only). The last row is the lead
against itself.

`<output>/manifest.yaml` (or `manifest.json` when PyYAML is absent): `skill`, `version`,
`lead`, `chromosome`, `window_bp`, `super_pop`, `panel_id`, `panel_version`,
`plink_version`, `n_partners_requested`, `n_partners_returned`, `fetched_at_utc`,
`outputs`, and `notes` (the cached region path, plus one note per partner whose r² plink
reported as `NA`).

`<output>/report.md`: the same facts as a six-line Markdown list.

## Gotchas

1. **r² is only meaningful for an ancestry-matched panel.** EUR LD against an East Asian
   GWAS gives wrong blocks and a misleading plot. Choose `super_pop` from the study's stated
   ancestry, and say which panel was used (with its sample count) in every reply. Rules for
   mapping a study to a super-population, and the cases where none fits, are in
   [references/ancestry_matching.md](references/ancestry_matching.md).

2. **Alleles in a variant id must be in ASCII order, or the variant is not found.** plink
   names each panel variant with `--set-missing-var-ids '@:#:$1:$2'`, and `$1`/`$2` are the
   two alleles sorted (A < C < G < T), not REF then ALT. The client passes your ids verbatim.
   Measured 2026-09-13 with plink v1.90p on a synthetic VCF and on the live demo: a lead
   spelled `T_G` fails with `plink exited with code 5 ... No valid variants specified by
   --ld-snp`; a partner spelled `G_A` is silently absent; the same partner as `A_G` is
   returned. Sort the alleles in every id before calling (r² does not depend on which allele
   is called reference) and keep your own map back to the original spelling.

3. **`n_partners_returned` counts the lead's own row.** plink reports the lead against
   itself with r² = 1 and the parser keeps it as a pair. The demo shows `5 / 4`: three
   partners and the self-row. To count partners actually measured, drop the pair whose id
   equals the lead.

4. **A partner that is absent, monomorphic, or on the wrong build vanishes without a
   note.** Only the lead fails loudly (Gotcha 2's error also fires for a lead that is not in
   the panel). A partner not in 1000G, one monomorphic in the chosen super-population, or one
   given at a GRCh37 position (the panel is GRCh38 with no liftover) produces no row and no
   note. Compare returned ids against requested ids; do not read absence as r² = 0.

5. **There is no minor-allele-frequency filter.** Every extracted partner gets an r²,
   however rare. In a 503-sample super-population an allele at frequency 0.005 sits on about
   five haplotypes, so its r² rests on five observations. Apply a MAF floor downstream when
   the r² feeds a decision (instrument pruning) rather than a plot colour.

6. **plink 1.9's `--keep` wants the sample id twice.** With `--vcf`, plink 1.9 sets
   FID = IID = sample id ([plink 1.9 input docs](https://www.cog-genomics.org/plink/1.9/data#vcf)),
   so the keep file is `<sample>\t<sample>`. A plink 2 style `0\t<sample>` file keeps nobody
   and plink stops with "No people remaining after --keep". Do not swap in plink 2.

7. **The caches do not know the panel version.** The region VCF is cached by
   `(chromosome, start, end)` and the CLI result by `(lead, super_pop, window, partner set)`.
   The 2019-03-12 release has not moved, so a warm cache is exact today; if it ever moves,
   delete `~/.clawbio/locuscompare_cache/1000g/` by hand, because nothing will re-fetch. The
   cached region file is plain gzip, not BGZF: plink reads it, but a tabix or htslib reader
   that seeks will not.

8. **Admixed cohorts do not fit one super-population.** AMR, and the ASW and ACB samples
   inside AFR, carry the panel's own admixture, not the study's. When `super_pop` is AMR or
   the study is African American, say so in the reply; when a cohort matches none of the
   five groups, say that rather than picking one.

Details of the plink call and the parse, with the measurements behind Gotchas 2 to 5, are in
[references/r2_computation.md](references/r2_computation.md).

## Operational notes

- **Coordinates are GRCh38, 1-based**, matching Open Targets, GWAS Catalog harmonised files
  and the eQTL Catalogue. Contig names in the panel are bare (`1`); a `chr` prefix on input
  is stripped.
- **Network on first call per window.** The region VCF (the 2,548 samples of the 2019-03-12
  release, so one cached window serves every super-population; the keep files select from
  the 2,504 the 2013 panel file lists, see `references/1000g_panel.md`) and the 55 KB
  sample panel land in `~/.clawbio/locuscompare_cache/1000g/` (parent directory
  overridable with `LOCUSCOMPARE_CACHE_DIR`). The CLI additionally caches its result JSON in
  `~/.clawbio/ld_1000g_region_compute_cache/` (`LD_1000G_RESULT_CACHE_DIR`); `--no-cache`
  bypasses that result cache, not the region cache.
- **Timing.** With the region cached, the demo rerun with `--no-cache` took 0.77 s
  wall-clock on an Apple-silicon laptop (one measurement, 2026-09-13). The first call is
  dominated by the fetch.
- **Attribution.** Quote the panel and the tool in any methods text: 1000 Genomes Phase 3
  GRCh38 (Auton 2015; Clarke 2017), plink 1.9 (Chang 2015). A ready-made sentence is in
  [references/1000g_panel.md](references/1000g_panel.md).

## Safety

**Not for clinical decisions.** The output is a research-grade LD estimate from a public
reference panel, for visualisation and analysis planning; it is not a basis for diagnosis or
treatment.

**Reference-panel LD is not the study's LD.** The five super-populations are approximations
of any real cohort. For a population absent from 1000G, an admixed cohort, or rare variants,
treat the r² as a visual aid; for instrument pruning in Mendelian randomisation, prefer the
GWAS's own LD reference when one exists.

## Agent boundary

The skill computes r² between a lead and a partner list in a window from one 1000G
super-population. The agent should:

- **Use r² for plots and for LD pruning**, and name the threshold when calling something
  "in LD" (r² above 0.8, 0.6 or 0.2 mean different things; state which).
- **Expand the panel in the reply**: `LD reference: 1000G Phase 3 EUR (n = 503 samples)`,
  never a bare `EUR`, and quote the plink version from the manifest.
- **Flag an ancestry mismatch** between the study and `super_pop` explicitly, rather than
  presenting the r² as if it were the study's.
- **Read absence as absence**, not as r² = 0: check which requested partners came back
  (Gotchas 2 to 4), and say which did not and why when that can be determined.
- **Not silently sort ids on the user's behalf without saying so.** Sorting alleles is the
  right fix for Gotcha 2; the reply should still show the user's original spelling.

## References

Read only when the question needs them:

- [references/r2_computation.md](references/r2_computation.md): the exact plink command,
  the id rewrite and its ASCII-order consequence, what is dropped silently, the caches, and
  the timing measurement.
- [references/1000g_panel.md](references/1000g_panel.md): the two files fetched, measured
  super-population sizes, the shape of a fetched region, licence and a methods sentence.
- [references/ancestry_matching.md](references/ancestry_matching.md): mapping a study's
  stated ancestry to a super-population, the Finnish and admixed special cases, and when no
  1000G panel is the right choice.

## Sources & licensing

- 1000 Genomes Phase 3 GRCh38 phased release (2019-03-12), EBI 1000 Genomes FTP: open
  access with attribution. 1000 Genomes Project Consortium, Auton A, et al. A global
  reference for human genetic variation. Nature. 2015;526(7571):68-74.
  doi:10.1038/nature15393. PMID: 26432245. Clarke L, et al. The international Genome sample
  resource (IGSR). Nucleic Acids Res. 2017;45(D1):D854-D859. doi:10.1093/nar/gkw829.
  PMID: 27638885.
- plink 1.9: **GPL-3**, called as a separate process and never bundled or linked, so the
  skill's own code stays MIT. Chang CC, et al. Second-generation PLINK: rising to the
  challenge of larger and richer datasets. Gigascience. 2015;4:7.
  doi:10.1186/s13742-015-0047-8. PMID: 25722852.
- Skill code: **MIT**.
