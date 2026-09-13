# Expected output: r² between the SORT1 lead and five nearby variants

Input: `input.json` (identical to `default.json`): lead `1_109274968_G_T` (rs12740374), five
partners within 5 kb, chromosome 1, a 1 Mb window, EUR.

## What you get

`OnDemandLDResult` with:

- `panel_id`: `1000g_phase3_v5b_grch38_basic`
- `panel_version`: `5b_remote_2019_03_12`
- `super_pop`: `EUR`
- `plink_version`: whatever `plink --version` prints, for example `PLINK v1.90p 64-bit (6 Sep 2023)`
- `n_partners_requested`: 5
- `n_partners_returned`: 4, which is three partners plus the lead's own row (r² = 1). The
  other two partners, `1_109270398_G_A` and `1_109274857_G_C`, have their alleles out of
  ASCII order and never match the panel id plink assigns; they are absent without a note.
  See the Gotchas in the skill's SKILL.md.
- `pairs`: `OnDemandLDPair(partner_variant_id, r2, dprime=None)` rows

Output of `bash run_example.sh` on 2026-09-13:

```
panel: 1000g_phase3_v5b_grch38_basic (EUR)
plink: PLINK v1.90p 64-bit (6 Sep 2023)
partners requested / returned: 5 / 4
pairs:
  1_109272630_A_G	r² = 0.479
  1_109274570_A_G	r² = 1.000
  1_109274623_C_T	r² = 0.553
  1_109274968_G_T	r² = 1.000
```

The values are fixed by the panel and the sample set; they change only if the 2019-03-12
release is replaced.

## Reproducing

```bash
bash run_example.sh
```

Network required on the first call (the demo window is a 3.15 MB region VCF plus a 55 KB
sample panel, both from the EBI 1000 Genomes FTP); later calls read the cache under
`~/.clawbio/locuscompare_cache/1000g/`. Requires plink 1.9 on PATH or `PLINK_BIN`.
