# r² computation: what the client asks plink 1.9 to do, and what falls out silently

The client (`scripts/ondemand_client.py`) computes r² between one lead and a partner list by
shelling out to plink 1.9 on a region VCF. Everything below was read from that code or
measured by running it; the measurements name their setup.

## The plink command

For a lead at position `P` on chromosome `C` with `window_bp = W`, the client fetches the
region `C:(P - W/2)-(P + W/2)` (every sample in the release, 2,548 columns) and runs:

```
plink --vcf <region.vcf.gz> \
      --keep <keep.txt> \
      --set-missing-var-ids '@:#:$1:$2' \
      --extract <extract.txt> \
      --r2 --ld-snp <lead as chr:pos:a1:a2> \
      --ld-window-kb <W / 1000> --ld-window 99999 --ld-window-r2 0 \
      --out <tmp>/ld_out
```

| Piece | What it does |
|---|---|
| `--keep keep.txt` | Restricts to the chosen super-population. Rows are `<sample>\t<sample>` because plink 1.9 with `--vcf` sets FID = IID = sample id ([plink 1.9 input docs](https://www.cog-genomics.org/plink/1.9/data#vcf)); a `0\t<sample>` row (the plink 2 convention) keeps nobody. |
| `--set-missing-var-ids '@:#:$1:$2'` | Gives every variant whose ID is `.` the name `chr:pos:allele1:allele2`. In the 20190312 release every ID is `.` (26,405 of 26,405 rows in the 1 Mb SORT1 demo window, checked 2026-09-13), so this names every variant. |
| `--extract extract.txt` | Keeps only the lead and the requested partners (ids with `_` swapped for `:`). Nothing else in the window is computed. |
| `--r2 --ld-snp <lead>` | r² between the lead and every remaining variant. |
| `--ld-window-kb`, `--ld-window 99999`, `--ld-window-r2 0` | No distance or r² floor inside the window: every extracted partner is reported, including r² = 0. |

plink writes `ld_out.ld` with whitespace-separated columns `CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2`;
`_parse_ld` keeps rows where `SNP_A` or `SNP_B` equals the lead, drops a row whose `R2` is
empty or `NA` with a note `missing r² for partner <id>`, and maps ids back to `chr_pos_a1_a2`.

## Allele order in variant ids

`$1` and `$2` are the two alleles in ASCII order, not REF then ALT. The client does not
reorder anything: it passes the caller's `chr_pos_ref_alt` with `_` replaced by `:`. So a
caller id matches the panel id only when its alleles already sort (A < C < G < T):
`1_109274968_G_T` matches, `1_109270398_G_A` never does.

Measured 2026-09-13 with plink v1.90p (6 Sep 2023) on a synthetic 5-variant, 40-sample VCF
routed through the client, and again on the live SORT1 demo:

| Case | Result |
|---|---|
| Lead `1_400_T_G` (VCF REF T, ALT G) | `OnDemandLDError: plink exited with code 5 ... No valid variants specified by --ld-snp` |
| Same variant requested as `1_400_G_T` | found; r² computed |
| Partner requested as `1_50_G_A` (VCF REF G, ALT A) | silently absent from the result |
| Same partner requested as `1_50_A_G` | returned as `1_50_A_G`, r² = 0.894 |
| Live demo, 5 partners, 2 of them out of order (`G_A`, `G_C`) | 3 partners returned, plus the self-row |

r² does not depend on which allele is called reference, so sorting the alleles in every id
before the call loses nothing; keep your own map from sorted id back to the original id,
because the result carries the sorted spelling.

## What is returned, and what is dropped without a note

- **The lead's own row is a pair.** plink reports the lead against itself (r² = 1) and the
  parser keeps it, so `pairs` contains the lead and `n_partners_returned` counts it. The demo
  reports `5 / 4`: three partners plus the self-row. To count partners actually measured, drop
  the pair whose id is the lead.
- **Lead absent from the panel (or spelled out of order): loud.** plink exits 5 with
  `No valid variants specified by --ld-snp` and the client raises `OnDemandLDError`. It does
  not return zeros.
- **Partner absent from the panel: silent.** It is not in the extract set, so it is not in
  the `.ld` file. Nothing is noted.
- **Partner monomorphic in the chosen super-population: silent.** Measured on the synthetic
  VCF: a partner with genotype `0|0` in every sample produced no row.
- **Partner on the wrong genome build: silent.** The panel is GRCh38 with no liftover; a
  GRCh37 position is simply absent and looks exactly like a panel gap.
- **No minor-allele-frequency filter.** Every extracted partner gets an r², however rare.
  In a 503-sample super-population an allele at frequency 0.005 is carried on about five
  haplotypes (503 × 2 × 0.005), so its r² is an estimate from five observations. Filter on
  MAF downstream if the use is inferential (instrument pruning) rather than a plot colour.

Because the drops are silent, compare the set of returned partner ids (minus the lead) with
the set requested; `n_partners_returned` alone cannot tell "absent" from "measured".

## Caches

Two caches, both under the home directory by default:

| What | Path | Key | Override |
|---|---|---|---|
| Region VCF (all samples) plus the sample panel TSV | `~/.clawbio/locuscompare_cache/1000g/chr<C>_<start>_<end>.vcf.gz` | `(chromosome, start, end)` | `LOCUSCOMPARE_CACHE_DIR` (parent of `1000g/`) |
| CLI result | `~/.clawbio/ld_1000g_region_compute_cache/<lead>__<super_pop>__win<W>__<sha1 of sorted partners>.json` | lead, super-pop, window, partner set | `LD_1000G_RESULT_CACHE_DIR`; `--no-cache` bypasses it |

The region file is written with Python's `gzip` (plain gzip, not BGZF; the 20-byte header
of the demo's cached file was checked). plink reads it; a tabix or htslib reader that needs
to seek will not. Neither key carries the panel version. The 20190312 release has not moved,
but if it ever does, a warm cache serves the old panel with no signal: delete the `1000g/`
directory on a panel change.

The region file holds the 2,548 samples of the 2019-03-12 release; the super-population
filter is applied by plink at compute time from the 2013 panel file's 2,504 samples, so one
cached window serves every super-population (45 VCF samples are in no keep file; the AFR keep
file names NA18498, which the VCF lacks, so AFR computes on 660; see
`references/1000g_panel.md`). Size for the 1 Mb SORT1 demo window on chromosome 1: 3.15 MB,
26,405 variants (1,779 indels, no multi-allelic rows).

## Timing

Single measurement, Apple-silicon laptop, 2026-09-13: the demo (1 Mb window, 5 partners,
EUR) rerun with the region already cached and `--no-cache` took 0.77 s wall-clock end to
end, including the plink subprocess. The first call is dominated by the fetch.

## Citations

- Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK:
  rising to the challenge of larger and richer datasets. Gigascience. 2015;4:7.
  doi:10.1186/s13742-015-0047-8. PMID: 25722852.
- plink 1.9 LD documentation: https://www.cog-genomics.org/plink/1.9/ld
