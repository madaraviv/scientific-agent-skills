# Reading the LocusCompare panel

The bottom-left panel (`-log10 p` exposure against `-log10 p` outcome, coloured by r² to the
lead) is the diagnostic centre of the figure. This note says what its two canonical patterns
mean, what they do not prove, and what to quote next to them.

## The two patterns

**A clean diagonal.** Points run along the identity line, the lead sits at or near the top-right
corner, the orange and green (r² above 0.6) points cluster near it, and the grey points sit
near the origin. This is the shape one shared causal variant produces: every variant's
association strength in one trait is inherited from its LD to the same causal variant, so it
tracks the other trait. In the colocalisation vocabulary of Giambartolomei 2014 this is the
picture of hypothesis H4.

**Two clusters.** One cloud is strong in the exposure and weak in the outcome, the other the
reverse, and the high-r² points split between them. Two distinct causal variants in LD, each
driving one trait, produce it: hypothesis H3.

In between, a bent or fanned diagonal is common when one study is much better powered than the
other, or when the two cohorts' LD differs from the reference used for colouring.

## What the pattern does not prove

- A clean diagonal supports H4 but does not establish it. Two causal variants in very tight LD
  give a diagonal too; only the colocalisation posterior, computed from the full summary
  statistics, weighs that.
- Two clusters suggest H3 but a small, noisy study can scatter a true diagonal into apparent
  clusters.
- The lead variant need not be the most significant point in either track. A lead chosen as
  the most probable shared variant (for example from the product of two fine-mapping
  posteriors) is often not the minimum-p variant of either study. That is expected; a lead far
  from the top-right corner of the LocusCompare panel is what deserves a second look.
- LD colouring is from a reference panel, not from either cohort. Orange points far from the
  diagonal can mean the cohorts' LD differs from the reference (ancestry), or that the locus
  carries more than one signal.
- With only a handful of joined points the diagonal is undefined and clustering is not
  interpretable; read `n_pairs` in the report before reading the picture.

## What to quote next to the figure

1. The colocalisation posterior (H4, and H3 when the tool reports it), with its source and
   the threshold that source uses. This figure does not compute either.
2. The number of variants the tool colocalised over, and the widths of the two credible sets
   when available.
3. `n_pairs` and `n_palindromic_excluded` from `manifest.yaml`.
4. The LD panel and super-population, and whether they match the two cohorts' ancestry.
5. The sample sizes of the two studies.
6. Every caveat in the caption: grey points, the credible-set-filter note for non-gene-expression
   eQTL Catalogue datasets, a missing gene track, and the fetchers' own warnings.

## What not to say

- "The plot proves the colocalisation." It supports or weakens a posterior; it proves nothing.
- "The two clusters mean H3." Say that the pattern suggests distinct signals and that the
  posterior for H3 should be inspected.
- "The lead variant is causal." It is the variant the upstream tool chose to centre on; a
  mechanism needs functional evidence outside this skill.
- Anything about the dashed slope in the effect-size panel as an effect estimate. It is the
  ratio of the lead's two betas drawn for orientation, without a standard error.

## Sources

- Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation between
  pairs of genetic association studies using summary statistics. PLoS Genet.
  2014;10(5):e1004383. doi:10.1371/journal.pgen.1004383. PMID: 24830394. (H3 / H4.)
- Wallace C. A more accurate method for colocalisation analysis allowing for multiple causal
  variants. PLoS Genet. 2021;17(9):e1009440. doi:10.1371/journal.pgen.1009440.
  PMID: 34587156. (Multiple causal variants at one locus.)
- Liu B, Gloudemans MJ, Rao AS, Ingelsson E, Montgomery SB. Abundant associations with gene
  expression complicate GWAS follow-up. Nat Genet. 2019;51(5):768-769.
  doi:10.1038/s41588-019-0404-0. PMID: 31043754. (The LocusCompare figure.)
