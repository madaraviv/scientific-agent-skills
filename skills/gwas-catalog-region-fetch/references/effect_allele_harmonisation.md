# Effect-allele harmonisation across sources

A GWAS beta is the change in the trait per copy of the study's effect allele. Two sources can
report the same variant with opposite effect alleles: "T raises LDL by 0.05" and "C lowers LDL
by 0.05" are the same result with opposite signs. Comparing a GWAS beta with an eQTL or pQTL
beta, or with another GWAS, without first putting both on the same allele silently inverts
directions. This is the check to run before a Wald ratio, a colocalisation with effect
directions, or a meta-analysis.

## What this skill gives you

Each `variants.tsv` row carries `allele_a` (other allele) and `allele_b` (effect allele), and
`beta` is per copy of `allele_b`. The `variant_id` is `<chromosome>_<position>_<allele_a>_<allele_b>`
(or the file's `hm_variant_id` where one exists, which has the same shape). The GWAS Catalog
harmoniser has already oriented the alleles to the forward strand of the Ensembl reference and,
where it swapped them to match the reference pair, flipped the sign of beta, inverted the odds
ratio and taken `1 - EAF` (see `harmonised_files.md`). So within one harmonised file the
alleles are oriented to one reference; the script does not record which allele of the pair is
the reference's alternate, and a palindromic row can still be on either strand (`hm_code` 5 to
8 are assumptions, not inferences).

## The rule

For each variant shared between source A and source B:

1. Match on chromosome and position, then on the unordered allele pair. If the pairs differ
   (a multi-allelic site reported with different alternates, or an indel spelled differently),
   the variant does not match; do not force it.
2. If B's effect allele equals A's effect allele, keep B's beta as is.
3. If B's effect allele equals A's other allele, use `-beta` for B (and `1 / OR`, `1 - EAF`).
4. If the alleles are the reverse complement of A's (a strand difference), first
   reverse-complement B's alleles, then apply 2 or 3.
5. For palindromic variants (A/T, G/C) step 4 is ambiguous. Use the effect-allele frequencies
   from both sources to decide the strand when they are far from 0.5, and drop the variant when
   they are not; the usual practice in two-sample MR is to drop palindromes whose frequency is
   within a chosen band around 0.5 (TwoSampleMR's `harmonise_data` does this with a
   configurable threshold; Hemani 2018).

## Failure modes

- **Strand flip without allele flip.** A palindromic variant reported on opposite strands in two
  sources looks identical (`A/T` in both) and is opposite biology. Only frequency can tell, and
  only away from 0.5.
- **Multi-allelic sites.** A position with alleles A, C and G can be reported as A/C in one
  source and A/G in another. Matching on position alone assigns one alternate's effect to the
  other.
- **Indel spelling.** `A/AG` and `-/G` describe the same insertion; left-alignment and the
  padding base differ between sources. Do not hand-edit these; normalise both sides with the
  same tool or drop them.
- **Odds ratio versus beta.** A binary-trait file may carry an odds ratio and no beta. In the
  `hm_` layout `hm_odds_ratio` is on the effect allele; the natural log of it is the effect on
  the log-odds scale. This skill passes the odds ratio through on the Python object and leaves
  `beta` empty; it does not derive one.
- **A beta-less file.** One of 17 sampled files (GCST90399910) has `z_score` and no `beta`.
  Its rows have an empty `beta` cell; a z-score is not a beta and should not be used as one.

## Reference implementation

Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal inference
across the human phenome. eLife. 2018;7:e34408. doi:10.7554/eLife.34408. PMID: 29846171.
TwoSampleMR's `harmonise_data()` implements steps 1 to 5 above; a Python port needs the same
five decisions, in that order.
