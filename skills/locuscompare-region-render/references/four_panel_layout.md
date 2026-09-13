# The four-panel figure, panel by panel

What `render_full_locuscompare` in `scripts/regional_plot.py` draws, from which data, in the
order it appears on the page. Figure width is 13 in; height is the sum of the panel heights
(2.4 in per Manhattan track, 0.45 in per stacked gene row with a 0.9 in floor, 4.6 in for the
scatter row, plus title and caption); saved at 160 dpi with tight bounding box.

## Top stack (one shared genomic axis)

**1. Outcome Manhattan (GWAS).** `-log10 p` against position (Mb) for every variant of the
outcome slice, not only the joined ones, so the track fills the window the way a LocusZoom plot
does. Points are coloured by r² to the lead when the LD client returned one for that
`variant_id`; a variant with no r² is drawn in the grey `0.0 - 0.2` bin, which is also what
"not in LD with the lead" looks like. The lead is a purple diamond with its `variant_id`
annotated. Where the source tested no variants inside the requested window (a cis window
narrower than the plot, sparse coverage) the empty stretch is shaded light grey and labelled
`shaded = no source data tested`.

**2. Exposure Manhattan (QTL).** Same as panel 1 for the exposure slice. The panel titles are
the one-line labels built by the composer: `Outcome (GWAS): <trait label> (<accession>)` and
`Exposure (eQTL): <gene> | <study> | <sample group> | <quant method> (<dataset id>)` (or the
protein, panel and ancestry for a pQTL exposure).

**3. Gene track.** Genes in the window, packed into rows so overlapping spans do not collide.
Each gene is a line across its span with strand arrows every 4% of the window (pointing
right for `+`, left for `-`), exon boxes when the source provides exons (Ensembl REST does; a
synthetic TSV has none), and its symbol. The gene named by `exposure.gene_symbol`
is drawn bold and red. Only `gene_track.biotypes` (default `protein_coding`) are shown. When
the track is empty the row collapses and the exposure Manhattan carries the position axis.

The x-axis is forced to `lead_position +/- window_bp // 2`, so the lead is centred even when
the data end before the window edge.

## Bottom row (two scatters, joined pairs only)

Both scatters use the joined, allele-reconciled pairs with palindromic pairs removed, and are
downsampled past 5,000 points (below).

**4. LocusCompare panel.** x = `-log10 p` exposure, y = `-log10 p` outcome, dotted identity
line `y = x`, lead as the purple diamond. Titled `LocusCompare panel: clean diagonal supports
H4; two clusters supports H3`. Points along the diagonal are variants whose association
strength tracks across the two traits, which is what one shared causal variant plus its LD
partners produce.

**5. Effect-size panel.** x = `β` exposure, y = `β` outcome (after any allele flip), zero
lines, lead as the diamond, and a dashed line through the origin with slope
`β_outcome / β_exposure` at the lead, labelled `WR slope (lead) = <value>`. The line is an
orientation aid: it shows which direction the outcome moves per unit of exposure at the lead,
and whether the LD partners fall along the same line. No uncertainty is drawn.

## Colours

The LocusZoom five-bin r² palette, applied in every panel:

| r² | Colour | Hex |
|---|---|---|
| 0.8 to 1.0 | orange | `#FF7F0E` |
| 0.6 to 0.8 | green | `#2CA02C` |
| 0.4 to 0.6 | light blue | `#87CEEB` |
| 0.2 to 0.4 | navy | `#1F3A93` |
| 0.0 to 0.2, or no r² | grey | `#999999` |
| lead | purple diamond | `#9B30FF` |

Bins are drawn weakest first, so strong-LD points sit on top.

## Points that are not drawn

- A variant whose p-value is `None`, `0`, negative or non-finite has no `-log10 p` and is
  skipped in panels 1, 2 and 4 (it is not floored to a tiny value).
- A pair missing either beta is skipped in panel 5.
- Palindromic pairs are skipped in panels 4 and 5 only.

## Downsampling (panels 4 and 5)

When more than 5,000 joined non-palindromic pairs remain, every pair in the `0.8 - 1.0` bin
is kept and each other bin keeps a share of the remaining budget proportional to its size,
sampled with a fixed seed so the same input gives the same figure. `scatter_downsampled` and
`scatter_downsample_target` in the manifest record it. The Manhattan tracks are never
downsampled.

## Title and caption

Title: `Regional LocusCompare: <n joined> variants joined (+/-<kb> kb of <lead>)`, or the
`title` field of the render input. The monospace caption under the figure has six lines:
`LD:` (panel, super-population and plink version, or `no LD reference (plot rendered without
LD coloring)`), `Window:`, `Exposure:` (source, release and study id), `Outcome:`, the
provenance line (`OT release: ... | Rendered: <UTC timestamp>`), and `Caveats:` (the curated
caveats plus `<n> palindromic-ambiguous variants excluded from LD panel`, or `none`).
