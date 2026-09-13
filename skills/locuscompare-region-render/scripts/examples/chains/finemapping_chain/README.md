# Chain: from two fine-mapping runs

After fine-mapping each side of a candidate colocalisation, you want to see whether the two
signals line up. Fine-mapping answers which variants are credible for one signal; this figure
shows whether two signals share one causal architecture. A fine-mapping tool that writes the
harmonised slice it worked on (in the format of the skill's `references/input_schema.md`) can
feed both sides here directly, so nothing is re-fetched.

## Workflow

1. Fine-map the exposure locus and the outcome locus over the same GRCh38 window.
2. Point `exposure.sumstats_path` and `outcome.sumstats_path` at the two harmonised slices
   (the paths in `config.yaml` are placeholders relative to the config file).
3. Render:

```bash
python cli.py --input examples/chains/finemapping_chain/config.yaml --output runs/locuscompare_chain/
```

## What you should see

If the two credible sets colocalise: a diagonal in the LocusCompare panel and both Manhattan
peaks at the same position. If they do not: two clusters, or peaks at different positions
despite a high posterior; that is the case the visual check exists for.

## Notes

- Both slices must use the same window, the same build (GRCh38) and the same `variant_id`
  form; the join is on `variant_id`.
- The `fine_mapping_outputs:` block in `config.yaml` is documentation only; the CLI does not
  read it and posterior inclusion probabilities are not drawn.
