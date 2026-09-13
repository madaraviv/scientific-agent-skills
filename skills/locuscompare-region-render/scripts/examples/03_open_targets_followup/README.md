# 03: from an Open Targets colocalisation row (LDLR x fatty-acid composition)

You found this row in an Open Targets `colocalisation` query (`ot_row.json` holds it
verbatim):

```
left_studyId:  GCST90499945  (saturated fatty acids to total fatty acids %; Zoodsma 2025)
right_studyId: gtex_tx_skin_sun_exposed_enst00000559340
h4:            0.965
numberColocalisingVariants: 31
lead_variant:  19_11113815_A_G  (LDLR locus, chr19:11.1 Mb)
```

LDLR clears LDL particles from circulation and is the target of statins (indirectly) and of
PCSK9 inhibitors; a colocalisation between an LDLR expression QTL and a plasma
fatty-acid-composition GWAS is a plausible regulatory link worth looking at.

## Resolving the row into a config

| Open Targets field | Value | Config |
|---|---|---|
| `right_studyId` (the QTL side) | `gtex_tx_skin_sun_exposed_enst00000559340` | `exposure.fetch.dataset_id: QTD000318` (GTEx skin sun-exposed, transcript quantification); `molecular_trait_id: ENSG00000130164` (the LDLR gene id) |
| `left_studyId` (the GWAS side) | `GCST90499945` | `outcome.fetch.accession: GCST90499945` |
| `lead_variant` | `19_11113815_A_G` | `lead.variant_id` |

The QTL study id names a transcript; pass the parent gene id as `molecular_trait_id`, which
the eQTL Catalogue fetcher applies to its `gene_id` column so every LDLR transcript in the
window is returned. The `<study>_<quantification>_<tissue>` prefix maps to a dataset id through
the table bundled with the `eqtl-catalogue-region-fetch` skill. `provenance.ot_release` and
the `ot_*` keys are carried into the caption and manifest.

## Run

```bash
python cli.py --demo 03_open_targets_followup --output runs/03_open_targets_followup/
```

## What you should see

Four panels anchored at chr19:11113815 with LDLR bold red on the gene track. Because the
exposure is a transcript-level dataset the fetcher reads the credible-set-filtered file and
the manifest carries the caveat `sumstats are credible-set-filtered (eQTL Catalogue .cc.tsv.gz;
quant_method=tx)`: the exposure track is sparser than a gene-expression run.
