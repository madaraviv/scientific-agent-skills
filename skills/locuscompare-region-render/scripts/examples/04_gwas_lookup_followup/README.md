# 04: from a per-variant lookup

You looked up an rsID across GWAS and QTL databases, saw a strong GWAS and eQTL pair at it,
and want the regional view. A per-variant lookup says which studies the variant matters in; the
regional figure shows whether the two signals share one causal architecture across the locus.

## Translating a lookup result into the config

| Lookup output | Config field |
|---|---|
| the eQTL dataset (eQTL Catalogue `QTD...` id) | `exposure.fetch.dataset_id` |
| the eQTL gene id (`ENSG...`) | `exposure.fetch.molecular_trait_id` |
| the tissue label | `exposure.trait_label` (label only) |
| the GWAS accession (`GCST...`) | `outcome.fetch.accession` |
| the trait name | `outcome.trait_label` (label only) |
| the variant's GRCh38 chromosome, position and alleles | `lead.chromosome`, `lead.position_bp`, `lead.variant_id` (`<chr>_<pos>_<ref>_<alt>`) |

`provenance.gwas_lookup_run_dir` and `provenance.rsid_queried` are free-text provenance; the
run directory is printed in the caption.

## Run

```bash
python cli.py --demo 04_gwas_lookup_followup --output runs/locuscompare_rs12740374/
```

This config resolves to the same SORT1 x cholesterol-VLDL pair as demo 02 (QTD000276 x
GCST90269602, lead rs12740374), so the figure is the same; the difference is where the
identifiers came from. To compare several pairs from one lookup, write one config per pair.
