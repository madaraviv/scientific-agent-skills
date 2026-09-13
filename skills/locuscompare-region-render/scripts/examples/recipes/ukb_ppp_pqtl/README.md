# Recipe: UKB-PPP plasma pQTL

Convert a UKB-PPP plasma protein-QTL TSV into the canonical locuscompare
schema. **Registration-gated source: see Authentication.**

## Source

UKB-PPP: <https://www.synapse.org/Synapse:syn51365303> (Synapse-hosted)

- **Format**: tab-separated, per-protein cis + trans sumstats
- **Coordinate system**: GRCh38
- **License**: registration with UK Biobank required, plus a separate
  Synapse data-use agreement for UKB-PPP. The harmonised output carries a
  `# license: yellow ukb_ppp` header line so downstream tooling can see it.
- **Why this matters**: with the optional `ukb-ppp-region-fetch` skill absent
  from this collection, this recipe is the way to put a plasma pQTL on the
  exposure side (through `sumstats_path`).

## Authentication

1. Apply for UKB access at <https://www.ukbiobank.ac.uk/enable-your-research>
2. Apply for UKB-PPP-specific access at <https://www.synapse.org/Synapse:syn51365303>
3. Set up Synapse credentials: `pip install 'synapseclient>=4.0,<5' && synapse login`
4. Download a protein's pQTL TSV: `synapse get syn51365303` (browse subfolders)

## Native columns (cis-pQTL TSV)

```
CHROM  POS  ID  REF  ALT  A1  TEST  OBS_CT  BETA  SE  P  ERRCODE
```

## Quick start

```bash
# Once UKB-PPP access is approved:
synapse get syn51365303/CYS_C_OID20128_v1.tsv.gz

# Harmonise to canonical schema
bash examples/recipes/ukb_ppp_pqtl/harmonise.sh \
    CYS_C_OID20128_v1.tsv.gz \
    cys_c_pqtl.canonical.tsv

# Drop into config:
# exposure:
#   trait_label: "Cystatin C plasma protein (UKB-PPP OID20128)"
#   sumstats_path: "cys_c_pqtl.canonical.tsv.gz"
```

## Column mapping

| UKB-PPP column | locuscompare canonical |
|---|---|
| `CHROM` | `chromosome` |
| `POS` | `position_bp` |
| `REF` | `allele_a` |
| `ALT` | `allele_b` |
| `BETA` | `beta` |
| `SE` | `se` |
| `P` | `p` |
| `OBS_CT` | `n` |
| (synthesised) | `variant_id` = chr_pos_ref_alt |

## Caveats

- **A1 vs ALT**: PLINK2 distinguishes the "tested" allele (A1) from the
  reference (REF). For UKB-PPP, A1 is set to ALT, so the canonical `allele_b
  is effect allele` convention holds. Recipe asserts `A1 == ALT` per row;
  flips beta if not.
- **License visibility**: the harmoniser writes a `# license: yellow ukb_ppp`
  comment in the output's header; the renderer does not act on it, so keeping
  the file inside the terms of your access agreement is on you.
- **Multi-protein analyses**: each protein has its own TSV. For a multi-
  protein cis-coloc analysis, run this recipe per protein and treat each
  produced TSV as a separate exposure slice.

## Sign-check status

UKB-PPP standard release has been sign-checked against ALT as the effect
allele. Recipe asserts `A1 == ALT`; if not, the row is flagged.

## Citation

> Sun et al. (2023) *Plasma proteomic associations with genetics and health
> in the UK Biobank.* Nature 622, 329-338. doi:10.1038/s41586-023-06592-6
