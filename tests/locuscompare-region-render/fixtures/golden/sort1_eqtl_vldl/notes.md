# SORT1 eQTL x cholesterol-VLDL golden fixture

Two offline cassettes (`inputs/*.json.gz`) captured 2026-05-15: the eQTL
Catalogue QTD000276 slice (GTEx minor salivary gland, gene expression,
ENSG00000134243) and the GWAS Catalog harmonised GCST90269602 slice
(cholesterol in medium VLDL), both +/-500 kb of the SORT1 lead
`1_109274968_G_T` (rs12740374; Musunuru 2010). The cassettes are replayed
through stand-in clients, so the test exercises the composer's join,
harmonisation, caveat and manifest logic on real-shaped data without network.

Locked: study identifiers, window, lead, `n_pairs` (2547) and
`n_palindromic_excluded` (333), the downsample state and the caveat list with
`ld_client=None`. Not locked: `fetched_at`, `plot_artifact`.

The upstream copy of this fixture also carried a UKB-PPP pQTL half. That half
needs the optional `ukb-ppp-region-fetch` skill, which this suite does not
require, so it is not shipped here.
