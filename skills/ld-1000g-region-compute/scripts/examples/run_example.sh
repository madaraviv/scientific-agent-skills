#!/usr/bin/env bash
# Example invocation: r² between the SORT1 lead (rs646776, 1_109274968_G_T) and
# the five partner variants in input.json, on the 1000G Phase 3 EUR panel.
# Network required on first call (~5-50 MB region VCF fetch from EBI 1000G FTP).
# Requires plink (1.9) on PATH (or PLINK_BIN env var pointing to it).

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python3 -c "
import json
from pathlib import Path
import sys

# Make the skill's modules importable: they live one level above examples/.
sys.path.insert(0, str(Path('${DIR}').parent))

from ondemand_client import OnDemand1000GLDClient

with open('${DIR}/input.json') as f:
    cfg = json.load(f)

client = OnDemand1000GLDClient(super_pop=cfg.get('super_pop', 'EUR'))
result = client.r2_with_lead(
    lead=cfg['lead'],
    partners=[p for p in cfg['partners'] if p != cfg['lead']],
    chromosome=cfg['chromosome'],
    window_bp=cfg['window_bp'],
)
print(f'panel: {result.panel_id} ({result.super_pop})')
print(f'plink: {result.plink_version}')
print(f'partners requested / returned: {result.n_partners_requested} / {result.n_partners_returned}')
print(f'pairs:')
for p in result.pairs:
    print(f'  {p.partner_variant_id}\\tr² = {p.r2:.3f}')
"
