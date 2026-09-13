"""Live smoke test: the full matplotlib pipeline draws a gene track for SORT1 and
writes it to disk. Gated on RUN_LIVE_TESTS=1 (or `pytest -m live`). For a real
four-panel render from live sources use the bundled demo
`scripts/examples/02_eqtl_catalogue_x_gwas_catalog/` with the sibling skills
installed."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

from regional_plot import GeneTrackEntry, render_gene_track  # noqa: E402

pytestmark = pytest.mark.live


def test_live_render_gene_track_sort1_smoke(tmp_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sort1 = GeneTrackEntry(
        gene_symbol="SORT1", start=109_817_590, end=109_900_000, strand="+",
        exons=[(109_817_590, 109_817_838), (109_850_000, 109_851_000)], biotype="protein_coding",
    )
    fig, ax = plt.subplots(figsize=(8, 1.5))
    render_gene_track(ax, genes=[sort1], xlim_bp=(109_700_000, 109_900_000),
                      lead_position=109_817_590, track_label="Genes (GENCODE v39, synthetic test)")
    out = tmp_path / "sort1_gene_track.png"
    fig.savefig(out)
    plt.close(fig)
    assert out.is_file() and out.stat().st_size > 0
