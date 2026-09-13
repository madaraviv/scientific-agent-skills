"""Unit tests for `_fetchers/gencode_ondemand.py` (gene track via Ensembl REST):
cache-root resolution honours LOCUSCOMPARE_CACHE_DIR at call time, a cache hit
short-circuits the HTTP request, and the REST payload is parsed into the Gene /
Exon shape the renderer consumes. The REST GET is monkeypatched; no network."""
from __future__ import annotations

import json
import sys
from pathlib import Path

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

from _fetchers import gencode_ondemand  # noqa: E402
from _fetchers.gencode_ondemand import (  # noqa: E402
    Exon,
    Gene,
    _default_cache_dir,
    fetch_region_genes_remote,
)


def test_default_cache_dir_falls_back_to_home_when_env_unset(monkeypatch):
    monkeypatch.delenv("LOCUSCOMPARE_CACHE_DIR", raising=False)
    assert _default_cache_dir() == Path.home() / ".clawbio" / "locuscompare_cache" / "gencode"


def test_default_cache_dir_honours_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("LOCUSCOMPARE_CACHE_DIR", str(tmp_path))
    assert _default_cache_dir() == tmp_path / "gencode"


def test_default_cache_dir_is_evaluated_at_call_time(monkeypatch, tmp_path):
    monkeypatch.setenv("LOCUSCOMPARE_CACHE_DIR", str(tmp_path / "first"))
    assert _default_cache_dir() == tmp_path / "first" / "gencode"
    monkeypatch.setenv("LOCUSCOMPARE_CACHE_DIR", str(tmp_path / "second"))
    assert _default_cache_dir() == tmp_path / "second" / "gencode"


def _stub_rest_payload() -> list[dict]:
    """Minimal Ensembl REST `overlap/region` shape: one gene + two exons."""
    return [
        {"feature_type": "gene", "gene_id": "ENSG00000134243", "external_name": "SORT1",
         "biotype": "protein_coding", "start": 109_200_000, "end": 109_300_000, "strand": -1,
         "canonical_transcript": "ENST00000256637.5"},
        {"feature_type": "gene", "gene_id": "ENSG00000999999", "external_name": "PSEUDOFAKE",
         "biotype": "lncRNA", "start": 109_210_000, "end": 109_220_000, "strand": 1,
         "canonical_transcript": "ENSTPSEUDO.1"},
        {"feature_type": "exon", "Parent": "ENST00000256637.5",
         "start": 109_210_000, "end": 109_211_000, "rank": 1},
        {"feature_type": "exon", "Parent": "ENST00000256637.5",
         "start": 109_280_000, "end": 109_281_500, "rank": 2},
    ]


class _FakeResp:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def test_fetch_region_genes_remote_caches_to_explicit_dir(tmp_path, monkeypatch):
    """First call writes the cache; the second call makes no HTTP request."""
    payload = _stub_rest_payload()
    calls = {"n": 0}

    def fake_get(url, params, headers, timeout):
        calls["n"] += 1
        return _FakeResp(payload)

    monkeypatch.setattr(gencode_ondemand.requests, "get", fake_get)

    genes, meta, notes = fetch_region_genes_remote(
        chromosome="1", start_bp=109_200_000, end_bp=109_300_000,
        biotypes=("protein_coding",), cache_dir=tmp_path,
    )
    assert calls["n"] == 1
    cache_file = tmp_path / "1_109200000_109300000.json"
    assert cache_file.is_file()
    assert json.loads(cache_file.read_text()) == payload
    assert meta["release_label"] == "Ensembl REST (GRCh38)"
    assert any("fetched from Ensembl REST" in n for n in notes)
    assert {g.gene_symbol for g in genes} == {"SORT1"}  # only protein_coding survives

    sort1 = genes[0]
    assert isinstance(sort1, Gene)
    assert (sort1.strand, sort1.start, sort1.end) == ("-", 109_200_000, 109_300_000)
    assert len(sort1.exons) == 2 and all(isinstance(e, Exon) for e in sort1.exons)
    assert {e.exon_number for e in sort1.exons} == {1, 2}

    genes2, _, notes2 = fetch_region_genes_remote(
        chromosome="1", start_bp=109_200_000, end_bp=109_300_000,
        biotypes=("protein_coding",), cache_dir=tmp_path,
    )
    assert calls["n"] == 1, "cache should suppress the second HTTP call"
    assert any("cache" in n for n in notes2)
    assert {g.gene_symbol for g in genes2} == {"SORT1"}


def test_fetch_region_genes_remote_uses_env_override(tmp_path, monkeypatch):
    monkeypatch.setattr(gencode_ondemand.requests, "get",
                        lambda *a, **kw: _FakeResp(_stub_rest_payload()))
    monkeypatch.setenv("LOCUSCOMPARE_CACHE_DIR", str(tmp_path))
    fetch_region_genes_remote(chromosome="chr1", start_bp=109_200_000, end_bp=109_300_000,
                              biotypes=("protein_coding",))
    assert (tmp_path / "gencode" / "1_109200000_109300000.json").is_file()


def test_fetch_region_genes_remote_filters_genes_outside_window(tmp_path, monkeypatch):
    payload = [
        {"feature_type": "gene", "gene_id": "ENSGINWINDOW", "external_name": "IN_WINDOW",
         "biotype": "protein_coding", "start": 500, "end": 1500, "strand": 1, "canonical_transcript": ""},
        {"feature_type": "gene", "gene_id": "ENSGOUTSIDE", "external_name": "OUT_OF_WINDOW",
         "biotype": "protein_coding", "start": 5000, "end": 6000, "strand": 1, "canonical_transcript": ""},
    ]
    monkeypatch.setattr(gencode_ondemand.requests, "get", lambda *a, **kw: _FakeResp(payload))
    genes, _, _ = fetch_region_genes_remote(chromosome="1", start_bp=1000, end_bp=2000,
                                            biotypes=("protein_coding",), cache_dir=tmp_path)
    assert {g.gene_symbol for g in genes} == {"IN_WINDOW"}
