"""Tests for `scripts/examples/01_synthetic_demo/generate_synthetic_fixtures.py`.

The generator is run from a copy in tmp_path (so the shipped fixtures stay
untouched) and must emit the four files with the schema the renderer consumes,
include the lead in both sumstats, be byte-stable across runs (seeded RNG), and
reproduce the shipped fixtures exactly."""
from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

import pytest

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "locuscompare-region-render"
DEMO_DIR = SKILL_ROOT / "scripts" / "examples" / "01_synthetic_demo"
GENERATOR_PATH = DEMO_DIR / "generate_synthetic_fixtures.py"
FIXTURE_NAMES = ("exposure.tsv", "outcome.tsv", "ld_matrix.tsv", "genes.tsv")

np = pytest.importorskip("numpy", reason="the fixture generator needs numpy")


def _load_generator_module(target_dir: Path):
    copy_path = target_dir / "generate_synthetic_fixtures.py"
    shutil.copyfile(GENERATOR_PATH, copy_path)
    spec = importlib.util.spec_from_file_location(f"_gen_test_{target_dir.name}", copy_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_generator_emits_four_files_with_expected_schema(tmp_path):
    module = _load_generator_module(tmp_path)
    module.generate()
    for name in FIXTURE_NAMES:
        assert (tmp_path / name).is_file(), f"generator missed {name}"
    expected_sumstats = "variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\tp"
    n = module.N_VARIANTS
    for name in ("exposure.tsv", "outcome.tsv"):
        lines = (tmp_path / name).read_text().splitlines()
        assert lines[0] == expected_sumstats
        assert len(lines) == n + 1
    ld_lines = (tmp_path / "ld_matrix.tsv").read_text().splitlines()
    assert ld_lines[0] == "partner_variant_id\tr2"
    assert len(ld_lines) == n  # header + (n - 1) partners; the lead is excluded
    gene_lines = (tmp_path / "genes.tsv").read_text().splitlines()
    assert gene_lines[0] == "gene_symbol\tstart\tend\tstrand\tbiotype"
    assert len(gene_lines) == 4
    assert any(line.startswith("DEMOGENE_B") for line in gene_lines[1:])


def test_generator_includes_the_lead_variant_in_both_sumstats(tmp_path):
    module = _load_generator_module(tmp_path)
    module.generate()
    lead_id = f"1_{module.LEAD_POSITION}_{module.LEAD_REF}_{module.LEAD_ALT}"
    for name in ("exposure.tsv", "outcome.tsv"):
        ids = [line.split("\t", 1)[0] for line in (tmp_path / name).read_text().splitlines()[1:]]
        assert lead_id in ids, f"{name} missing lead variant {lead_id}"


def test_generator_is_deterministic_and_matches_the_shipped_fixtures(tmp_path):
    first, second = tmp_path / "first", tmp_path / "second"
    first.mkdir()
    second.mkdir()
    _load_generator_module(first).generate()
    _load_generator_module(second).generate()
    for name in FIXTURE_NAMES:
        assert (first / name).read_bytes() == (second / name).read_bytes(), f"{name} differs across runs"
        assert (first / name).read_bytes() == (DEMO_DIR / name).read_bytes(), (
            f"{name} differs from the shipped fixture; re-run the generator or revert the edit"
        )
