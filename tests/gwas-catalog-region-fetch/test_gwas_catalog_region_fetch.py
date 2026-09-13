"""Unit tests for the GWAS Catalog harmonised region fetcher and its CLI.

`pysam.TabixFile` is mocked and the HTTP session is injected, so the suite runs offline
(`env -i PATH=/usr/bin:/bin HOME=$HOME python -m pytest tests/gwas-catalog-region-fetch`).
The live smoke test (a real tabix slice from the EBI FTP) is in
test_live_gwas_catalog_region_fetch.py.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest
import requests

# Inject the skill's scripts dir so the bare-name module resolves under the
# kebab-cased skill directory (`skills/gwas-catalog-region-fetch/`), which is not
# a valid Python identifier.
SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "gwas-catalog-region-fetch"
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

import skill_contract  # noqa: E402

import gwas_catalog_region_fetch as cli_module  # noqa: E402
from gwas_catalog_region_fetch import (  # noqa: E402
    GWASCatalogClient,
    GWASCatalogFetchError,
    GWASCatalogRelease,
    RegionResult,
    RegionVariant,
    _maybe_float,
    _maybe_int,
    gcst_url_base,
    harmonised_file_url,
)

# Every argparse script under scripts/ answers --help in a subprocess.
CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)


# ----------------------- URL construction -----------------------


def test_gcst_url_base_buckets_correctly():
    """GCST90475990 sits in bucket GCST90475001-GCST90476000."""
    base = gcst_url_base("GCST90475990")
    assert "GCST90475001-GCST90476000" in base
    assert base.endswith("/GCST90475990/harmonised")


def test_gcst_url_base_low_accession():
    """GCST000001 sits in bucket GCST000001-GCST001000."""
    base = gcst_url_base("GCST000001")
    assert "GCST000001-GCST001000" in base


def test_gcst_url_base_bucket_boundaries():
    """The bucket is 1,000 wide and starts at ...001: the last member of a bucket and the
    first member of the next land in different buckets."""
    assert "GCST90269001-GCST90270000" in gcst_url_base("GCST90270000")
    assert "GCST90270001-GCST90271000" in gcst_url_base("GCST90270001")


def test_gcst_url_base_strips_trailing_slash_on_ftp_base():
    base = gcst_url_base("GCST90269602", ftp_base="https://example.org/sumstats/")
    assert base.startswith("https://example.org/sumstats/GCST90269001-GCST90270000/")


@pytest.mark.parametrize("bad", [
    "../etc/passwd",
    "GCST90475990; rm -rf /",
    "GCST",
    "90475990",
    "gcst90475990",
    "GCST 90475990",
    "not_a_gcst",
])
def test_gcst_url_base_rejects_malformed_accession(bad):
    """Reject anything that doesn't match ^GCST\\d+$ before URL formatting."""
    with pytest.raises(ValueError):
        gcst_url_base(bad)


def test_harmonised_file_url_uses_simple_name_when_head_is_200(offline_ftp_session):
    url = harmonised_file_url("GCST90475990", session=offline_ftp_session)
    assert url.endswith("/GCST90475990/harmonised/GCST90475990.h.tsv.gz")


class _ListingSession:
    """HEAD answers `head_status`; GET answers the directory listing HTML."""

    def __init__(self, head_status: int, listing_html: str, get_status: int = 200):
        self.head_status = head_status
        self.listing_html = listing_html
        self.get_status = get_status
        self.calls: list[tuple[str, str]] = []

    def head(self, url, **kwargs):
        self.calls.append(("HEAD", url))
        return MagicMock(status_code=self.head_status)

    def get(self, url, **kwargs):
        self.calls.append(("GET", url))
        return MagicMock(status_code=self.get_status, text=self.listing_html)


def test_harmonised_file_url_falls_back_to_the_prefixed_name_from_the_listing():
    """The PMID-GCST-EFO spelling (GCST90019016 on the live FTP) is found by listing the
    directory when the simple name answers 404."""
    listing = (
        '<a href="34927100-GCST90019016-EFO_0000676-Build37.f.tsv.gz">raw</a>'
        '<a href="34927100-GCST90019016-EFO_0000676.h.tsv.gz">harm</a>'
        '<a href="34927100-GCST90019016-EFO_0000676.h.tsv.gz-meta.yaml">meta</a>'
    )
    sess = _ListingSession(head_status=404, listing_html=listing)
    url = harmonised_file_url("GCST90019016", session=sess)
    assert url.endswith("/GCST90019016/harmonised/34927100-GCST90019016-EFO_0000676.h.tsv.gz")
    assert [c[0] for c in sess.calls] == ["HEAD", "GET"]
    assert sess.calls[1][1].endswith("/harmonised/")


def test_harmonised_file_url_takes_the_first_listed_h_tsv_gz():
    """With several .h.tsv.gz links the first in listing order wins, silently."""
    listing = (
        '<a href="GCST90267286_buildGRCh37.h.tsv.gz">a</a>'
        '<a href="GCST90267286_other.h.tsv.gz">b</a>'
    )
    url = harmonised_file_url("GCST90267286", session=_ListingSession(404, listing))
    assert url.endswith("/GCST90267286_buildGRCh37.h.tsv.gz")


def test_harmonised_file_url_raises_when_listing_has_no_harmonised_file():
    listing = '<a href="md5sum.txt">md5</a><a href="GCST000001.tsv.gz">raw</a>'
    with pytest.raises(GWASCatalogFetchError, match="no harmonised .h.tsv.gz found"):
        harmonised_file_url("GCST000001", session=_ListingSession(404, listing))


def test_harmonised_file_url_raises_when_listing_is_not_200():
    with pytest.raises(GWASCatalogFetchError):
        harmonised_file_url("GCST000001", session=_ListingSession(404, "", get_status=503))


def test_harmonised_file_url_head_exception_falls_through_to_listing():
    """A transport error on HEAD is swallowed; the listing still resolves the name."""

    class _Sess(_ListingSession):
        def head(self, url, **kwargs):
            raise requests.ConnectionError("no route")

    sess = _Sess(404, '<a href="GCST000001.h.tsv.gz">h</a>')
    assert harmonised_file_url("GCST000001", session=sess).endswith("GCST000001.h.tsv.gz")


def test_harmonised_file_url_listing_exception_is_a_fetch_error():
    class _Sess(_ListingSession):
        def get(self, url, **kwargs):
            raise requests.ConnectionError("no route")

    with pytest.raises(GWASCatalogFetchError, match="could not list harmonised dir"):
        harmonised_file_url("GCST000001", session=_Sess(404, ""))


# ----------------------- TabixFile mocking -----------------------


def _mock_tabix(rows: list[str], header_cols: list[str] | None):
    """Build a MagicMock standing in for pysam.TabixFile.

    `header_cols=None` models a GWAS-SSF file, whose column line is not `#`-prefixed and
    so is invisible to tabix (`.header` empty).
    """
    tbx = MagicMock()
    tbx.header = [] if header_cols is None else ["\t".join(header_cols)]
    tbx.fetch.return_value = rows
    tbx.close = MagicMock()
    return tbx


HM_HEADER = [
    "hm_variant_id", "hm_rsid", "hm_chrom", "hm_pos",
    "hm_other_allele", "hm_effect_allele",
    "hm_beta", "hm_odds_ratio",
    "hm_effect_allele_frequency", "standard_error", "p_value",
]

# The GWAS-SSF layout as read from GCST90691573 on 2026-09-13: harmonised values in the
# standard columns, only hm_coordinate_conversion and hm_code prefixed.
SSF_HEADER = [
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid",
    "low_confidence", "hm_coordinate_conversion", "hm_code",
]


def test_fetch_region_parses_hm_layout_rows(monkeypatch, offline_ftp_session):
    rows = [
        "\t".join(["2_36910110_C_T", "rs10748691", "2", "36910110", "C", "T",
                   "-0.05", "NA", "0.32", "0.008", "1.2e-9"]),
        "\t".join(["2_36932656_A_G", "rs99999",   "2", "36932656", "A", "G",
                   "-0.04", "NA", "0.31", "0.008", "5.4e-9"]),
    ]
    tbx = _mock_tabix(rows, HM_HEADER)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)

    client = GWASCatalogClient(session=offline_ftp_session)
    result = client.fetch_region(
        accession="GCST90475990", chromosome="2",
        start_bp=36_410_000, end_bp=37_410_000,
    )
    assert result.accession == "GCST90475990"
    assert result.n_variants == 2
    v0 = result.variants[0]
    assert isinstance(v0, RegionVariant)
    assert v0.variant_id == "2_36910110_C_T"
    assert v0.chromosome == "2"
    assert v0.position == 36_910_110
    assert v0.ref == "C" and v0.alt == "T"
    assert v0.beta == pytest.approx(-0.05)
    assert v0.se == pytest.approx(0.008)
    assert v0.p_value == pytest.approx(1.2e-9)
    assert v0.effect_allele_frequency == pytest.approx(0.32)
    assert v0.odds_ratio is None
    assert v0.raw["hm_rsid"] == "rs10748691"
    tbx.close.assert_called_once()
    # 1-based inclusive input becomes a 0-based half-open tabix query.
    tbx.fetch.assert_called_once_with("2", 36_409_999, 37_410_000)


def _gzip_header_response(cols: list[str]):
    """A `Range` GET response carrying the gzipped first line of a GWAS-SSF file."""
    body = gzip.compress(("\t".join(cols) + "\n1\t2\tA\tG\n").encode())
    return MagicMock(status_code=206, content=body)


class _RangeSession:
    """HEAD 200 for name resolution; GET serves the gzipped file head for the column names."""

    def __init__(self, cols: list[str]):
        self.cols = cols
        self.get_calls: list[dict] = []

    def head(self, url, **kwargs):
        return MagicMock(status_code=200)

    def get(self, url, **kwargs):
        self.get_calls.append(kwargs)
        return _gzip_header_response(self.cols)


def test_fetch_region_parses_ssf_layout_via_http_range_header(monkeypatch):
    """A GWAS-SSF file has no tabix-visible header; the column names come from a 32 KB
    Range read, the standard columns are used, and the variant id is built."""
    rows = [
        "\t".join(["1", "108775337", "T", "C", "-0.00913978", "0.00497078",
                   "0.339325", "0.0659628", "rs1", "False", "rs", "10"]),
        "\t".join(["1", "108775456", "A", "T", "NA", "0.105754",
                   "0.00114873", "0.583193", "rs2", "False", "lo", "11"]),
    ]
    tbx = _mock_tabix(rows, header_cols=None)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)
    sess = _RangeSession(SSF_HEADER)

    result = GWASCatalogClient(session=sess).fetch_region("GCST90691573", "1", 108_774_968, 109_774_968)
    assert len(sess.get_calls) == 1
    assert sess.get_calls[0]["headers"] == {"Range": "bytes=0-32767"}
    assert result.n_variants == 2
    v0, v1 = result.variants
    assert v0.variant_id == "1_108775337_C_T"
    assert (v0.ref, v0.alt) == ("C", "T")
    assert v0.beta == pytest.approx(-0.00913978)
    assert v0.effect_allele_frequency == pytest.approx(0.339325)
    assert v0.raw["hm_code"] == "10"
    assert v1.beta is None          # "NA" -> None, row kept
    assert v1.variant_id == "1_108775456_T_A"


def test_fetch_region_raises_when_range_read_is_refused(monkeypatch):
    tbx = _mock_tabix([], header_cols=None)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)

    class _Sess:
        def head(self, url, **kwargs):
            return MagicMock(status_code=200)

        def get(self, url, **kwargs):
            return MagicMock(status_code=403, content=b"")

    with pytest.raises(GWASCatalogFetchError, match="could not fetch column-name row"):
        GWASCatalogClient(session=_Sess()).fetch_region("GCST000001", "1", 1, 1000)


def test_fetch_region_skips_rows_with_missing_essentials(monkeypatch, offline_ftp_session):
    rows = [
        "\t".join(["2_36910110_C_T", "rs1", "2", "36910110", "C", "T",
                   "-0.05", "NA", "0.32", "0.008", "1.2e-9"]),
        # missing chrom -> drop
        "\t".join(["2_X_C_T", "rs2", "", "36910200", "C", "T",
                   "0.01", "NA", "0.5", "0.005", "0.1"]),
    ]
    tbx = _mock_tabix(rows, HM_HEADER)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)

    client = GWASCatalogClient(session=offline_ftp_session)
    result = client.fetch_region("GCST90475990", "2", 36_410_000, 37_410_000)
    assert result.n_variants == 1
    assert result.notes == []


def test_fetch_region_skips_rows_with_wrong_column_count_and_notes_it(monkeypatch, offline_ftp_session):
    rows = [
        "\t".join(["2_36910110_C_T", "rs1", "2", "36910110", "C", "T",
                   "-0.05", "NA", "0.32", "0.008", "1.2e-9"]),
        "\t".join(["short", "row"]),
    ]
    tbx = _mock_tabix(rows, HM_HEADER)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)

    result = GWASCatalogClient(session=offline_ftp_session).fetch_region("GCST90475990", "2", 1, 10**8)
    assert result.n_variants == 1
    assert result.notes == ["row column count 2 != header 11; skipping"]


def test_fetch_region_handles_chr_prefix_retry(monkeypatch, offline_ftp_session):
    """Tabix indexes can have 'chr2' or '2'; client retries with the alt form."""
    rows = ["\t".join(["2_1_C_T", "rs1", "2", "1", "C", "T",
                       "-0.05", "NA", "0.32", "0.008", "1.2e-9"])]
    tbx = _mock_tabix(rows, HM_HEADER)

    calls: list[str] = []

    def flaky_fetch(chrom, *args, **kwargs):
        calls.append(chrom)
        if len(calls) == 1:
            raise ValueError("could not create iterator for region")
        return rows
    tbx.fetch = flaky_fetch

    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)
    client = GWASCatalogClient(session=offline_ftp_session)
    result = client.fetch_region("GCST90475990", "chr2", 1, 1000)
    assert calls == ["2", "chr2"]
    assert result.n_variants == 1


def test_fetch_region_absent_contig_propagates_value_error(monkeypatch, offline_ftp_session):
    """When both spellings are rejected the second ValueError is not wrapped: callers
    looping over chromosomes must catch ValueError as well as GWASCatalogFetchError."""
    tbx = _mock_tabix([], HM_HEADER)

    def reject(chrom, *args, **kwargs):
        raise ValueError(f"could not create iterator for region '{chrom}'")
    tbx.fetch = reject
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)

    with pytest.raises(ValueError):
        GWASCatalogClient(session=offline_ftp_session).fetch_region("GCST90475990", "26", 1, 1000)
    tbx.close.assert_called_once()


def test_fetch_region_empty_window_returns_zero_variants(monkeypatch, offline_ftp_session):
    tbx = _mock_tabix([], HM_HEADER)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)
    result = GWASCatalogClient(session=offline_ftp_session).fetch_region("GCST90475990", "2", 1, 1000)
    assert result.n_variants == 0
    assert result.variants == []


def test_fetch_region_records_release(monkeypatch, offline_ftp_session):
    tbx = _mock_tabix([], HM_HEADER)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)
    client = GWASCatalogClient(session=offline_ftp_session)
    result = client.fetch_region("GCST90475990", "2", 1, 1000)
    assert result.release.accession == "GCST90475990"
    assert result.release.harmonised_url.endswith("GCST90475990.h.tsv.gz")
    assert result.release.fetched_at_utc.endswith("Z")
    assert result.release.harmoniser_version is None
    assert (result.region_start_bp, result.region_end_bp) == (1, 1000)


def test_fetch_region_raises_on_unopenable_index(monkeypatch, offline_ftp_session):
    """The state EBI publishes for GCST90019016 (file served, no .tbi): pysam cannot open
    the index and the script raises its generic fetch error."""
    def boom(url):
        raise OSError("could not open index for `%s`" % url)
    monkeypatch.setattr("pysam.TabixFile", boom)

    client = GWASCatalogClient(session=offline_ftp_session)
    with pytest.raises(GWASCatalogFetchError, match="could not open tabix index for"):
        client.fetch_region("GCST99999999", "2", 1, 1000)


def test_fetch_region_raises_when_pysam_is_missing(monkeypatch, offline_ftp_session):
    import builtins
    real_import = builtins.__import__

    def no_pysam(name, *args, **kwargs):
        if name == "pysam":
            raise ImportError("No module named 'pysam'")
        return real_import(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", no_pysam)

    with pytest.raises(GWASCatalogFetchError, match="pysam is required"):
        GWASCatalogClient(session=offline_ftp_session).fetch_region("GCST000001", "1", 1, 10)


def test_to_dict_round_trips_through_the_cache_reader(monkeypatch, offline_ftp_session):
    rows = ["\t".join(["2_1_C_T", "rs1", "2", "1", "C", "T",
                       "-0.05", "1.2", "0.32", "0.008", "1.2e-9"])]
    tbx = _mock_tabix(rows, HM_HEADER)
    monkeypatch.setattr("pysam.TabixFile", lambda url: tbx)
    result = GWASCatalogClient(session=offline_ftp_session).fetch_region("GCST90475990", "2", 1, 1000)

    again = cli_module._region_result_from_cache(json.loads(json.dumps(result.to_dict())))
    assert again.to_dict() == result.to_dict()
    assert again.variants[0].odds_ratio == pytest.approx(1.2)


# ----------------------- scalar parsing -----------------------


@pytest.mark.parametrize("raw, expected", [
    ("0.5", 0.5), ("1e-9", 1e-9), ("NA", None), ("", None), (None, None),
    ("nan", None), ("abc", None),
])
def test_maybe_float(raw, expected):
    assert _maybe_float(raw) == expected


@pytest.mark.parametrize("raw, expected", [
    ("12", 12), ("12.0", 12), ("NA", None), ("", None), (None, None), ("x", None),
])
def test_maybe_int(raw, expected):
    assert _maybe_int(raw) == expected


# ----------------------- CLI -----------------------


def _canned_result(accession: str = "GCST90269602") -> RegionResult:
    release = GWASCatalogRelease(
        accession=accession,
        harmonised_url=f"https://example.org/{accession}.h.tsv.gz",
        fetched_at_utc="2026-09-13T08:27:47Z",
    )
    variants = [
        RegionVariant("1_108775337_C_T", "1", 108775337, "C", "T",
                      -0.00913978, 0.00497078, 0.0659628, None, 0.339325, {"hm_code": "10"}),
        RegionVariant("1_108775456_T_A", "1", 108775456, "T", "A",
                      None, 0.105754, 0.583193, None, None, {"hm_code": "11"}),
    ]
    return RegionResult(
        accession=accession, chromosome="1",
        region_start_bp=108774968, region_end_bp=109774968,
        n_variants=len(variants), variants=variants, release=release,
        notes=["row column count 2 != header 11; skipping"],
    )


@pytest.fixture
def fake_client(monkeypatch, tmp_path):
    """Replace the client the CLI builds and point the result cache at tmp_path."""
    calls: list[dict] = []

    class _Client:
        def __init__(self, *args, **kwargs):
            pass

        def fetch_region(self, **kwargs):
            calls.append(kwargs)
            return _canned_result(kwargs["accession"])

    monkeypatch.setattr(cli_module, "GWASCatalogClient", _Client)
    monkeypatch.setattr(cli_module, "DEFAULT_CACHE_DIR", tmp_path / "cache")
    return calls


def _read_manifest(out: Path) -> dict:
    y = out / "manifest.yaml"
    if y.is_file():
        yaml = pytest.importorskip("yaml")
        return yaml.safe_load(y.read_text())
    return json.loads((out / "manifest.json").read_text())


def test_cli_list_demos_names_the_bundled_configs(capsys):
    assert cli_module.main(["--list-demos"]) == 0
    out = capsys.readouterr().out
    for name in ("default", "sort1_cholesterol_vldl", "il6r_crp", "tcf7l2_hba1c", "input"):
        assert name in out
    assert "default (default)" in out


def test_cli_requires_output_and_a_config_source(tmp_path):
    with pytest.raises(SystemExit):
        cli_module.main([])
    with pytest.raises(SystemExit):
        cli_module.main(["--demo"])


def test_cli_writes_tsv_manifest_and_report(fake_client, tmp_path, capsys):
    cfg = tmp_path / "region.json"
    cfg.write_text(json.dumps({
        "accession": "GCST90269602", "chromosome": "1",
        "start_bp": 108774968, "end_bp": 109774968,
    }))
    out = tmp_path / "out"
    assert cli_module.main(["--input", str(cfg), "--output", str(out)]) == 0
    assert fake_client == [{"accession": "GCST90269602", "chromosome": "1",
                            "start_bp": 108774968, "end_bp": 109774968}]

    lines = (out / "variants.tsv").read_text().splitlines()
    assert lines[:3] == ["# locuscompare-schema-version: 1.0", "# source: gwas_catalog",
                         "# accession: GCST90269602"]
    assert lines[3].split("\t") == ["variant_id", "chromosome", "position_bp", "allele_a",
                                    "allele_b", "beta", "se", "p", "eaf", "study_id"]
    assert lines[4].split("\t") == ["1_108775337_C_T", "1", "108775337", "C", "T",
                                    "-0.00913978", "0.00497078", "0.0659628", "0.339325",
                                    "GCST90269602"]
    # None -> empty cell, not "None".
    assert lines[5].split("\t") == ["1_108775456_T_A", "1", "108775456", "T", "A",
                                    "", "0.105754", "0.583193", "", "GCST90269602"]
    assert len(lines) == 6

    manifest = _read_manifest(out)
    assert manifest["skill"] == "gwas-catalog-region-fetch"
    assert manifest["accession"] == "GCST90269602"
    assert manifest["region"] == {"chromosome": "1", "start_bp": 108774968, "end_bp": 109774968}
    assert manifest["n_variants"] == 2
    assert manifest["release"]["harmonised_url"] == "https://example.org/GCST90269602.h.tsv.gz"
    assert manifest["release"]["harmoniser_version"] is None
    assert manifest["outputs"] == {"variants_tsv": "variants.tsv"}
    assert "notes" not in manifest

    report = (out / "report.md").read_text()
    assert "`GCST90269602`" in report
    assert "chr1:108,774,968-109,774,968" in report
    assert "**Variants returned:** 2" in report

    stdout = capsys.readouterr().out
    assert "gwas-catalog-region-fetch: 2 variants ->" in stdout
    assert "accession GCST90269602" in stdout


def test_cli_accepts_yaml_config(fake_client, tmp_path):
    pytest.importorskip("yaml")
    cfg = tmp_path / "region.yaml"
    cfg.write_text('accession: GCST90691573\nchromosome: "1"\nstart_bp: 153925508\nend_bp: 154925508\n')
    assert cli_module.main(["--input", str(cfg), "--output", str(tmp_path / "out")]) == 0
    assert fake_client[0]["accession"] == "GCST90691573"
    assert _read_manifest(tmp_path / "out")["accession"] == "GCST90691573"


def test_cli_result_cache_hit_skips_the_client(fake_client, tmp_path):
    cfg = tmp_path / "region.json"
    cfg.write_text(json.dumps({"accession": "GCST90269602", "chromosome": "chr1",
                               "start_bp": 108774968, "end_bp": 109774968}))
    assert cli_module.main(["--input", str(cfg), "--output", str(tmp_path / "a")]) == 0
    cache_files = sorted(p.name for p in (tmp_path / "cache").iterdir())
    assert cache_files == ["GCST90269602__chr1_108774968_109774968.json"]

    assert cli_module.main(["--input", str(cfg), "--output", str(tmp_path / "b")]) == 0
    assert len(fake_client) == 1, "second run must be served from the cache"
    assert (tmp_path / "b" / "variants.tsv").read_text() == (tmp_path / "a" / "variants.tsv").read_text()

    assert cli_module.main(["--input", str(cfg), "--output", str(tmp_path / "c"), "--no-cache"]) == 0
    assert len(fake_client) == 2, "--no-cache must bypass the cache"


def test_cli_demo_resolves_the_bundled_default(fake_client, tmp_path, capsys):
    assert cli_module.main(["--demo", "--output", str(tmp_path / "out")]) == 0
    assert fake_client[0] == {"accession": "GCST90269602", "chromosome": "1",
                              "start_bp": 108774968, "end_bp": 109774968}
    assert "using bundled demo default.json" in capsys.readouterr().err


def test_cli_demo_by_name_reads_that_config(fake_client, tmp_path):
    assert cli_module.main(["--demo", "tcf7l2_hba1c", "--output", str(tmp_path / "out")]) == 0
    assert fake_client[0]["accession"] == "GCST90691576"
    assert fake_client[0]["chromosome"] == "10"


def test_cli_demo_unknown_name_lists_the_available_ones(fake_client, tmp_path):
    with pytest.raises(FileNotFoundError, match="Available: .*sort1_cholesterol_vldl"):
        cli_module.main(["--demo", "nope", "--output", str(tmp_path / "out")])


def test_cli_rejects_unknown_config_extension(fake_client, tmp_path):
    cfg = tmp_path / "region.toml"
    cfg.write_text("accession = 'GCST90269602'\n")
    with pytest.raises(ValueError, match="unsupported config extension"):
        cli_module.main(["--input", str(cfg), "--output", str(tmp_path / "out")])


def test_bundled_examples_are_the_documented_set():
    """The configs the SKILL.md table names exist, parse, and carry the four required keys."""
    examples = SKILL_ROOT / "scripts" / "examples"
    names = sorted(p.name for p in examples.iterdir())
    assert names == ["default.json", "expected_output.md", "il6r_crp.yaml", "input.json",
                     "run_example.sh", "sort1_cholesterol_vldl.json", "tcf7l2_hba1c.json"]
    for p in cli_module._list_demos():
        if p.suffix == ".yaml":
            pytest.importorskip("yaml")
        cfg = cli_module._load_config(p)
        assert {"accession", "chromosome", "start_bp", "end_bp"} <= set(cfg), p.name
        assert cli_module._GCST_RE.fullmatch(cfg["accession"]), p.name
