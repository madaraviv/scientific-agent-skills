"""pytest configuration."""

import pytest


@pytest.fixture
def offline_ftp_session():
    """A `requests.Session` stand-in whose HEAD and GET both answer 200 with no body.

    `harmonised_file_url()` resolves the harmonised file name by HEAD-ing the simple
    `<GCST>.h.tsv.gz` form before any tabix call, so a test that only mocks
    `pysam.TabixFile` still reaches ftp.ebi.ac.uk. Injecting this session through the
    `session=` argument keeps the suite offline and makes the simple-name branch
    deterministic.
    """

    class _Resp:
        status_code = 200
        text = ""
        content = b""

    class _Sess:
        def head(self, url, **kwargs):
            return _Resp()

        def get(self, url, **kwargs):
            return _Resp()

    return _Sess()


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "live: integration test that reads the real EBI GWAS Catalog FTP via tabix",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    import os
    run_live = os.environ.get("RUN_LIVE_TESTS") == "1" or "live" in (config.getoption("-m") or "")
    if run_live:
        return
    skip_live = pytest.mark.skip(reason="live test (set RUN_LIVE_TESTS=1 or pytest -m live)")
    for item in items:
        if "live" in item.keywords:
            item.add_marker(skip_live)
