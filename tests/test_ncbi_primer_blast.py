"""Tests cho NCBIPrimerBlast — NCBI Primer-BLAST API interaction.

Tests use mocking to avoid real network calls.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import requests

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import DesignTarget, PrimerPair, ProcessedSequence
from tta_primer_design.modules.ncbi_primer_blast import (
    NCBIPrimerBlast,
    _extract_oligo_data,
)


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def api(config: AppConfig) -> NCBIPrimerBlast:
    return NCBIPrimerBlast(config)


@pytest.fixture
def target() -> DesignTarget:
    return DesignTarget(
        target_id="ACTB",
        accession="NM_001101",
        organism="Homo sapiens",
        design_mode="qpcr",
    )


@pytest.fixture
def processed_seq() -> ProcessedSequence:
    return ProcessedSequence(
        sequence="ATCGATCGATCGATCGATCG" * 10,  # 200 bp
        gc_content=0.5,
    )


# ---------------------------------------------------------------------------
# HTML fixtures
# ---------------------------------------------------------------------------

_HTML_WITH_JOB_KEY = """
<html>
<form>
<input type="hidden" name="JOB_KEY" value="abc123XYZ">
</form>
</html>
"""

_HTML_JOB_KEY_URL_STYLE = """
<html>
<a href="?CMD=Get&JOB_KEY=def456ABC">Check status</a>
</html>
"""

_HTML_RUNNING = """<html><body>Your job is still running. Please wait.</body></html>"""
_HTML_COMPLETE = """
<html>
<body>
Primer pairs found. Graphic view available.
<table>
  <tr>
    <td>Forward primer</td>
    <td>ATCGATCGATCGATCGATCG</td>
    <td>+</td>
    <td>20</td>
    <td>1</td>
    <td>20</td>
    <td>60.0</td>
    <td>50%</td>
  </tr>
  <tr>
    <td>Reverse primer</td>
    <td>GCTAGCTAGCTAGCTAGCTA</td>
    <td>-</td>
    <td>20</td>
    <td>180</td>
    <td>161</td>
    <td>59.5</td>
    <td>55%</td>
  </tr>
  <tr>
    <td>Product length</td>
    <td>180</td>
  </tr>
</table>
</body>
</html>
"""

_HTML_NO_PRIMERS = """<html><body>No primers found for the given parameters.</body></html>"""
_HTML_ERROR = """<html><body>Error: Internal Server Error occurred.</body></html>"""


# ---------------------------------------------------------------------------
# Unit tests — constants
# ---------------------------------------------------------------------------


class TestConstants:
    """Test class constants and initialization."""

    def test_base_url(self, api: NCBIPrimerBlast) -> None:
        assert api.BASE_URL == "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"

    def test_poll_interval(self, api: NCBIPrimerBlast) -> None:
        assert api.POLL_INTERVAL > 0

    def test_max_wait(self, api: NCBIPrimerBlast) -> None:
        assert api.MAX_WAIT > api.POLL_INTERVAL


# ---------------------------------------------------------------------------
# Unit tests — _extract_job_key
# ---------------------------------------------------------------------------


class TestExtractJobKey:
    """Test _extract_job_key static method."""

    def test_extract_hidden_input(self) -> None:
        key = NCBIPrimerBlast._extract_job_key(_HTML_WITH_JOB_KEY)
        assert key == "abc123XYZ"

    def test_extract_url_style(self) -> None:
        key = NCBIPrimerBlast._extract_job_key(_HTML_JOB_KEY_URL_STYLE)
        assert key == "def456ABC"

    def test_returns_none_if_missing(self) -> None:
        key = NCBIPrimerBlast._extract_job_key("<html><body>No job here</body></html>")
        assert key is None

    def test_extract_value_before_name(self) -> None:
        html = '<input value="xyz789" name="JOB_KEY" type="hidden">'
        key = NCBIPrimerBlast._extract_job_key(html)
        assert key == "xyz789"


# ---------------------------------------------------------------------------
# Unit tests — _classify_status
# ---------------------------------------------------------------------------


class TestClassifyStatus:
    """Test _classify_status static method."""

    def test_running(self) -> None:
        assert NCBIPrimerBlast._classify_status(_HTML_RUNNING) == "running"

    def test_complete(self) -> None:
        assert NCBIPrimerBlast._classify_status(_HTML_COMPLETE) == "complete"

    def test_no_primers(self) -> None:
        assert NCBIPrimerBlast._classify_status(_HTML_NO_PRIMERS) == "no_primers"

    def test_error(self) -> None:
        assert NCBIPrimerBlast._classify_status(_HTML_ERROR) == "error"

    def test_unknown_defaults_to_running(self) -> None:
        # Unknown HTML should default to "running" (safe assumption)
        assert NCBIPrimerBlast._classify_status("<html>Initializing...</html>") == "running"


# ---------------------------------------------------------------------------
# Unit tests — _extract_oligo_data
# ---------------------------------------------------------------------------


class TestExtractOligoData:
    """Test _extract_oligo_data helper."""

    def test_extract_sequence(self) -> None:
        cells = ["Forward primer", "ATCGATCGATCGATCGATCG", "+", "20", "1", "20", "60.0", "50%"]
        data = _extract_oligo_data(cells)
        assert data.get("sequence") == "ATCGATCGATCGATCGATCG"

    def test_extract_tm(self) -> None:
        cells = ["Forward primer", "ATCGATCGATCGATCGATCG", "+", "20", "1", "20", "60.0", "50%"]
        data = _extract_oligo_data(cells)
        assert data.get("tm") == "60.0"

    def test_extract_gc(self) -> None:
        cells = ["Forward primer", "ATCGATCGATCGATCGATCG", "+", "20", "1", "20", "60.0", "50%"]
        data = _extract_oligo_data(cells)
        assert data.get("gc") == "50"

    def test_no_sequence_in_cells(self) -> None:
        cells = ["Forward primer", "not-a-sequence", "60.0"]
        data = _extract_oligo_data(cells)
        assert "sequence" not in data


# ---------------------------------------------------------------------------
# Unit tests — parse_html_results
# ---------------------------------------------------------------------------


class TestParseHtmlResults:
    """Test parse_html_results."""

    def test_parse_complete_html(self, api: NCBIPrimerBlast) -> None:
        pairs = api.parse_html_results(_HTML_COMPLETE)
        assert isinstance(pairs, list)

    def test_parse_no_primers_html(self, api: NCBIPrimerBlast) -> None:
        pairs = api.parse_html_results(_HTML_NO_PRIMERS)
        assert pairs == []

    def test_parse_returns_primer_pairs(self, api: NCBIPrimerBlast) -> None:
        pairs = api.parse_html_results(_HTML_COMPLETE)
        for pair in pairs:
            assert isinstance(pair, PrimerPair)

    def test_parse_fallback_regex(self, api: NCBIPrimerBlast) -> None:
        html = """
        <html><body>
        Forward primer: ATCGATCGATCGATCGATCG
        Reverse primer: GCTAGCTAGCTAGCTAGCTA
        </body></html>
        """
        pairs = api.parse_html_results(html)
        assert isinstance(pairs, list)

    def test_parse_empty_html(self, api: NCBIPrimerBlast) -> None:
        pairs = api.parse_html_results("<html></html>")
        assert pairs == []


# ---------------------------------------------------------------------------
# Unit tests — _build_params
# ---------------------------------------------------------------------------


class TestBuildParams:
    """Test _build_params."""

    def test_includes_sequence(
        self, api: NCBIPrimerBlast, target: DesignTarget, processed_seq: ProcessedSequence
    ) -> None:
        params = api._build_params(target, processed_seq)
        assert "SEQUENCE_TEMPLATE" in params
        assert params["SEQUENCE_TEMPLATE"] == processed_seq.sequence

    def test_includes_organism(
        self, api: NCBIPrimerBlast, target: DesignTarget, processed_seq: ProcessedSequence
    ) -> None:
        params = api._build_params(target, processed_seq)
        assert "ORGANISM" in params
        # Homo sapiens → 9606
        assert params["ORGANISM"] == "9606"

    def test_custom_params_applied(
        self, api: NCBIPrimerBlast, processed_seq: ProcessedSequence
    ) -> None:
        target = DesignTarget(
            target_id="ACTB",
            custom_params={"PRIMER_OPT_SIZE": "22"},
        )
        params = api._build_params(target, processed_seq)
        assert params["PRIMER_OPT_SIZE"] == "22"

    def test_unknown_organism_passed_through(
        self, api: NCBIPrimerBlast, processed_seq: ProcessedSequence
    ) -> None:
        target = DesignTarget(target_id="X", organism="Danio rerio")
        params = api._build_params(target, processed_seq)
        assert params["ORGANISM"] == "Danio rerio"


# ---------------------------------------------------------------------------
# Integration tests (mocked network)
# ---------------------------------------------------------------------------


class TestSubmitJob:
    """Test submit_job with mocked HTTP."""

    def test_submit_returns_job_key(self, api: NCBIPrimerBlast) -> None:
        mock_resp = MagicMock()
        mock_resp.text = _HTML_WITH_JOB_KEY
        mock_resp.raise_for_status = MagicMock()

        with patch.object(api._session, "post", return_value=mock_resp):
            key = api.submit_job({"SEQUENCE_TEMPLATE": "ATCG"})
            assert key == "abc123XYZ"

    def test_submit_raises_on_missing_job_key(self, api: NCBIPrimerBlast) -> None:
        mock_resp = MagicMock()
        mock_resp.text = "<html>No job key here</html>"
        mock_resp.raise_for_status = MagicMock()

        with patch.object(api._session, "post", return_value=mock_resp):
            with pytest.raises(ValueError, match="JOB_KEY"):
                api.submit_job({"SEQUENCE_TEMPLATE": "ATCG"})

    def test_submit_raises_on_http_error(self, api: NCBIPrimerBlast) -> None:
        mock_resp = MagicMock()
        mock_resp.raise_for_status.side_effect = requests.HTTPError("500 Server Error")

        with patch.object(api._session, "post", return_value=mock_resp):
            with pytest.raises(requests.HTTPError):
                api.submit_job({"SEQUENCE_TEMPLATE": "ATCG"})


class TestPollStatus:
    """Test poll_status with mocked HTTP."""

    def test_poll_running(self, api: NCBIPrimerBlast) -> None:
        mock_resp = MagicMock()
        mock_resp.text = _HTML_RUNNING
        mock_resp.raise_for_status = MagicMock()

        with patch.object(api._session, "get", return_value=mock_resp):
            status = api.poll_status("abc123")
            assert status == "running"

    def test_poll_complete(self, api: NCBIPrimerBlast) -> None:
        mock_resp = MagicMock()
        mock_resp.text = _HTML_COMPLETE
        mock_resp.raise_for_status = MagicMock()

        with patch.object(api._session, "get", return_value=mock_resp):
            status = api.poll_status("abc123")
            assert status == "complete"

    def test_poll_error_on_exception(self, api: NCBIPrimerBlast) -> None:
        with patch.object(api._session, "get", side_effect=requests.ConnectionError("timeout")):
            status = api.poll_status("abc123")
            assert status == "error"


class TestDesignPrimers:
    """Test design_primers with mocked HTTP."""

    def test_design_primers_complete(
        self,
        api: NCBIPrimerBlast,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> None:
        with (
            patch.object(api, "submit_job", return_value="job_abc"),
            patch.object(api, "poll_status", return_value="complete"),
            patch.object(api, "fetch_results", return_value=_HTML_COMPLETE),
            patch.object(api, "parse_html_results", return_value=[]) as mock_parse,
        ):
            result = api.design_primers(target, processed_seq)
            assert isinstance(result, list)
            mock_parse.assert_called_once()

    def test_design_primers_no_primers(
        self,
        api: NCBIPrimerBlast,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> None:
        with (
            patch.object(api, "submit_job", return_value="job_xyz"),
            patch.object(api, "poll_status", return_value="no_primers"),
        ):
            result = api.design_primers(target, processed_seq)
            assert result == []

    def test_design_primers_error_raises(
        self,
        api: NCBIPrimerBlast,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> None:
        with (
            patch.object(api, "submit_job", return_value="job_err"),
            patch.object(api, "poll_status", return_value="error"),
        ):
            with pytest.raises(RuntimeError):
                api.design_primers(target, processed_seq)

    def test_design_primers_timeout(
        self,
        api: NCBIPrimerBlast,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> None:
        # Use a very small MAX_WAIT so the loop exits quickly
        api.MAX_WAIT = 5
        api.POLL_INTERVAL = 10  # Larger than MAX_WAIT → loop exits immediately

        with (
            patch.object(api, "submit_job", return_value="job_wait"),
            patch.object(api, "poll_status", return_value="running"),
            patch("time.sleep"),
        ):
            with pytest.raises(TimeoutError):
                api.design_primers(target, processed_seq)
