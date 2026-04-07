"""Tests cho ncbi_primer_blast module — submit/poll/parse và NCBIPrimerBlast class."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget, Oligo, PrimerPair, ProcessedSequence
from tta_primer_design.modules.ncbi_primer_blast import (
    _DEFAULT_PARAMS,
    NCBIPrimerBlast,
    _detect_status,
    _extract_job_key,
    _parse_primer_text_block,
)

# ---------------------------------------------------------------------------
# HTML fixtures
# ---------------------------------------------------------------------------

_HTML_WITH_HIDDEN_INPUT = """
<html><body>
<form>
  <input type="hidden" name="JOB_KEY" value="ABC123XYZ" />
</form>
</body></html>
"""

_HTML_WITH_JOB_KEY_URL = """
<html><body>
<a href="?CMD=Get&JOB_KEY=XYZ987ABC">check status</a>
</body></html>
"""

_HTML_WITH_JOB_ID_TEXT = """
<html><body>
<p>Job id# JOB456DEF is being processed.</p>
</body></html>
"""

_HTML_NO_JOB_KEY = """
<html><body><p>Something went wrong.</p></body></html>
"""

_HTML_RUNNING = """
<html><body>
<p>Your job is still running. Please wait.</p>
</body></html>
"""

_HTML_COMPLETE = """
<html><body>
<pre>
Primer pair 1
  Sequence (5'->3')        Start  Stop  Len   Tm    GC%
Forward primer  GCAAGGAATGGTTTCAGAAATCCA   10    33   24  60.1  45.8
Reverse primer  CAGGACTCCATGTCGTCCA        200   218   19  59.8  57.9
Product length  209
</pre>
</body></html>
"""

_HTML_MULTIPLE_PAIRS = """
<html><body>
<pre>
Primer pair 1
  Sequence (5'->3')        Start  Stop  Len   Tm    GC%
Forward primer  GCAAGGAATGGTTTCAGAAATCCA   10    33   24  60.1  45.8
Reverse primer  CAGGACTCCATGTCGTCCA        200   218   19  59.8  57.9
Product length  209

Primer pair 2
  Sequence (5'->3')        Start  Stop  Len   Tm    GC%
Forward primer  ATCGATCGATCGATCGATCG       50    69   20  61.0  50.0
Reverse primer  TAGCTAGCTAGCTAGCTAGC       300   319   20  60.5  50.0
Product length  270
</pre>
</body></html>
"""

_HTML_NO_PRIMERS = """
<html><body>
<p>No primers were found for the given sequence.</p>
</body></html>
"""

_HTML_ERROR = """
<html><body>
<p>An error occurred while processing your request. Invalid sequence.</p>
</body></html>
"""

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _make_target(organism: str = "Homo sapiens") -> DesignTarget:
    return DesignTarget(
        target_id="test_target",
        input_type="sequence",
        organism=organism,
    )


def _make_processed_seq(seq: str = "ATCGATCGATCGATCG" * 10) -> ProcessedSequence:
    return ProcessedSequence(sequence=seq)


def _make_pair(pair_id: str = "pair_01") -> PrimerPair:
    return PrimerPair(
        pair_id=pair_id,
        left_primer=Oligo(sequence="GCAAGGAATGGTTTCAGAAATCCA"),
        right_primer=Oligo(sequence="CAGGACTCCATGTCGTCCA"),
    )


# ---------------------------------------------------------------------------
# _extract_job_key
# ---------------------------------------------------------------------------


class TestExtractJobKey:
    """Tests for _extract_job_key()."""

    def test_extracts_from_hidden_input(self) -> None:
        key = _extract_job_key(_HTML_WITH_HIDDEN_INPUT)
        assert key == "ABC123XYZ"

    def test_extracts_from_url_pattern(self) -> None:
        key = _extract_job_key(_HTML_WITH_JOB_KEY_URL)
        assert key == "XYZ987ABC"

    def test_extracts_from_job_id_text(self) -> None:
        key = _extract_job_key(_HTML_WITH_JOB_ID_TEXT)
        assert key == "JOB456DEF"

    def test_raises_when_no_key_found(self) -> None:
        with pytest.raises(ValueError, match="Không tìm thấy job key"):
            _extract_job_key(_HTML_NO_JOB_KEY)

    def test_hidden_input_takes_priority_over_url(self) -> None:
        html = _HTML_WITH_HIDDEN_INPUT + _HTML_WITH_JOB_KEY_URL
        key = _extract_job_key(html)
        assert key == "ABC123XYZ"


# ---------------------------------------------------------------------------
# _detect_status
# ---------------------------------------------------------------------------


class TestDetectStatus:
    """Tests for _detect_status()."""

    def test_running_when_still_running_message(self) -> None:
        assert _detect_status(_HTML_RUNNING) == "running"

    def test_complete_when_primer_pair_present(self) -> None:
        assert _detect_status(_HTML_COMPLETE) == "complete"

    def test_complete_when_no_primers_message(self) -> None:
        assert _detect_status(_HTML_NO_PRIMERS) == "complete"

    def test_error_when_error_message(self) -> None:
        assert _detect_status(_HTML_ERROR) == "error"

    def test_running_for_unknown_html(self) -> None:
        assert _detect_status("<html><body></body></html>") == "running"

    def test_error_takes_priority_over_complete(self) -> None:
        html = _HTML_ERROR + _HTML_COMPLETE
        assert _detect_status(html) == "error"


# ---------------------------------------------------------------------------
# _parse_primer_text_block
# ---------------------------------------------------------------------------


class TestParsePrimerTextBlock:
    """Tests for _parse_primer_text_block()."""

    def test_parses_forward_and_reverse(self) -> None:
        text = """
Forward primer  GCAAGGAATGGTTTCAGAAATCCA   10    33   24  60.1  45.8
Reverse primer  CAGGACTCCATGTCGTCCA        200   218   19  59.8  57.9
Product length  209
"""
        pair = _parse_primer_text_block(text, 0)
        assert pair is not None
        assert pair.left_primer.sequence == "GCAAGGAATGGTTTCAGAAATCCA"
        assert pair.right_primer.sequence == "CAGGACTCCATGTCGTCCA"

    def test_sets_amplicon_size(self) -> None:
        text = """
Forward primer  ATCGATCGATCG   10    21   12  60.0  50.0
Reverse primer  TAGCTAGCTAGC   150   161   12  60.0  50.0
Product length  152
"""
        pair = _parse_primer_text_block(text, 0)
        assert pair is not None
        assert pair.amplicon_size == 152

    def test_sets_pair_id_with_index(self) -> None:
        text = """
Forward primer  ATCGATCGATCG   10    21   12  60.0  50.0
Reverse primer  TAGCTAGCTAGC   150   161   12  60.0  50.0
Product length  152
"""
        pair = _parse_primer_text_block(text, 2)
        assert pair is not None
        assert pair.pair_id == "ncbi_pair_03"

    def test_converts_start_to_0based(self) -> None:
        text = """
Forward primer  ATCGATCGATCG   10    21   12  60.0  50.0
Reverse primer  TAGCTAGCTAGC   150   161   12  60.0  50.0
Product length  152
"""
        pair = _parse_primer_text_block(text, 0)
        assert pair is not None
        assert pair.left_primer.start == 9  # 10 - 1

    def test_returns_none_if_no_primer_data(self) -> None:
        pair = _parse_primer_text_block("Some random text without primer info", 0)
        assert pair is None

    def test_parses_tm_and_gc(self) -> None:
        text = """
Forward primer  GCAAGGAATGGT   10    21   12  60.1  45.8
Reverse primer  CAGGACTCCATG   150   161   12  59.8  57.9
Product length  152
"""
        pair = _parse_primer_text_block(text, 0)
        assert pair is not None
        assert pair.left_primer.tm == pytest.approx(60.1)
        assert pair.left_primer.gc_percent == pytest.approx(45.8)
        assert pair.right_primer.tm == pytest.approx(59.8)

    def test_parses_left_right_keyword(self) -> None:
        text = """
Left primer     ATCGATCGATCG   10    21   12  60.0  50.0
Right primer    TAGCTAGCTAGC   150   161   12  60.0  50.0
Product length  150
"""
        pair = _parse_primer_text_block(text, 0)
        assert pair is not None
        assert pair.left_primer.sequence == "ATCGATCGATCG"


# ---------------------------------------------------------------------------
# NCBIPrimerBlast.parse_html_results
# ---------------------------------------------------------------------------


class TestParseHtmlResults:
    """Tests for NCBIPrimerBlast.parse_html_results()."""

    def test_returns_empty_when_no_primers(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        pairs = api.parse_html_results(_HTML_NO_PRIMERS)
        assert pairs == []

    def test_parses_single_primer_pair(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        pairs = api.parse_html_results(_HTML_COMPLETE)
        assert len(pairs) == 1
        assert pairs[0].left_primer.sequence == "GCAAGGAATGGTTTCAGAAATCCA"
        assert pairs[0].right_primer.sequence == "CAGGACTCCATGTCGTCCA"
        assert pairs[0].amplicon_size == 209

    def test_parses_multiple_primer_pairs(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        pairs = api.parse_html_results(_HTML_MULTIPLE_PAIRS)
        assert len(pairs) == 2
        assert pairs[0].pair_id == "ncbi_pair_01"
        assert pairs[1].pair_id == "ncbi_pair_02"

    def test_returns_empty_for_empty_html(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        pairs = api.parse_html_results("")
        assert pairs == []

    def test_returns_list_of_primer_pairs(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        pairs = api.parse_html_results(_HTML_COMPLETE)
        assert all(isinstance(p, PrimerPair) for p in pairs)


# ---------------------------------------------------------------------------
# NCBIPrimerBlast.submit_job
# ---------------------------------------------------------------------------


class TestSubmitJob:
    """Tests for NCBIPrimerBlast.submit_job() with mocked requests."""

    def test_posts_to_base_url(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_WITH_HIDDEN_INPUT
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.post", return_value=mock_resp
        ) as mock_post:
            key = api.submit_job({"SEQUENCE_TEMPLATE": "ATCG"})

        mock_post.assert_called_once()
        call_url = mock_post.call_args[0][0]
        assert call_url == NCBIPrimerBlast.BASE_URL
        assert key == "ABC123XYZ"

    def test_merges_custom_params_with_defaults(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_WITH_HIDDEN_INPUT
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.post", return_value=mock_resp
        ) as mock_post:
            api.submit_job({"ORGANISM": "Mus musculus", "SEQUENCE_TEMPLATE": "ATCG"})

        sent_data = mock_post.call_args[1]["data"]
        # Custom param merged
        assert sent_data["ORGANISM"] == "Mus musculus"
        # Default params present
        assert "PRIMER_OPT_SIZE" in sent_data

    def test_cmd_is_always_request(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_WITH_HIDDEN_INPUT
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.post", return_value=mock_resp
        ) as mock_post:
            api.submit_job({"CMD": "Get"})  # caller tries to override CMD

        sent_data = mock_post.call_args[1]["data"]
        assert sent_data["CMD"] == "request"

    def test_raises_on_http_error(self) -> None:
        from requests import HTTPError

        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.raise_for_status.side_effect = HTTPError("500 error")

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.post", return_value=mock_resp
        ):
            with pytest.raises(HTTPError):
                api.submit_job({"SEQUENCE_TEMPLATE": "ATCG"})


# ---------------------------------------------------------------------------
# NCBIPrimerBlast.poll_status
# ---------------------------------------------------------------------------


class TestPollStatus:
    """Tests for NCBIPrimerBlast.poll_status()."""

    def test_returns_running(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_RUNNING
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.get", return_value=mock_resp
        ):
            status = api.poll_status("JOBKEY123")

        assert status == "running"

    def test_returns_complete(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_COMPLETE
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.get", return_value=mock_resp
        ):
            status = api.poll_status("JOBKEY123")

        assert status == "complete"

    def test_returns_error(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_ERROR
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.get", return_value=mock_resp
        ):
            status = api.poll_status("JOBKEY123")

        assert status == "error"

    def test_passes_job_key_in_params(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_RUNNING
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.get", return_value=mock_resp
        ) as mock_get:
            api.poll_status("MY_JOB_KEY")

        params = mock_get.call_args[1]["params"]
        assert params["JOB_KEY"] == "MY_JOB_KEY"
        assert params["CMD"] == "Get"


# ---------------------------------------------------------------------------
# NCBIPrimerBlast.fetch_results
# ---------------------------------------------------------------------------


class TestFetchResults:
    """Tests for NCBIPrimerBlast.fetch_results()."""

    def test_returns_html_text(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = _HTML_COMPLETE
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.get", return_value=mock_resp
        ):
            html = api.fetch_results("KEY123")

        assert html == _HTML_COMPLETE

    def test_passes_correct_params(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        mock_resp = MagicMock()
        mock_resp.text = ""
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.ncbi_primer_blast.requests.get", return_value=mock_resp
        ) as mock_get:
            api.fetch_results("KEY_ABC")

        params = mock_get.call_args[1]["params"]
        assert params["CMD"] == "Get"
        assert params["JOB_KEY"] == "KEY_ABC"


# ---------------------------------------------------------------------------
# NCBIPrimerBlast._build_params
# ---------------------------------------------------------------------------


class TestBuildParams:
    """Tests for NCBIPrimerBlast._build_params()."""

    def test_includes_organism_and_sequence(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        target = _make_target(organism="Mus musculus")
        processed = _make_processed_seq("ATCGATCG")
        params = api._build_params(target, processed)
        assert params["ORGANISM"] == "Mus musculus"
        assert params["SEQUENCE_TEMPLATE"] == "ATCGATCG"

    def test_includes_region_when_set(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        target = _make_target()
        target.region_include = (100, 200)
        processed = _make_processed_seq()
        params = api._build_params(target, processed)
        assert "SEQUENCE_TARGET" in params

    def test_no_sequence_target_when_no_region(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        target = _make_target()
        processed = _make_processed_seq()
        params = api._build_params(target, processed)
        assert "SEQUENCE_TARGET" not in params

    def test_custom_params_override(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        target = _make_target()
        target.custom_params = {"PRIMER_NUM_RETURN": "5"}
        processed = _make_processed_seq()
        params = api._build_params(target, processed)
        assert params["PRIMER_NUM_RETURN"] == "5"


# ---------------------------------------------------------------------------
# NCBIPrimerBlast.design_primers (integration with mocks)
# ---------------------------------------------------------------------------


class TestDesignPrimers:
    """Tests for NCBIPrimerBlast.design_primers() with mocked sub-methods."""

    def test_returns_primer_pairs_on_success(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        api.MAX_WAIT = 100  # speed up test
        target = _make_target()
        processed = _make_processed_seq()

        with patch.object(api, "submit_job", return_value="JOB001"):
            with patch.object(api, "poll_status", return_value="complete"):
                with patch.object(api, "fetch_results", return_value=_HTML_COMPLETE):
                    pairs = api.design_primers(target, processed)

        assert len(pairs) == 1
        assert pairs[0].left_primer.sequence == "GCAAGGAATGGTTTCAGAAATCCA"

    def test_raises_runtime_error_on_job_error(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        target = _make_target()
        processed = _make_processed_seq()

        with patch.object(api, "submit_job", return_value="JOB_ERR"):
            with patch.object(api, "poll_status", return_value="error"):
                with pytest.raises(RuntimeError, match="lỗi"):
                    api.design_primers(target, processed)

    def test_raises_timeout_when_max_wait_exceeded(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        api.MAX_WAIT = 0  # immediately timeout
        target = _make_target()
        processed = _make_processed_seq()

        with patch.object(api, "submit_job", return_value="JOB_SLOW"):
            with patch.object(api, "poll_status", return_value="running"):
                with pytest.raises(TimeoutError):
                    api.design_primers(target, processed)

    def test_polls_until_complete(self) -> None:
        api = NCBIPrimerBlast(AppConfig())
        api.POLL_INTERVAL = 0  # no wait in test
        api.MAX_WAIT = 999
        target = _make_target()
        processed = _make_processed_seq()

        poll_responses = iter(["running", "running", "complete"])

        with patch.object(api, "submit_job", return_value="JOB002"):
            with patch.object(api, "poll_status", side_effect=poll_responses):
                with patch.object(api, "fetch_results", return_value=_HTML_COMPLETE):
                    with patch("tta_primer_design.modules.ncbi_primer_blast.time.sleep"):
                        pairs = api.design_primers(target, processed)

        assert len(pairs) == 1


# ---------------------------------------------------------------------------
# NCBIPrimerBlast constants
# ---------------------------------------------------------------------------


class TestConstants:
    def test_base_url_is_ncbi(self) -> None:
        assert NCBIPrimerBlast.BASE_URL.startswith("https://www.ncbi.nlm.nih.gov/")

    def test_poll_interval_positive(self) -> None:
        assert NCBIPrimerBlast.POLL_INTERVAL > 0

    def test_max_wait_positive(self) -> None:
        assert NCBIPrimerBlast.MAX_WAIT > 0

    def test_default_params_has_required_keys(self) -> None:
        assert "PRIMER_OPT_SIZE" in _DEFAULT_PARAMS
        assert "PRIMER_OPT_TM" in _DEFAULT_PARAMS
        assert "PRIMER_PRODUCT_SIZE_RANGE" in _DEFAULT_PARAMS
