"""Tests cho snp_checker module — helper functions và SNPChecker class."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tta_primer_design.config import AppConfig
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.snp_checker import (
    SNPChecker,
    _build_snp_result,
    _extract_alleles,
    _extract_maf,
    _fetch_snp_summaries,
    _search_snps_by_accession_range,
)

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _make_pair(pair_id: str = "pair_01", with_probe: bool = False) -> PrimerPair:
    pair = PrimerPair(
        pair_id=pair_id,
        left_primer=Oligo(sequence="GCAAGGAATGGTTTCAGAAATCCA", start=9, length=24),
        right_primer=Oligo(sequence="CAGGACTCCATGTCGTCCA", start=199, length=19),
    )
    if with_probe:
        pair.probe = Oligo(sequence="TGCAGCCACACTTTCTACAATGAGC", start=100, length=25)
    return pair


def _make_seq_record(accession: str = "NM_001101") -> object:
    """Minimal BioPython-like SeqRecord stub."""
    rec = MagicMock()
    rec.id = accession
    return rec


def _make_snp_summary(
    snp_id: str = "12345",
    maf: float = 0.05,
    allele_origin: str = "A/G",
    chrpos: str = "7:5600020",
) -> dict:
    return {
        "snp_id": snp_id,
        "uid": snp_id,
        "maf": str(maf),
        "allele_origin": allele_origin,
        "chrpos": chrpos,
    }


# ---------------------------------------------------------------------------
# _extract_maf
# ---------------------------------------------------------------------------


class TestExtractMaf:
    """Tests for _extract_maf()."""

    def test_returns_maf_from_direct_field(self) -> None:
        summary = {"maf": "0.15"}
        assert _extract_maf(summary) == pytest.approx(0.15)

    def test_returns_maf_as_float(self) -> None:
        summary = {"maf": 0.25}
        assert _extract_maf(summary) == pytest.approx(0.25)

    def test_returns_zero_when_no_maf(self) -> None:
        assert _extract_maf({}) == 0.0

    def test_returns_from_global_mafs(self) -> None:
        summary = {"global_mafs": [{"freq": "0.08"}]}
        assert _extract_maf(summary) == pytest.approx(0.08)

    def test_prefers_direct_maf_over_global_mafs(self) -> None:
        summary = {"maf": "0.12", "global_mafs": [{"freq": "0.08"}]}
        assert _extract_maf(summary) == pytest.approx(0.12)

    def test_returns_zero_for_invalid_maf_string(self) -> None:
        summary = {"maf": "n/a"}
        assert _extract_maf(summary) == 0.0

    def test_returns_zero_for_none_maf(self) -> None:
        summary = {"maf": None}
        assert _extract_maf(summary) == 0.0


# ---------------------------------------------------------------------------
# _extract_alleles
# ---------------------------------------------------------------------------


class TestExtractAlleles:
    """Tests for _extract_alleles()."""

    def test_returns_allele_origin(self) -> None:
        summary = {"allele_origin": "A/G"}
        assert _extract_alleles(summary) == "A/G"

    def test_returns_empty_when_missing(self) -> None:
        assert _extract_alleles({}) == ""

    def test_returns_docsum_as_fallback(self) -> None:
        summary = {"docsum": "C/T"}
        assert _extract_alleles(summary) == "C/T"


# ---------------------------------------------------------------------------
# _build_snp_result
# ---------------------------------------------------------------------------


class TestBuildSnpResult:
    """Tests for _build_snp_result()."""

    def test_pass_when_no_snps(self) -> None:
        result = _build_snp_result(
            oligo_id="test_oligo",
            oligo_length=24,
            snp_summaries=[],
            genomic_start=100,
            oligo_genomic_start=100,
            maf_threshold=0.01,
            is_probe=False,
        )
        assert result.recommendation == "PASS"
        assert result.has_snp is False

    def test_fail_when_snp_at_3prime_above_threshold(self) -> None:
        # SNP at position 22 in oligo of length 24 → 3' end (pos 22 >= 24-3=21)
        summary = _make_snp_summary(snp_id="111", maf=0.05, chrpos="7:122")  # 122 - 100 = 22
        result = _build_snp_result(
            oligo_id="test",
            oligo_length=24,
            snp_summaries=[summary],
            genomic_start=100,
            oligo_genomic_start=100,
            maf_threshold=0.01,
            is_probe=False,
        )
        assert result.recommendation == "FAIL"
        assert result.has_snp is True
        assert len(result.critical_snps) == 1

    def test_warning_when_snp_low_maf_not_3prime(self) -> None:
        # SNP at position 5, low MAF → WARNING (not critical enough for FAIL)
        summary = _make_snp_summary(snp_id="222", maf=0.005, chrpos="7:105")
        result = _build_snp_result(
            oligo_id="test",
            oligo_length=24,
            snp_summaries=[summary],
            genomic_start=100,
            oligo_genomic_start=100,
            maf_threshold=0.01,
            is_probe=False,
        )
        assert result.recommendation == "WARNING"
        assert result.has_snp is True
        assert result.critical_snps == []

    def test_fail_when_maf_above_threshold(self) -> None:
        # SNP in middle but MAF above threshold → FAIL (critical)
        summary = _make_snp_summary(snp_id="333", maf=0.05, chrpos="7:110")
        result = _build_snp_result(
            oligo_id="test",
            oligo_length=24,
            snp_summaries=[summary],
            genomic_start=100,
            oligo_genomic_start=100,
            maf_threshold=0.01,
            is_probe=False,
        )
        assert result.recommendation == "FAIL"
        assert len(result.critical_snps) == 1

    def test_probe_all_snps_are_critical(self) -> None:
        # For probe, any SNP → critical → FAIL
        summary = _make_snp_summary(snp_id="444", maf=0.001, chrpos="7:105")
        result = _build_snp_result(
            oligo_id="probe",
            oligo_length=24,
            snp_summaries=[summary],
            genomic_start=100,
            oligo_genomic_start=100,
            maf_threshold=0.01,
            is_probe=True,
        )
        assert result.recommendation == "FAIL"
        assert len(result.critical_snps) == 1

    def test_sets_rsid_correctly(self) -> None:
        summary = _make_snp_summary(snp_id="99999", maf=0.1)
        result = _build_snp_result(
            oligo_id="test",
            oligo_length=20,
            snp_summaries=[summary],
            genomic_start=100,
            oligo_genomic_start=100,
            maf_threshold=0.01,
            is_probe=False,
        )
        assert result.snp_list[0].rsid == "rs99999"

    def test_is_3prime_flagged_correctly(self) -> None:
        # oligo_length=20, 3' end = positions 17,18,19 (0-based)
        summary = _make_snp_summary(snp_id="555", maf=0.0, chrpos="7:118")  # 118-100=18 → 3'
        result = _build_snp_result(
            oligo_id="test",
            oligo_length=20,
            snp_summaries=[summary],
            genomic_start=100,
            oligo_genomic_start=100,
            maf_threshold=0.01,
            is_probe=False,
        )
        assert result.snp_list[0].is_3prime is True


# ---------------------------------------------------------------------------
# SNPChecker.check_oligo
# ---------------------------------------------------------------------------


class TestSNPCheckerCheckOligo:
    """Tests for SNPChecker.check_oligo()."""

    def test_returns_pass_when_no_accession(self) -> None:
        checker = SNPChecker(AppConfig())
        result = checker.check_oligo(
            oligo_sequence="ATCGATCG",
            genomic_start=100,
            genomic_end=107,
            accession="",
        )
        assert result.recommendation == "PASS"
        assert result.has_snp is False

    def test_returns_pass_when_no_snps_found(self) -> None:
        checker = SNPChecker(AppConfig())
        with patch(
            "tta_primer_design.modules.snp_checker._search_snps_by_accession_range",
            return_value=[],
        ):
            result = checker.check_oligo(
                oligo_sequence="ATCGATCG",
                genomic_start=100,
                genomic_end=107,
                accession="NM_001101",
            )
        assert result.recommendation == "PASS"

    def test_returns_snp_result_when_snps_found(self) -> None:
        checker = SNPChecker(AppConfig())
        summary = _make_snp_summary("11111", maf=0.02, chrpos="7:103")
        with patch(
            "tta_primer_design.modules.snp_checker._search_snps_by_accession_range",
            return_value=["11111"],
        ):
            with patch(
                "tta_primer_design.modules.snp_checker._fetch_snp_summaries",
                return_value=[summary],
            ):
                result = checker.check_oligo(
                    oligo_sequence="ATCGATCGATCGATCG",
                    genomic_start=100,
                    genomic_end=115,
                    accession="NM_001101",
                )
        assert result.has_snp is True

    def test_probe_mode_fails_on_any_snp(self) -> None:
        checker = SNPChecker(AppConfig())
        summary = _make_snp_summary("22222", maf=0.0, chrpos="7:102")
        with patch(
            "tta_primer_design.modules.snp_checker._search_snps_by_accession_range",
            return_value=["22222"],
        ):
            with patch(
                "tta_primer_design.modules.snp_checker._fetch_snp_summaries",
                return_value=[summary],
            ):
                result = checker.check_oligo(
                    oligo_sequence="ATCGATCGATCGATCGATCGATCG",
                    genomic_start=100,
                    genomic_end=123,
                    accession="NM_001101",
                    is_probe=True,
                )
        assert result.recommendation == "FAIL"

    def test_uses_config_email_in_api_call(self) -> None:
        cfg = AppConfig()
        cfg.ncbi.email = "test@lab.org"
        checker = SNPChecker(cfg)

        with patch(
            "tta_primer_design.modules.snp_checker._search_snps_by_accession_range",
            return_value=[],
        ) as mock_search:
            checker.check_oligo("ATCG", 100, 103, accession="NM_001")

        call_kwargs = mock_search.call_args[1]
        assert call_kwargs["email"] == "test@lab.org"


# ---------------------------------------------------------------------------
# SNPChecker.check_all
# ---------------------------------------------------------------------------


class TestSNPCheckerCheckAll:
    """Tests for SNPChecker.check_all()."""

    def test_returns_all_pairs(self) -> None:
        checker = SNPChecker(AppConfig())
        pairs = [_make_pair("p1"), _make_pair("p2")]
        seq_rec = _make_seq_record()

        with patch.object(checker, "_check_pair_snps"):
            result = checker.check_all(pairs, seq_rec)

        assert len(result) == 2

    def test_returns_pairs_unchanged_when_no_accession(self) -> None:
        checker = SNPChecker(AppConfig())
        pairs = [_make_pair()]
        result = checker.check_all(pairs, None)
        assert len(result) == 1
        assert result[0].snp_flags == []

    def test_handles_seq_record_without_id(self) -> None:
        checker = SNPChecker(AppConfig())
        pairs = [_make_pair()]
        seq_rec = MagicMock(spec=[])  # no .id attribute
        result = checker.check_all(pairs, seq_rec)
        assert len(result) == 1

    def test_catches_exception_per_pair(self) -> None:
        checker = SNPChecker(AppConfig())
        pair_ok = _make_pair("ok")
        pair_bad = _make_pair("bad")
        seq_rec = _make_seq_record()

        def side_effect(pair: PrimerPair, accession: str) -> None:
            if pair.pair_id == "bad":
                raise RuntimeError("API error")

        with patch.object(checker, "_check_pair_snps", side_effect=side_effect):
            result = checker.check_all([pair_ok, pair_bad], seq_rec)

        # Should not raise; both pairs returned
        assert len(result) == 2

    def test_empty_list_returns_empty(self) -> None:
        checker = SNPChecker(AppConfig())
        assert checker.check_all([], _make_seq_record()) == []

    def test_accession_version_suffix_stripped(self) -> None:
        """Accession 'NM_001101.3' should be stripped to 'NM_001101'."""
        checker = SNPChecker(AppConfig())
        seq_rec = _make_seq_record("NM_001101.3")
        pair = _make_pair()

        captured: list[str] = []

        def capture(p: PrimerPair, acc: str) -> None:
            captured.append(acc)

        with patch.object(checker, "_check_pair_snps", side_effect=capture):
            checker.check_all([pair], seq_rec)

        assert captured[0] == "NM_001101"

    def test_snp_flags_appended_to_pair(self) -> None:
        checker = SNPChecker(AppConfig())
        pair = _make_pair()
        seq_rec = _make_seq_record()

        # Simulate SNP found in left primer
        summary = _make_snp_summary("77777", maf=0.05, chrpos="7:18")

        with patch(
            "tta_primer_design.modules.snp_checker._search_snps_by_accession_range",
            return_value=["77777"],
        ):
            with patch(
                "tta_primer_design.modules.snp_checker._fetch_snp_summaries",
                return_value=[summary],
            ):
                checker.check_all([pair], seq_rec)

        assert len(pair.snp_flags) > 0
        assert any("LEFT" in f or "RIGHT" in f for f in pair.snp_flags)

    def test_probe_snp_flag_added_when_probe_present(self) -> None:
        checker = SNPChecker(AppConfig())
        pair = _make_pair(with_probe=True)
        seq_rec = _make_seq_record()

        # No SNPs for left/right, one SNP for probe
        probe_summary = _make_snp_summary("88888", maf=0.001, chrpos="7:110")

        call_count = 0

        def mock_search(accession: str, start: int, end: int, **kwargs: object) -> list[str]:
            nonlocal call_count
            call_count += 1
            # Only return SNP for probe region (start=101, end=125)
            if start == 101:
                return ["88888"]
            return []

        def mock_summary(ids: list[str], **kwargs: object) -> list[dict]:
            if "88888" in ids:
                return [probe_summary]
            return []

        with patch(
            "tta_primer_design.modules.snp_checker._search_snps_by_accession_range",
            side_effect=mock_search,
        ):
            with patch(
                "tta_primer_design.modules.snp_checker._fetch_snp_summaries",
                side_effect=mock_summary,
            ):
                checker.check_all([pair], seq_rec)

        assert any("PROBE" in f for f in pair.snp_flags)


# ---------------------------------------------------------------------------
# SNPChecker constructor
# ---------------------------------------------------------------------------


class TestSNPCheckerConstructor:
    def test_default_maf_threshold(self) -> None:
        checker = SNPChecker(AppConfig())
        assert checker.maf_threshold == 0.01

    def test_custom_maf_threshold(self) -> None:
        checker = SNPChecker(AppConfig(), maf_threshold=0.05)
        assert checker.maf_threshold == 0.05


# ---------------------------------------------------------------------------
# _search_snps_by_accession_range (mocked HTTP)
# ---------------------------------------------------------------------------


class TestSearchSNPs:
    """Tests for _search_snps_by_accession_range() with mocked requests."""

    def test_returns_id_list(self) -> None:
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"esearchresult": {"idlist": ["12345", "67890"]}}
        mock_resp.raise_for_status = MagicMock()

        with patch("tta_primer_design.modules.snp_checker.requests.get", return_value=mock_resp):
            ids = _search_snps_by_accession_range("NM_001101", 100, 200, "test@test.com")

        assert ids == ["12345", "67890"]

    def test_returns_empty_list_when_no_results(self) -> None:
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"esearchresult": {"idlist": []}}
        mock_resp.raise_for_status = MagicMock()

        with patch("tta_primer_design.modules.snp_checker.requests.get", return_value=mock_resp):
            ids = _search_snps_by_accession_range("NM_001101", 100, 200, "test@test.com")

        assert ids == []

    def test_includes_api_key_when_provided(self) -> None:
        mock_resp = MagicMock()
        mock_resp.json.return_value = {"esearchresult": {"idlist": []}}
        mock_resp.raise_for_status = MagicMock()

        with patch(
            "tta_primer_design.modules.snp_checker.requests.get", return_value=mock_resp
        ) as mock_get:
            _search_snps_by_accession_range("NM_001101", 100, 200, "e@test.com", api_key="mykey")

        params = mock_get.call_args[1]["params"]
        assert params["api_key"] == "mykey"


# ---------------------------------------------------------------------------
# _fetch_snp_summaries (mocked HTTP)
# ---------------------------------------------------------------------------


class TestFetchSnpSummaries:
    """Tests for _fetch_snp_summaries() with mocked requests."""

    def test_returns_empty_when_no_ids(self) -> None:
        summaries = _fetch_snp_summaries([], email="test@test.com")
        assert summaries == []

    def test_returns_summaries_list(self) -> None:
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "result": {
                "uids": ["111"],
                "111": {"snp_id": "111", "maf": "0.1"},
            }
        }
        mock_resp.raise_for_status = MagicMock()

        with patch("tta_primer_design.modules.snp_checker.requests.get", return_value=mock_resp):
            summaries = _fetch_snp_summaries(["111"], email="test@test.com")

        assert len(summaries) == 1
        assert summaries[0]["snp_id"] == "111"
