"""Tests cho SNPChecker — SNP detection in primer/probe regions.

Tests use mocking to avoid real network calls.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.snp_checker import (
    SNPChecker,
    SNPInfo,
    SNPResult,
    _extract_genomic_pos,
)


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def checker(config: AppConfig) -> SNPChecker:
    return SNPChecker(config)


def _make_pair(
    pair_id: str = "pair_0",
    left_seq: str = "ATCGATCGATCGATCGATCG",
    right_seq: str = "GCTAGCTAGCTAGCTAGCTA",
    left_start: int = 100,
    right_start: int = 280,
) -> PrimerPair:
    left = Oligo(sequence=left_seq, start=left_start, tm=60.0)
    right = Oligo(sequence=right_seq, start=right_start, tm=60.0)
    return PrimerPair(pair_id=pair_id, left_primer=left, right_primer=right)


# Fake seq_record (BioPython-like)
class _FakeSeqRecord:
    def __init__(self, accession: str = "NM_001101") -> None:
        self.id = accession


# ---------------------------------------------------------------------------
# Unit tests — dataclasses
# ---------------------------------------------------------------------------


class TestDataclasses:
    """Test SNPInfo and SNPResult dataclasses."""

    def test_snp_info_defaults(self) -> None:
        snp = SNPInfo(rsid="rs12345")
        assert snp.maf == 0.0
        assert snp.is_3prime is False
        assert snp.position_in_oligo == 0

    def test_snp_result_defaults(self) -> None:
        result = SNPResult(oligo_id="primer_1")
        assert result.has_snp is False
        assert result.snp_list == []
        assert result.recommendation == "PASS"


# ---------------------------------------------------------------------------
# Unit tests — _determine_recommendation
# ---------------------------------------------------------------------------


class TestDetermineRecommendation:
    """Test _determine_recommendation logic."""

    def test_pass_no_snps(self, checker: SNPChecker) -> None:
        rec = checker._determine_recommendation([], is_probe=False)
        assert rec == "PASS"

    def test_warning_low_maf_primer(self, checker: SNPChecker) -> None:
        snp = SNPInfo(rsid="rs1", maf=0.005, is_3prime=False)
        rec = checker._determine_recommendation([snp], is_probe=False)
        assert rec == "WARNING"

    def test_fail_3prime_snp(self, checker: SNPChecker) -> None:
        snp = SNPInfo(rsid="rs2", maf=0.001, is_3prime=True)
        rec = checker._determine_recommendation([snp], is_probe=False)
        assert rec == "FAIL"

    def test_fail_high_maf_primer(self, checker: SNPChecker) -> None:
        snp = SNPInfo(rsid="rs3", maf=0.05, is_3prime=False)
        rec = checker._determine_recommendation([snp], is_probe=False)
        assert rec == "FAIL"

    def test_fail_any_maf_probe(self, checker: SNPChecker) -> None:
        snp = SNPInfo(rsid="rs4", maf=0.02, is_3prime=False)
        rec = checker._determine_recommendation([snp], is_probe=True)
        assert rec == "FAIL"

    def test_warning_low_maf_probe(self, checker: SNPChecker) -> None:
        snp = SNPInfo(rsid="rs5", maf=0.005, is_3prime=False)
        rec = checker._determine_recommendation([snp], is_probe=True)
        # MAF < threshold in probe → WARNING
        assert rec == "WARNING"


# ---------------------------------------------------------------------------
# Unit tests — check_oligo (mocked)
# ---------------------------------------------------------------------------


class TestCheckOligo:
    """Test check_oligo with mocked network calls."""

    def test_no_snps_returns_pass(self, checker: SNPChecker) -> None:
        checker._current_accession = "NM_001101"
        with patch.object(checker, "_search_snps_in_region", return_value=[]):
            result = checker.check_oligo("ATCGATCGATCGATCGATCG", 100, 120)
            assert result.recommendation == "PASS"
            assert result.has_snp is False

    def test_snp_found_returns_warning(self, checker: SNPChecker) -> None:
        checker._current_accession = "NM_001101"
        snp = SNPInfo(rsid="rs999", maf=0.005, position_in_oligo=5, is_3prime=False)
        with patch.object(checker, "_search_snps_in_region", return_value=["999"]):
            with patch.object(checker, "_fetch_snp_details", return_value=[snp]):
                result = checker.check_oligo("ATCGATCGATCGATCGATCG", 100, 120)
                assert result.has_snp is True
                assert result.recommendation == "WARNING"

    def test_snp_3prime_returns_fail(self, checker: SNPChecker) -> None:
        checker._current_accession = "NM_001101"
        snp = SNPInfo(rsid="rs888", maf=0.001, position_in_oligo=18, is_3prime=True)
        with patch.object(checker, "_search_snps_in_region", return_value=["888"]):
            with patch.object(checker, "_fetch_snp_details", return_value=[snp]):
                result = checker.check_oligo("ATCGATCGATCGATCGATCG", 100, 120)
                assert result.recommendation == "FAIL"
                assert snp in result.critical_snps

    def test_probe_snp_always_critical(self, checker: SNPChecker) -> None:
        checker._current_accession = "NM_001101"
        snp = SNPInfo(rsid="rs777", maf=0.05, position_in_oligo=5, is_3prime=False)
        with patch.object(checker, "_search_snps_in_region", return_value=["777"]):
            with patch.object(checker, "_fetch_snp_details", return_value=[snp]):
                result = checker.check_oligo("ATCGATCGATCGATCGATCG", 100, 120, is_probe=True)
                assert result.recommendation == "FAIL"
                assert snp in result.critical_snps

    def test_no_accession_returns_pass(self, checker: SNPChecker) -> None:
        checker._current_accession = None
        result = checker.check_oligo("ATCGATCGATCGATCGATCG", 100, 120)
        assert result.recommendation == "PASS"


# ---------------------------------------------------------------------------
# Unit tests — check_all
# ---------------------------------------------------------------------------


class TestCheckAll:
    """Test check_all with mocked network."""

    def test_check_all_sets_snp_flags(self, checker: SNPChecker) -> None:
        pairs = [_make_pair()]
        seq_record = _FakeSeqRecord("NM_001101")
        with patch.object(checker, "_search_snps_in_region", return_value=[]):
            result = checker.check_all(pairs, seq_record)
            assert result[0].snp_flags == []

    def test_check_all_sets_accession(self, checker: SNPChecker) -> None:
        pairs = [_make_pair()]
        seq_record = _FakeSeqRecord("NM_000999")
        with patch.object(checker, "_search_snps_in_region", return_value=[]):
            checker.check_all(pairs, seq_record)
            assert checker._current_accession == "NM_000999"

    def test_check_all_flags_fail(self, checker: SNPChecker) -> None:
        pairs = [_make_pair()]
        seq_record = _FakeSeqRecord("NM_001101")
        snp_fail = SNPInfo(rsid="rs100", maf=0.05, position_in_oligo=18, is_3prime=True)
        snp_result = SNPResult(
            oligo_id="test", has_snp=True, snp_list=[snp_fail], recommendation="FAIL"
        )
        with patch.object(checker, "check_oligo", return_value=snp_result):
            result = checker.check_all(pairs, seq_record)
            assert len(result[0].snp_flags) > 0
            assert any("FAIL" in f for f in result[0].snp_flags)

    def test_check_all_with_probe(self, checker: SNPChecker) -> None:
        pair = _make_pair()
        pair.probe = Oligo(sequence="CATGCATGCATGCATGCATG", start=130, tm=68.0)
        seq_record = _FakeSeqRecord("NM_001101")
        with patch.object(checker, "_search_snps_in_region", return_value=[]):
            result = checker.check_all([pair], seq_record)
            # probe was checked, no SNPs → no flags
            assert result[0].snp_flags == []

    def test_check_all_empty(self, checker: SNPChecker) -> None:
        result = checker.check_all([], _FakeSeqRecord())
        assert result == []

    def test_check_all_no_accession_in_record(self, checker: SNPChecker) -> None:
        pairs = [_make_pair()]
        seq_record = object()  # no .id attribute
        with patch.object(checker, "_search_snps_in_region", return_value=[]):
            result = checker.check_all(pairs, seq_record)
            assert checker._current_accession is None
            assert result[0].snp_flags == []


# ---------------------------------------------------------------------------
# Unit tests — _parse_refsnp_response
# ---------------------------------------------------------------------------


class TestParseRefSNPResponse:
    """Test _parse_refsnp_response."""

    def test_parse_basic_response(self, checker: SNPChecker) -> None:
        data = {
            "primary_snapshot_data": {
                "allele_annotations": [{"frequency": [{"minor_allele_freq": 0.05}]}],
                "placements_with_allele": [
                    {"alleles": [{"allele": {"spdi": {"inserted_sequence": "G", "position": 150}}}]}
                ],
            }
        }
        result = checker._parse_refsnp_response(data, "12345", 100, 20, is_probe=False)
        assert result is not None
        assert result.rsid == "rs12345"
        assert result.maf == 0.05

    def test_parse_3prime_detection(self, checker: SNPChecker) -> None:
        # SNP at position 118, oligo from 100 to 120 (length 20), so pos_in_oligo = 18 → 3'
        data = {
            "primary_snapshot_data": {
                "allele_annotations": [],
                "placements_with_allele": [
                    {"alleles": [{"allele": {"spdi": {"inserted_sequence": "A", "position": 118}}}]}
                ],
            }
        }
        result = checker._parse_refsnp_response(data, "5678", 100, 20, is_probe=False)
        assert result is not None
        assert result.is_3prime is True  # last 3 bases of 20-bp oligo

    def test_parse_empty_data(self, checker: SNPChecker) -> None:
        result = checker._parse_refsnp_response({}, "9999", 100, 20, is_probe=False)
        # Should not crash; might return None or an SNPInfo with defaults
        # (empty primary_snapshot_data)
        # Depending on implementation, either None or a default SNPInfo is acceptable
        assert result is None or isinstance(result, SNPInfo)


# ---------------------------------------------------------------------------
# Unit tests — _extract_genomic_pos
# ---------------------------------------------------------------------------


class TestExtractGenomicPos:
    """Test _extract_genomic_pos helper."""

    def test_extract_from_spdi(self) -> None:
        data: dict = {}
        placements = [
            {"alleles": [{"allele": {"spdi": {"position": 12345, "inserted_sequence": "A"}}}]}
        ]
        pos = _extract_genomic_pos(data, placements)
        assert pos == 12345

    def test_extract_empty(self) -> None:
        pos = _extract_genomic_pos({}, [])
        assert pos == 0

    def test_extract_missing_position(self) -> None:
        placements = [{"alleles": [{"allele": {"spdi": {"inserted_sequence": "A"}}}]}]
        pos = _extract_genomic_pos({}, placements)
        assert pos == 0
