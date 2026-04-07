"""Tests cho BlastSpecificity — BLAST-based primer specificity checking.

Tests use mocking to avoid real network calls.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.blast_specificity import (
    BlastHit,
    BlastSpecificity,
    OffTargetAmplicon,
    SpecificityResult,
    _count_3prime_mismatches,
    _xml_text,
)


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def checker(config: AppConfig) -> BlastSpecificity:
    return BlastSpecificity(config)


def _make_pair(
    pair_id: str = "pair_0",
    left_seq: str = "ATCGATCGATCGATCGATCG",
    right_seq: str = "GCTAGCTAGCTAGCTAGCTA",
) -> PrimerPair:
    left = Oligo(sequence=left_seq, tm=60.0)
    right = Oligo(sequence=right_seq, tm=60.0)
    return PrimerPair(pair_id=pair_id, left_primer=left, right_primer=right)


# Minimal BLAST XML result (one hit, plus strand)
_BLAST_XML_ONE_HIT = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput//EN"
  "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>gi|12345|ref|NM_001101.5|</Hit_id>
          <Hit_accession>NM_001101</Hit_accession>
          <Hit_def>Homo sapiens actin beta (ACTB), mRNA</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_align-len>20</Hsp_align-len>
              <Hsp_mismatch>0</Hsp_mismatch>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>100</Hsp_hit-from>
              <Hsp_hit-to>119</Hsp_hit-to>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_bit-score>40.1</Hsp_bit-score>
              <Hsp_midline>||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

_BLAST_XML_EMPTY = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits></Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


# ---------------------------------------------------------------------------
# Unit tests — dataclasses
# ---------------------------------------------------------------------------


class TestDataclasses:
    """Test BlastHit, OffTargetAmplicon, SpecificityResult dataclasses."""

    def test_blast_hit_defaults(self) -> None:
        hit = BlastHit(subject_id="NM_001")
        assert hit.identity == 0.0
        assert hit.mismatches == 0
        assert hit.mismatches_3prime == 0

    def test_specificity_result_defaults(self) -> None:
        result = SpecificityResult(primer_pair_id="pair_0")
        assert result.is_specific is True
        assert result.specificity_score == 100.0
        assert result.off_target_amplicons == []
        assert result.blast_hits_probe is None

    def test_off_target_amplicon(self) -> None:
        lh = BlastHit(subject_id="NM_001", subject_start=100, subject_end=119)
        rh = BlastHit(subject_id="NM_001", subject_start=300, subject_end=281)
        amp = OffTargetAmplicon(subject_id="NM_001", amplicon_size=200, left_hit=lh, right_hit=rh)
        assert amp.amplicon_size == 200


# ---------------------------------------------------------------------------
# Unit tests — helper functions
# ---------------------------------------------------------------------------


class TestHelpers:
    """Test module-level helper functions."""

    def test_xml_text_found(self) -> None:
        import xml.etree.ElementTree as ET

        root = ET.fromstring("<Hit><Hit_id>NM_001</Hit_id></Hit>")
        assert _xml_text(root, "Hit_id") == "NM_001"

    def test_xml_text_missing(self) -> None:
        import xml.etree.ElementTree as ET

        root = ET.fromstring("<Hit></Hit>")
        assert _xml_text(root, "Hit_id", "default") == "default"

    def test_count_3prime_mismatches_no_mismatch(self) -> None:
        midline = "||||||||||||||||||||"
        assert _count_3prime_mismatches(midline, 100, 119) == 0

    def test_count_3prime_mismatches_plus_strand(self) -> None:
        # last 3 chars have spaces (mismatches)
        midline = "||||||||||||||||||   "
        # s_start < s_end → plus strand, 3' is right side
        assert _count_3prime_mismatches(midline, 100, 119) == 3

    def test_count_3prime_mismatches_minus_strand(self) -> None:
        # first 3 chars have spaces (mismatches at 3' end for minus strand)
        midline = "   ||||||||||||||||||"
        # s_start > s_end → minus strand, 3' is left side
        assert _count_3prime_mismatches(midline, 119, 100) == 3

    def test_count_3prime_empty_midline(self) -> None:
        assert _count_3prime_mismatches("", 100, 119) == 0


# ---------------------------------------------------------------------------
# Unit tests — _parse_blast_xml
# ---------------------------------------------------------------------------


class TestParseBLASTXML:
    """Test _parse_blast_xml."""

    def test_parse_one_hit(self, checker: BlastSpecificity) -> None:
        hits = checker._parse_blast_xml(_BLAST_XML_ONE_HIT, 20)
        assert len(hits) == 1
        assert hits[0].subject_id == "gi|12345|ref|NM_001101.5|"
        assert hits[0].alignment_length == 20
        assert hits[0].mismatches == 0
        assert hits[0].subject_start == 100
        assert hits[0].subject_end == 119

    def test_parse_empty(self, checker: BlastSpecificity) -> None:
        hits = checker._parse_blast_xml(_BLAST_XML_EMPTY, 20)
        assert hits == []

    def test_parse_invalid_xml(self, checker: BlastSpecificity) -> None:
        hits = checker._parse_blast_xml("not valid xml", 20)
        assert hits == []

    def test_parse_identity_percentage(self, checker: BlastSpecificity) -> None:
        hits = checker._parse_blast_xml(_BLAST_XML_ONE_HIT, 20)
        assert hits[0].identity == 100.0  # 20/20 * 100


# ---------------------------------------------------------------------------
# Unit tests — _find_off_target_amplicons
# ---------------------------------------------------------------------------


class TestFindOffTargetAmplicons:
    """Test _find_off_target_amplicons."""

    def test_no_off_target_same_subject(self, checker: BlastSpecificity) -> None:
        """When primers hit on same strand, no off-target amplicon."""
        left_hit = BlastHit(
            subject_id="NM_001", subject_start=100, subject_end=119, mismatches_3prime=0
        )
        right_hit = BlastHit(
            subject_id="NM_001", subject_start=300, subject_end=319, mismatches_3prime=0
        )  # both plus strand → no amplicon
        result = checker._find_off_target_amplicons([left_hit], [right_hit])
        assert result == []

    def test_off_target_detected(self, checker: BlastSpecificity) -> None:
        """Left hit plus strand + right hit minus strand at reasonable distance."""
        left_hit = BlastHit(
            subject_id="NM_002", subject_start=100, subject_end=119, mismatches_3prime=0
        )
        right_hit = BlastHit(
            subject_id="NM_002",
            subject_start=350,  # minus strand: start > end
            subject_end=331,
            mismatches_3prime=0,
        )
        result = checker._find_off_target_amplicons([left_hit], [right_hit])
        assert len(result) == 1
        assert result[0].amplicon_size == 350 - 100

    def test_off_target_too_large(self, checker: BlastSpecificity) -> None:
        """Off-target amplicon larger than max size should be ignored."""
        left_hit = BlastHit(
            subject_id="NM_003", subject_start=100, subject_end=119, mismatches_3prime=0
        )
        right_hit = BlastHit(
            subject_id="NM_003",
            subject_start=10000,  # too far away
            subject_end=9981,
            mismatches_3prime=0,
        )
        result = checker._find_off_target_amplicons([left_hit], [right_hit])
        assert result == []

    def test_3prime_mismatch_skipped(self, checker: BlastSpecificity) -> None:
        """Hits with 3' mismatches should be ignored (primer won't extend)."""
        left_hit = BlastHit(
            subject_id="NM_004",
            subject_start=100,
            subject_end=119,
            mismatches_3prime=1,  # 3' mismatch
        )
        right_hit = BlastHit(
            subject_id="NM_004", subject_start=350, subject_end=331, mismatches_3prime=0
        )
        result = checker._find_off_target_amplicons([left_hit], [right_hit])
        assert result == []

    def test_different_subjects_no_amplicon(self, checker: BlastSpecificity) -> None:
        """Hits on different subjects cannot form an amplicon."""
        left_hit = BlastHit(
            subject_id="NM_001", subject_start=100, subject_end=119, mismatches_3prime=0
        )
        right_hit = BlastHit(
            subject_id="NM_002", subject_start=350, subject_end=331, mismatches_3prime=0
        )
        result = checker._find_off_target_amplicons([left_hit], [right_hit])
        assert result == []


# ---------------------------------------------------------------------------
# Unit tests — _calculate_specificity_score
# ---------------------------------------------------------------------------


class TestCalculateSpecificityScore:
    """Test _calculate_specificity_score."""

    def test_perfect_score_no_off_targets(self, checker: BlastSpecificity) -> None:
        score = checker._calculate_specificity_score([], 1, 1)
        assert score == 100.0

    def test_score_reduced_by_off_target(self, checker: BlastSpecificity) -> None:
        lh = BlastHit(subject_id="X", subject_start=100, subject_end=119)
        rh = BlastHit(subject_id="X", subject_start=350, subject_end=331)
        off = OffTargetAmplicon(subject_id="X", amplicon_size=250, left_hit=lh, right_hit=rh)
        score = checker._calculate_specificity_score([off], 2, 2)
        assert score < 100.0

    def test_score_zero_for_many_off_targets(self, checker: BlastSpecificity) -> None:
        lh = BlastHit(subject_id="X", subject_start=100, subject_end=119)
        rh = BlastHit(subject_id="X", subject_start=350, subject_end=331)
        off = OffTargetAmplicon(subject_id="X", amplicon_size=250, left_hit=lh, right_hit=rh)
        score = checker._calculate_specificity_score([off] * 10, 5, 5)
        assert score == 0.0

    def test_many_hits_reduce_score(self, checker: BlastSpecificity) -> None:
        score = checker._calculate_specificity_score([], 50, 50)
        assert score < 100.0
        assert score >= 80.0


# ---------------------------------------------------------------------------
# Integration tests (mocked network) — check_pair & check_all
# ---------------------------------------------------------------------------


class TestCheckPair:
    """Test check_pair with mocked BLAST calls."""

    def test_check_pair_no_off_targets(self, checker: BlastSpecificity) -> None:
        pair = _make_pair()
        with patch.object(checker, "_blast_sequence", return_value=[]) as mock_blast:
            result = checker.check_pair(pair)
            assert isinstance(result, SpecificityResult)
            assert result.primer_pair_id == "pair_0"
            assert result.is_specific is True
            assert result.specificity_score == 100.0
            assert mock_blast.call_count == 2  # left + right primer

    def test_check_pair_with_off_target(self, checker: BlastSpecificity) -> None:
        pair = _make_pair()
        left_hit = BlastHit(
            subject_id="NM_999", subject_start=100, subject_end=119, mismatches_3prime=0
        )
        right_hit = BlastHit(
            subject_id="NM_999", subject_start=350, subject_end=331, mismatches_3prime=0
        )
        with patch.object(
            checker,
            "_blast_sequence",
            side_effect=[[left_hit], [right_hit]],
        ):
            result = checker.check_pair(pair)
            assert result.is_specific is False
            assert len(result.off_target_amplicons) == 1

    def test_check_pair_blasts_probe(self, checker: BlastSpecificity) -> None:
        pair = _make_pair()
        pair.probe = Oligo(sequence="ATCGATCGATCGATCGATCG", tm=68.0)
        with patch.object(checker, "_blast_sequence", return_value=[]) as mock_blast:
            result = checker.check_pair(pair)
            assert mock_blast.call_count == 3  # left + right + probe
            assert result.blast_hits_probe == []


class TestCheckAll:
    """Test check_all with mocked network."""

    def test_check_all_sets_specificity_result(self, checker: BlastSpecificity) -> None:
        pairs = [_make_pair("pair_0"), _make_pair("pair_1")]
        with patch.object(checker, "_blast_sequence", return_value=[]):
            result = checker.check_all(pairs)
            for pair in result:
                assert pair.specificity_result is not None
                assert isinstance(pair.specificity_result, SpecificityResult)

    def test_check_all_handles_exception(self, checker: BlastSpecificity) -> None:
        pairs = [_make_pair()]
        with patch.object(checker, "check_pair", side_effect=RuntimeError("network error")):
            result = checker.check_all(pairs)
            # Should not raise; should set specificity_result to failed result
            assert result[0].specificity_result is not None
            assert result[0].specificity_result.is_specific is False

    def test_check_all_empty(self, checker: BlastSpecificity) -> None:
        result = checker.check_all([])
        assert result == []
