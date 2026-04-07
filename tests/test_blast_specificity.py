"""Tests cho blast_specificity module — helper functions and BlastSpecificity class."""

from __future__ import annotations

from unittest.mock import patch

from tta_primer_design.config import AppConfig
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.blast_specificity import (
    BlastHit,
    BlastSpecificity,
    OffTargetAmplicon,
    SpecificityResult,
    _calc_3prime_mismatches,
    _find_off_target_amplicons,
)

# ---------------------------------------------------------------------------
# MockHSP helper
# ---------------------------------------------------------------------------


class MockHSP:
    def __init__(self, query: str, sbjct: str, query_end: int) -> None:
        self.query = query
        self.sbjct = sbjct
        self.query_end = query_end


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------


def _make_blast_hit(
    subject_id: str = "seq1",
    subject_start: int = 100,
    subject_end: int = 120,
    mismatches: int = 0,
) -> BlastHit:
    return BlastHit(
        subject_id=subject_id,
        subject_start=subject_start,
        subject_end=subject_end,
        mismatches=mismatches,
    )


def _make_pair(pair_id: str = "pair_01") -> PrimerPair:
    return PrimerPair(
        pair_id=pair_id,
        left_primer=Oligo(sequence="GCAAGGAATGGTTTCAGAAATCCA"),
        right_primer=Oligo(sequence="CAGGACTCCATGTCGTCCA"),
    )


def _make_spec_result(off_targets: int = 0) -> SpecificityResult:
    left_hits = [_make_blast_hit("seq1", 100, 120, 0)]
    right_hits = [_make_blast_hit("seq1", 300, 280, 0)]
    off = [
        OffTargetAmplicon(
            subject_id="seq1",
            amplicon_size=200,
            left_hit=left_hits[0],
            right_hit=right_hits[0],
        )
    ] * off_targets
    return SpecificityResult(
        primer_pair_id="pair_01",
        blast_hits_left=left_hits,
        blast_hits_right=right_hits,
        off_target_amplicons=off,
        specificity_score=max(0.0, 100.0 - off_targets * 20.0),
        is_specific=(off_targets == 0),
    )


# ---------------------------------------------------------------------------
# _calc_3prime_mismatches
# ---------------------------------------------------------------------------


class TestCalc3PrimeMismatches:
    """Tests for _calc_3prime_mismatches()."""

    def test_returns_zero_when_3prime_not_in_alignment(self) -> None:
        hsp = MockHSP(query="ATCGATCG", sbjct="ATCGATCG", query_end=8)
        # query_len=12 > query_end=8 → 3' end not in alignment
        result = _calc_3prime_mismatches(hsp, query_len=12)
        assert result == 0

    def test_returns_zero_for_perfect_match(self) -> None:
        seq = "ATCGATCG"
        hsp = MockHSP(query=seq, sbjct=seq, query_end=len(seq))
        result = _calc_3prime_mismatches(hsp, query_len=len(seq))
        assert result == 0

    def test_counts_one_mismatch_at_3prime(self) -> None:
        # Last base differs: G vs C
        query = "ATCGATCG"
        sbjct = "ATCGATCC"
        hsp = MockHSP(query=query, sbjct=sbjct, query_end=len(query))
        result = _calc_3prime_mismatches(hsp, query_len=len(query), prime_len=3)
        assert result == 1

    def test_counts_multiple_mismatches(self) -> None:
        # Last 3 bases all differ: ATG → CGT
        query = "ATCGATCGATG"
        sbjct = "ATCGATCGCGT"
        hsp = MockHSP(query=query, sbjct=sbjct, query_end=len(query))
        result = _calc_3prime_mismatches(hsp, query_len=len(query), prime_len=3)
        assert result == 3

    def test_only_counts_prime_len_bases(self) -> None:
        # First base differs (5' end), last 3 are perfect
        query = "ATCGATCGATCG"
        sbjct = "TTCGATCGATCG"
        hsp = MockHSP(query=query, sbjct=sbjct, query_end=len(query))
        result = _calc_3prime_mismatches(hsp, query_len=len(query), prime_len=3)
        assert result == 0

    def test_returns_int(self) -> None:
        hsp = MockHSP(query="ATCG", sbjct="ATCG", query_end=4)
        result = _calc_3prime_mismatches(hsp, query_len=4)
        assert isinstance(result, int)


# ---------------------------------------------------------------------------
# _find_off_target_amplicons
# ---------------------------------------------------------------------------


class TestFindOffTargetAmplicons:
    """Tests for _find_off_target_amplicons()."""

    def test_returns_empty_when_no_matching_subject_ids(self) -> None:
        left = [_make_blast_hit("seq1", 100, 120)]
        right = [_make_blast_hit("seq2", 300, 280)]  # different subject
        result = _find_off_target_amplicons(left, right)
        assert result == []

    def test_finds_amplicon_on_same_subject(self) -> None:
        left = [_make_blast_hit("seq1", subject_start=100, subject_end=120)]
        # Reverse strand: subject_start > subject_end
        right = [BlastHit(subject_id="seq1", subject_start=300, subject_end=280, mismatches=0)]
        result = _find_off_target_amplicons(left, right, max_size=2000)
        assert len(result) == 1
        assert result[0].subject_id == "seq1"
        assert result[0].amplicon_size == 201  # 300 - 100 + 1

    def test_skips_amplicon_exceeding_max_size(self) -> None:
        left = [_make_blast_hit("seq1", subject_start=100, subject_end=120)]
        right = [BlastHit(subject_id="seq1", subject_start=3000, subject_end=2980, mismatches=0)]
        result = _find_off_target_amplicons(left, right, max_size=500)
        assert result == []

    def test_skips_left_hit_with_too_many_mismatches(self) -> None:
        # left hit has 3 mismatches, max_mismatches=2 → skip
        left = [_make_blast_hit("seq1", subject_start=100, subject_end=120, mismatches=3)]
        right = [BlastHit(subject_id="seq1", subject_start=300, subject_end=280, mismatches=0)]
        result = _find_off_target_amplicons(left, right, max_mismatches=2)
        assert result == []

    def test_skips_right_hit_with_too_many_mismatches(self) -> None:
        left = [_make_blast_hit("seq1", subject_start=100, subject_end=120, mismatches=0)]
        right = [BlastHit(subject_id="seq1", subject_start=300, subject_end=280, mismatches=3)]
        result = _find_off_target_amplicons(left, right, max_mismatches=2)
        assert result == []

    def test_skips_left_hit_on_reverse_strand(self) -> None:
        # Left hit on reverse strand (subject_start > subject_end) should be ignored
        left = [BlastHit(subject_id="seq1", subject_start=120, subject_end=100, mismatches=0)]
        right = [BlastHit(subject_id="seq1", subject_start=300, subject_end=280, mismatches=0)]
        result = _find_off_target_amplicons(left, right)
        assert result == []

    def test_multiple_hits_same_subject(self) -> None:
        left = [
            _make_blast_hit("seq1", 100, 120),
            _make_blast_hit("seq1", 200, 220),
        ]
        right = [BlastHit(subject_id="seq1", subject_start=500, subject_end=480, mismatches=0)]
        result = _find_off_target_amplicons(left, right, max_size=5000)
        assert len(result) == 2

    def test_returns_off_target_amplicon_objects(self) -> None:
        left = [_make_blast_hit("seq1", 100, 120)]
        right = [BlastHit(subject_id="seq1", subject_start=300, subject_end=280, mismatches=0)]
        result = _find_off_target_amplicons(left, right)
        assert isinstance(result[0], OffTargetAmplicon)

    def test_empty_inputs(self) -> None:
        assert _find_off_target_amplicons([], []) == []


# ---------------------------------------------------------------------------
# BlastSpecificity.check_pair
# ---------------------------------------------------------------------------


class TestBlastSpecificityCheckPair:
    """Tests for BlastSpecificity.check_pair() with mocked BLAST."""

    def test_returns_specificity_result(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        pair = _make_pair()

        mock_hits: list[BlastHit] = [
            BlastHit(subject_id="NM_001101", mismatches=0, subject_start=100, subject_end=120)
        ]
        with patch(
            "tta_primer_design.modules.blast_specificity._blast_oligo_ncbi",
            return_value=mock_hits,
        ):
            result = checker.check_pair(pair)

        assert isinstance(result, SpecificityResult)
        assert result.primer_pair_id == "pair_01"

    def test_blast_called_for_left_and_right(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        pair = _make_pair()

        with patch(
            "tta_primer_design.modules.blast_specificity._blast_oligo_ncbi",
            return_value=[],
        ) as mock_blast:
            checker.check_pair(pair)

        assert mock_blast.call_count == 2

    def test_blast_called_three_times_with_probe(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        pair = _make_pair()
        pair.probe = Oligo(sequence="TGCAGCCACACTTTCTACAATGAGC")

        with patch(
            "tta_primer_design.modules.blast_specificity._blast_oligo_ncbi",
            return_value=[],
        ) as mock_blast:
            checker.check_pair(pair)

        assert mock_blast.call_count == 3

    def test_result_contains_blast_hits(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        pair = _make_pair()

        left_hits = [BlastHit(subject_id="NM_001", mismatches=0)]
        right_hits = [BlastHit(subject_id="NM_002", mismatches=1)]

        with patch(
            "tta_primer_design.modules.blast_specificity._blast_oligo_ncbi",
            side_effect=[left_hits, right_hits],
        ):
            result = checker.check_pair(pair)

        assert result.blast_hits_left == left_hits
        assert result.blast_hits_right == right_hits

    def test_no_off_targets_is_specific(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        pair = _make_pair()

        with patch(
            "tta_primer_design.modules.blast_specificity._blast_oligo_ncbi",
            return_value=[],
        ):
            result = checker.check_pair(pair)

        assert result.is_specific is True
        assert result.off_target_amplicons == []
        assert result.specificity_score == 100.0


# ---------------------------------------------------------------------------
# BlastSpecificity.check_all
# ---------------------------------------------------------------------------


class TestBlastSpecificityCheckAll:
    """Tests for BlastSpecificity.check_all() exception handling."""

    def test_returns_all_pairs(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        pairs = [_make_pair("p1"), _make_pair("p2")]

        with patch.object(checker, "check_pair", return_value=_make_spec_result()):
            result = checker.check_all(pairs)

        assert len(result) == 2

    def test_sets_specificity_result_on_pairs(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        pair = _make_pair()

        spec = _make_spec_result()
        with patch.object(checker, "check_pair", return_value=spec):
            checker.check_all([pair])

        assert pair.specificity_result is spec

    def test_catches_exception_per_pair_without_crashing(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        good = _make_pair("good")
        bad = _make_pair("bad")

        def side_effect(pair: PrimerPair, **kwargs: object) -> SpecificityResult:
            if pair.pair_id == "bad":
                raise RuntimeError("BLAST API error")
            return _make_spec_result()

        with patch.object(checker, "check_pair", side_effect=side_effect):
            result = checker.check_all([good, bad])

        # Should not raise — both pairs returned, bad pair has no spec result
        assert len(result) == 2
        assert good.specificity_result is not None
        assert bad.specificity_result is None

    def test_empty_list(self) -> None:
        cfg = AppConfig()
        checker = BlastSpecificity(cfg)
        result = checker.check_all([])
        assert result == []
