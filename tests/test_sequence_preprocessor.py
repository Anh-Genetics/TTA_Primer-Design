"""Tests cho SequencePreprocessor module."""

from __future__ import annotations

import pytest

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget, ProcessedSequence
from tta_primer_design.modules.sequence_preprocessor import SequencePreprocessor


@pytest.fixture()
def preprocessor() -> SequencePreprocessor:
    return SequencePreprocessor(AppConfig())


@pytest.fixture()
def simple_target() -> DesignTarget:
    return DesignTarget(target_id="ACTB", input_type="sequence")


class TestValidateSequence:
    """Tests for validate_sequence()."""

    def test_valid_sequence_returns_true(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.validate_sequence("ATCGN") is True

    def test_valid_all_bases(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.validate_sequence("AATTCCGGNN") is True

    def test_lowercase_valid(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.validate_sequence("atcgn") is True

    def test_raises_for_empty_string(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError, match="empty"):
            preprocessor.validate_sequence("")

    def test_raises_for_invalid_character_x(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError, match="invalid"):
            preprocessor.validate_sequence("ATCGX")

    def test_raises_for_invalid_character_u(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError, match="invalid"):
            preprocessor.validate_sequence("ATCGU")

    def test_raises_for_spaces(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError):
            preprocessor.validate_sequence("ATC G")

    def test_raises_for_numbers(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError):
            preprocessor.validate_sequence("ATC1G")


class TestCalculateGcContent:
    """Tests for calculate_gc_content()."""

    def test_equal_gc_at(self, preprocessor: SequencePreprocessor) -> None:
        # "ATCG" → 2 GC / 4 = 0.5
        assert preprocessor.calculate_gc_content("ATCG") == pytest.approx(0.5)

    def test_all_gc(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("GCGC") == pytest.approx(1.0)

    def test_all_at(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("ATAT") == pytest.approx(0.0)

    def test_empty_string_returns_zero(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("") == pytest.approx(0.0)

    def test_lowercase_handled(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("atcg") == pytest.approx(0.5)

    def test_with_n_bases(self, preprocessor: SequencePreprocessor) -> None:
        # "ATCGN" → G + C = 2, total = 5 → 0.4
        result = preprocessor.calculate_gc_content("ATCGN")
        assert result == pytest.approx(0.4)

    def test_returns_fraction_not_percent(self, preprocessor: SequencePreprocessor) -> None:
        result = preprocessor.calculate_gc_content("ATCGATCG")
        assert 0.0 <= result <= 1.0


class TestMaskLowComplexity:
    """Tests for mask_low_complexity()."""

    def test_masks_poly_a_run(self, preprocessor: SequencePreprocessor) -> None:
        # 16 A's exceeds default threshold in any window
        seq = "AAAAAAAAAAAAAAAA"  # 16 A's
        masked = preprocessor.mask_low_complexity(seq, window=12, threshold=0.7)
        assert "N" in masked

    def test_does_not_mask_complex_sequence(self, preprocessor: SequencePreprocessor) -> None:
        seq = "ATCGATCGATCGATCG"  # balanced sequence
        masked = preprocessor.mask_low_complexity(seq, window=12, threshold=0.9)
        # Should not mask all positions
        non_n = masked.replace("N", "")
        assert len(non_n) > 0

    def test_returns_same_length(self, preprocessor: SequencePreprocessor) -> None:
        seq = "ATCGATCGAAAAAAAAATCG"
        masked = preprocessor.mask_low_complexity(seq)
        assert len(masked) == len(seq)

    def test_mixed_sequence_partial_mask(self, preprocessor: SequencePreprocessor) -> None:
        # Poly-A region in the middle
        seq = "ATCGATCG" + "A" * 16 + "ATCGATCG"
        masked = preprocessor.mask_low_complexity(seq, window=12, threshold=0.7)
        assert "N" in masked


class TestProcess:
    """Tests for process()."""

    def test_returns_processed_sequence(
        self, preprocessor: SequencePreprocessor, simple_target: DesignTarget
    ) -> None:
        result = preprocessor.process("ATCGATCGATCG", simple_target)
        assert isinstance(result, ProcessedSequence)

    def test_converts_to_uppercase(
        self, preprocessor: SequencePreprocessor, simple_target: DesignTarget
    ) -> None:
        result = preprocessor.process("atcgatcgatcg", simple_target)
        assert result.sequence == result.sequence.upper()

    def test_gc_content_correct(
        self, preprocessor: SequencePreprocessor, simple_target: DesignTarget
    ) -> None:
        result = preprocessor.process("ATCG", simple_target)
        assert result.gc_content == pytest.approx(0.5)

    def test_applies_excluded_regions_from_target(self, preprocessor: SequencePreprocessor) -> None:
        target = DesignTarget(
            target_id="TEST",
            input_type="sequence",
            region_exclude=[(10, 20), (30, 40)],
        )
        result = preprocessor.process("A" * 100, target)
        assert result.excluded_regions == [(10, 20), (30, 40)]

    def test_included_region_passed_through(self, preprocessor: SequencePreprocessor) -> None:
        target = DesignTarget(
            target_id="TEST",
            input_type="sequence",
            region_include=(5, 80),
        )
        result = preprocessor.process("A" * 100, target)
        assert result.included_region == (5, 80)

    def test_no_excluded_regions_when_none(
        self, preprocessor: SequencePreprocessor, simple_target: DesignTarget
    ) -> None:
        result = preprocessor.process("ATCGATCGATCG", simple_target)
        assert result.excluded_regions == []

    def test_complexity_score_set(
        self, preprocessor: SequencePreprocessor, simple_target: DesignTarget
    ) -> None:
        result = preprocessor.process("ATCGATCGATCG", simple_target)
        assert result.complexity_score >= 0.0

    def test_raises_for_invalid_sequence(
        self, preprocessor: SequencePreprocessor, simple_target: DesignTarget
    ) -> None:
        with pytest.raises(ValueError):
            preprocessor.process("ATCGXYZ", simple_target)

    def test_accepts_biopython_seqrecord_like_object(
        self, preprocessor: SequencePreprocessor, simple_target: DesignTarget
    ) -> None:
        class FakeSeqRecord:
            def __init__(self, seq_str: str) -> None:
                self.seq = seq_str

            def __str__(self) -> str:
                return str(self.seq)

        # Without biopython, falls back to str()
        result = preprocessor.process(FakeSeqRecord("ATCGATCG"), simple_target)
        assert isinstance(result, ProcessedSequence)
