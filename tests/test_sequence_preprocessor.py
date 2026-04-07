"""Tests for SequencePreprocessor."""

from __future__ import annotations

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import DesignTarget, ProcessedSequence
from tta_primer_design.modules.sequence_preprocessor import (
    SequencePreprocessor,
    _calculate_complexity_score,
    _extract_exon_junctions,
    _extract_sequence_str,
)

SEQ_SIMPLE = "ATCGATCGATCGATCGATCGATCGATCGATCG"
SEQ_GC50 = "GCGCGCATATATATATAT"  # roughly 50% GC
SEQ_LOWCOMP = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # all-A, low complexity
SEQ_MIXED = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def preprocessor(config: AppConfig) -> SequencePreprocessor:
    return SequencePreprocessor(config)


@pytest.fixture
def target() -> DesignTarget:
    return DesignTarget(target_id="T1")


# ---------------------------------------------------------------------------
# validate_sequence
# ---------------------------------------------------------------------------


class TestValidateSequence:
    def test_valid_atcg(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.validate_sequence("ATCGATCG") is True

    def test_valid_with_n(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.validate_sequence("ATCGNNNATCG") is True

    def test_lowercase_accepted(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.validate_sequence("atcgatcg") is True

    def test_empty_raises(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError, match="empty"):
            preprocessor.validate_sequence("")

    def test_invalid_char_raises(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError, match="invalid characters"):
            preprocessor.validate_sequence("ATCGX")

    def test_invalid_space_raises(self, preprocessor: SequencePreprocessor) -> None:
        with pytest.raises(ValueError, match="invalid characters"):
            preprocessor.validate_sequence("ATCG ATCG")


# ---------------------------------------------------------------------------
# calculate_gc_content
# ---------------------------------------------------------------------------


class TestCalculateGcContent:
    def test_all_gc(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("GCGCGCGC") == pytest.approx(1.0)

    def test_all_at(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("ATATATATAT") == pytest.approx(0.0)

    def test_fifty_percent(self, preprocessor: SequencePreprocessor) -> None:
        gc = preprocessor.calculate_gc_content("ATGCATGC")
        assert gc == pytest.approx(0.5)

    def test_empty_returns_zero(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("") == pytest.approx(0.0)

    def test_case_insensitive(self, preprocessor: SequencePreprocessor) -> None:
        assert preprocessor.calculate_gc_content("gcgc") == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# mask_low_complexity
# ---------------------------------------------------------------------------


class TestMaskLowComplexity:
    def test_all_same_base_fully_masked(self, preprocessor: SequencePreprocessor) -> None:
        masked = preprocessor.mask_low_complexity("A" * 20)
        assert all(c == "N" for c in masked)

    def test_complex_seq_unchanged(self, preprocessor: SequencePreprocessor) -> None:
        masked = preprocessor.mask_low_complexity(SEQ_MIXED)
        assert "N" not in masked

    def test_same_length(self, preprocessor: SequencePreprocessor) -> None:
        seq = SEQ_LOWCOMP + SEQ_MIXED
        masked = preprocessor.mask_low_complexity(seq)
        assert len(masked) == len(seq)

    def test_custom_threshold(self, preprocessor: SequencePreprocessor) -> None:
        # With very low threshold, all windows flagged
        masked = preprocessor.mask_low_complexity("ATCGATCG", threshold=0.1)
        assert "N" in masked

    def test_shorter_than_window(self, preprocessor: SequencePreprocessor) -> None:
        # sequence shorter than window — should not raise
        masked = preprocessor.mask_low_complexity("AAA", window=12)
        assert len(masked) == 3

    def test_lowercase_input(self, preprocessor: SequencePreprocessor) -> None:
        masked = preprocessor.mask_low_complexity("a" * 20)
        assert all(c == "N" for c in masked)


# ---------------------------------------------------------------------------
# process()
# ---------------------------------------------------------------------------


class TestProcess:
    def test_returns_processed_sequence(
        self, preprocessor: SequencePreprocessor, target: DesignTarget
    ) -> None:
        ps = preprocessor.process(SEQ_MIXED, target)
        assert isinstance(ps, ProcessedSequence)

    def test_sequence_same_length(
        self, preprocessor: SequencePreprocessor, target: DesignTarget
    ) -> None:
        ps = preprocessor.process(SEQ_MIXED, target)
        assert len(ps.sequence) == len(SEQ_MIXED)

    def test_gc_content_computed(
        self, preprocessor: SequencePreprocessor, target: DesignTarget
    ) -> None:
        ps = preprocessor.process(SEQ_MIXED, target)
        assert 0.0 <= ps.gc_content <= 1.0

    def test_complexity_score_range(
        self, preprocessor: SequencePreprocessor, target: DesignTarget
    ) -> None:
        ps = preprocessor.process(SEQ_MIXED, target)
        assert 0.0 <= ps.complexity_score <= 1.0

    def test_no_regions_by_default(
        self, preprocessor: SequencePreprocessor, target: DesignTarget
    ) -> None:
        ps = preprocessor.process(SEQ_MIXED, target)
        assert ps.excluded_regions == []
        assert ps.target_regions == []
        assert ps.included_region is None

    def test_excluded_regions_from_target(self, preprocessor: SequencePreprocessor) -> None:
        t = DesignTarget(target_id="T2", region_exclude=[(5, 15)])
        ps = preprocessor.process(SEQ_MIXED, t)
        assert (5, 15) in ps.excluded_regions

    def test_included_region_from_target(self, preprocessor: SequencePreprocessor) -> None:
        t = DesignTarget(target_id="T3", region_include=(10, 40))
        ps = preprocessor.process(SEQ_MIXED, t)
        assert ps.included_region is not None
        # Stored as (start, length)
        assert ps.included_region[0] == 10

    def test_invalid_sequence_raises(
        self, preprocessor: SequencePreprocessor, target: DesignTarget
    ) -> None:
        with pytest.raises(ValueError):
            preprocessor.process("ATCGXYZ", target)

    def test_biopython_seq_record_string(
        self, preprocessor: SequencePreprocessor, target: DesignTarget
    ) -> None:
        """Test với object giả lập SeqRecord (có .seq attribute)."""

        class FakeSeqRecord:
            def __init__(self, s: str) -> None:
                self.seq = s

        ps = preprocessor.process(FakeSeqRecord(SEQ_MIXED), target)
        assert isinstance(ps, ProcessedSequence)

    def test_exon_junctions_with_exon_coords(self, preprocessor: SequencePreprocessor) -> None:
        class FakeRecord:
            def __init__(self, seq: str) -> None:
                self.seq = seq
                self.exon_coords = [(0, 20), (20, 40), (40, 50)]

        t = DesignTarget(target_id="T4", exon_junction=True)
        ps = preprocessor.process(FakeRecord(SEQ_MIXED), t)
        assert ps.exon_junctions == [20, 40]


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


class TestExtractSequenceStr:
    def test_plain_str(self) -> None:
        assert _extract_sequence_str("atcg") == "ATCG"

    def test_object_with_seq_attr(self) -> None:
        class FakeRecord:
            seq = "gcta"

        assert _extract_sequence_str(FakeRecord()) == "GCTA"

    def test_unsupported_type_raises(self) -> None:
        with pytest.raises(TypeError):
            _extract_sequence_str(12345)  # type: ignore[arg-type]


class TestCalculateComplexityScore:
    def test_uniform_sequence_low(self) -> None:
        score = _calculate_complexity_score("A" * 40)
        assert score < 0.3

    def test_diverse_sequence_high(self) -> None:
        score = _calculate_complexity_score("ATCG" * 10)
        assert score > 0.9

    def test_empty_returns_zero(self) -> None:
        assert _calculate_complexity_score("") == pytest.approx(0.0)


class TestExtractExonJunctions:
    def test_no_exon_coords_no_features(self) -> None:
        junctions = _extract_exon_junctions("plain_string")
        assert junctions == []

    def test_with_exon_coords_attr(self) -> None:
        class FakeRecord:
            exon_coords = [(0, 100), (100, 200), (200, 300)]

        junctions = _extract_exon_junctions(FakeRecord())
        assert junctions == [100, 200]

    def test_single_exon_no_junction(self) -> None:
        class FakeRecord:
            exon_coords = [(0, 300)]

        junctions = _extract_exon_junctions(FakeRecord())
        assert junctions == []
