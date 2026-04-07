"""Tests cho data models — DesignTarget, Oligo, PrimerPair, DesignResult."""

from __future__ import annotations

from tta_primer_design.models import (
    DesignResult,
    DesignTarget,
    Oligo,
    PrimerPair,
    ProcessedSequence,
)


class TestDesignTarget:
    """Test DesignTarget dataclass."""

    def test_minimal_creation(self) -> None:
        t = DesignTarget(target_id="ACTB")
        assert t.target_id == "ACTB"
        assert t.input_type == "accession"
        assert t.organism == "Homo sapiens"
        assert t.design_mode == "qpcr"
        assert t.exon_junction is False

    def test_with_accession(self) -> None:
        t = DesignTarget(target_id="ACTB", accession="NM_001101")
        assert t.accession == "NM_001101"

    def test_with_sequence(self) -> None:
        t = DesignTarget(target_id="SEQ1", input_type="sequence", sequence="ATCG")
        assert t.sequence == "ATCG"

    def test_optional_fields_default_none(self) -> None:
        t = DesignTarget(target_id="X")
        assert t.accession is None
        assert t.sequence is None
        assert t.gene_name is None
        assert t.region_include is None
        assert t.region_exclude is None
        assert t.custom_params is None


class TestOligo:
    """Test Oligo dataclass."""

    def test_minimal_creation(self) -> None:
        o = Oligo(sequence="ATCGATCGATCG")
        assert o.sequence == "ATCGATCGATCG"
        assert o.length == 12  # auto-computed

    def test_explicit_length(self) -> None:
        o = Oligo(sequence="ATCG", length=4)
        assert o.length == 4

    def test_default_numerics(self) -> None:
        o = Oligo(sequence="ATCG")
        assert o.tm == 0.0
        assert o.gc_percent == 0.0
        assert o.penalty == 0.0


class TestPrimerPair:
    """Test PrimerPair dataclass."""

    def _make_pair(self) -> PrimerPair:
        left = Oligo(sequence="ATCGATCGATCGATCGATCG", tm=60.0)
        right = Oligo(sequence="GCTAGCTAGCTAGCTAGCTA", tm=60.0)
        return PrimerPair(pair_id="pair_0", left_primer=left, right_primer=right)

    def test_creation(self) -> None:
        pair = self._make_pair()
        assert pair.pair_id == "pair_0"
        assert pair.probe is None
        assert pair.amplicon_size == 0
        assert pair.score == 0.0

    def test_snp_flags_default_empty(self) -> None:
        pair = self._make_pair()
        assert pair.snp_flags == []


class TestProcessedSequence:
    """Test ProcessedSequence dataclass."""

    def test_minimal_creation(self) -> None:
        ps = ProcessedSequence(sequence="ATCGATCG")
        assert ps.sequence == "ATCGATCG"
        assert ps.excluded_regions == []
        assert ps.target_regions == []
        assert ps.exon_junctions == []
        assert ps.gc_content == 0.0

    def test_with_regions(self) -> None:
        ps = ProcessedSequence(
            sequence="ATCGATCG",
            excluded_regions=[(10, 20), (50, 60)],
            target_regions=[(30, 40)],
        )
        assert len(ps.excluded_regions) == 2
        assert len(ps.target_regions) == 1


class TestDesignResult:
    """Test DesignResult dataclass."""

    def test_success_result(self) -> None:
        target = DesignTarget(target_id="ACTB")
        result = DesignResult(target=target)
        assert result.status == "success"
        assert result.primer_pairs == []
        assert result.error is None

    def test_failed_result(self) -> None:
        target = DesignTarget(target_id="ACTB")
        result = DesignResult(target=target, status="failed", error="Network error")
        assert result.status == "failed"
        assert result.error == "Network error"
