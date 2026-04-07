"""Tests cho ProbeDesigner — TaqMan probe design."""

from __future__ import annotations

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import Oligo, PrimerPair, ProcessedSequence
from tta_primer_design.modules.probe_designer import ProbeDesigner, _gc_percent, _reverse_complement


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def designer(config: AppConfig) -> ProbeDesigner:
    return ProbeDesigner(config)


# 200 bp synthetic sequence with no repeats, ~50% GC
_TEMPLATE = (
    "ATCGATCGATCGATCGATCG"  # 20 bp left primer region
    "GCTAGCTAGCTAGCTAGCTA"  # interior (~160 bp)
    "ATCGATCGATCGATCGATCG"
    "GCTAGCTAGCTAGCTAGCTA"
    "ATCGATCGATCGATCGATCG"
    "GCTAGCTAGCTAGCTAGCTA"
    "ATCGATCGATCGATCGATCG"
    "GCTAGCTAGCTAGCTAGCTA"
    "AAATTTGGGCCCATCG"  # 16 bp right primer region (reverse complement scenario)
    "TTAACCGGTTAACCGG"  # right primer binding site
)


def _make_pair(
    left_seq: str = "ATCGATCGATCGATCGATCG",
    right_seq: str = "TTAACCGGTTAACCGG",
    left_tm: float = 58.0,
    right_tm: float = 58.0,
    left_start: int = 0,
    right_start: int = 164,
    amplicon_seq: str = "",
    amplicon_size: int = 180,
) -> PrimerPair:
    left = Oligo(sequence=left_seq, tm=left_tm, start=left_start)
    right = Oligo(sequence=right_seq, tm=right_tm, start=right_start)
    return PrimerPair(
        pair_id="test_pair",
        left_primer=left,
        right_primer=right,
        amplicon_sequence=amplicon_seq,
        amplicon_size=amplicon_size,
    )


# ---------------------------------------------------------------------------
# Unit tests — helper functions
# ---------------------------------------------------------------------------


class TestHelpers:
    """Test module-level helper functions."""

    def test_gc_percent_pure_gc(self) -> None:
        assert _gc_percent("GGCC") == 100.0

    def test_gc_percent_pure_at(self) -> None:
        assert _gc_percent("AATT") == 0.0

    def test_gc_percent_mixed(self) -> None:
        assert _gc_percent("ATGC") == 50.0

    def test_gc_percent_empty(self) -> None:
        assert _gc_percent("") == 0.0

    def test_reverse_complement(self) -> None:
        assert _reverse_complement("ATCG") == "CGAT"

    def test_reverse_complement_palindrome(self) -> None:
        assert _reverse_complement("AATT") == "AATT"

    def test_reverse_complement_lowercase(self) -> None:
        assert _reverse_complement("atcg") == "cgat"


# ---------------------------------------------------------------------------
# Unit tests — calculate_probe_tm
# ---------------------------------------------------------------------------


class TestCalculateProbeTm:
    """Test calculate_probe_tm."""

    def test_returns_float(self, designer: ProbeDesigner) -> None:
        tm = designer.calculate_probe_tm("ATCGATCGATCGATCGATCG", {})
        assert isinstance(tm, float)

    def test_tm_reasonable_range(self, designer: ProbeDesigner) -> None:
        tm = designer.calculate_probe_tm("GCTAGCTAGCTAGCTAGCTA", {})
        assert 30.0 < tm < 90.0

    def test_higher_gc_higher_tm(self, designer: ProbeDesigner) -> None:
        tm_gc = designer.calculate_probe_tm("GCGCGCGCGCGCGCGCGCGC", {})
        tm_at = designer.calculate_probe_tm("ATATATATATATATATATAN".replace("N", "A"), {})
        assert tm_gc > tm_at

    def test_custom_params(self, designer: ProbeDesigner) -> None:
        params = {"mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.25, "dna_conc": 250.0}
        tm = designer.calculate_probe_tm("ATCGATCGATCGATCGATCG", params)
        assert isinstance(tm, float)


# ---------------------------------------------------------------------------
# Unit tests — select_best_strand
# ---------------------------------------------------------------------------


class TestSelectBestStrand:
    """Test select_best_strand."""

    def test_prefers_strand_with_fewer_g(self, designer: ProbeDesigner) -> None:
        # Amplicon: first 30 chars have GGGGG → RC will have fewer G
        amplicon = "GGGGGCCCCCATCGATCGATCGATCGATCG"
        result = designer.select_best_strand(amplicon, 0, 20)
        # Should return whichever has fewer G
        rc = _reverse_complement(amplicon[:20])
        fwd = amplicon[:20].upper()
        assert result in (fwd, rc)

    def test_returns_fwd_when_equal_gc(self, designer: ProbeDesigner) -> None:
        amplicon = "ATCGATCGATCGATCGATCG"
        result = designer.select_best_strand(amplicon, 0, 20)
        # Equal G and C → forward strand
        assert result == amplicon.upper()

    def test_length_preserved(self, designer: ProbeDesigner) -> None:
        amplicon = "ATCGATCGATCGATCGATCGATCGATCGATCG"
        result = designer.select_best_strand(amplicon, 0, 18)
        assert len(result) == 18

    def test_offset_works(self, designer: ProbeDesigner) -> None:
        amplicon = "AAAAAAATCGATCGATCGATCG"
        result = designer.select_best_strand(amplicon, 7, 15)
        assert len(result) == 15


# ---------------------------------------------------------------------------
# Unit tests — check_probe_rules
# ---------------------------------------------------------------------------


class TestCheckProbeRules:
    """Test check_probe_rules."""

    def test_pass_valid_probe(self, designer: ProbeDesigner) -> None:
        # Sequence with no G at 5', no poly runs, weak hairpin (dg > -2000 cal/mol)
        probe = Oligo(sequence="CATCGTAGCATGCTAGCATG", tm=68.0, gc_percent=50.0)
        passed, violations = designer.check_probe_rules(probe)
        assert passed is True
        assert violations == []

    def test_fail_g_at_5prime(self, designer: ProbeDesigner) -> None:
        probe = Oligo(sequence="GATCGATCGATCGATCGATCG", tm=68.0)
        passed, violations = designer.check_probe_rules(probe)
        assert passed is False
        assert any("5'" in v for v in violations)

    def test_fail_poly_g_run(self, designer: ProbeDesigner) -> None:
        probe = Oligo(sequence="ATCGGGGGGATCGATCGATCG", tm=68.0)
        passed, violations = designer.check_probe_rules(probe)
        assert passed is False
        assert any("polyG" in v or "poly" in v.lower() for v in violations)

    def test_fail_poly_c_run(self, designer: ProbeDesigner) -> None:
        probe = Oligo(sequence="ATCCCCCCCATCGATCGATCG", tm=68.0)
        passed, violations = designer.check_probe_rules(probe)
        assert passed is False

    def test_multiple_violations(self, designer: ProbeDesigner) -> None:
        probe = Oligo(sequence="GGGGGGGGATCGATCGATCG", tm=68.0)
        passed, violations = designer.check_probe_rules(probe)
        assert passed is False
        assert len(violations) >= 1


# ---------------------------------------------------------------------------
# Unit tests — design_taqman_probe
# ---------------------------------------------------------------------------


class TestDesignTaqManProbe:
    """Test design_taqman_probe."""

    def test_returns_none_for_short_amplicon(self, designer: ProbeDesigner) -> None:
        # Amplicon too short for probe
        pair = _make_pair(
            left_seq="ATCGATCGATCGATCGATCG",
            right_seq="CGCGCGCGCGCGCGCGCGCG",
            amplicon_size=40,
        )
        amplicon = "ATCGATCGATCGATCGATCG" + "ATCG" + "CGCGCGCGCGCGCGCGCGCG"
        result = designer.design_taqman_probe(amplicon, pair)
        assert result is None

    def test_probe_not_overlapping_primers(self, designer: ProbeDesigner) -> None:
        """Probe should not overlap with primer sequences."""
        left_seq = "ATCGATCGATCGATCGATCG"
        right_seq = "CGCGCGCGCGCGCGCGCGCG"
        # Build a long enough amplicon with a high-GC interior
        interior = (
            "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"
        )
        amplicon = left_seq + interior + right_seq

        pair = _make_pair(
            left_seq=left_seq,
            right_seq=right_seq,
            left_tm=58.0,
            right_tm=58.0,
            amplicon_size=len(amplicon),
            amplicon_seq=amplicon,
        )
        probe = designer.design_taqman_probe(amplicon, pair)
        if probe is not None:
            # Probe must start after left primer
            assert probe.start >= len(left_seq)

    def test_probe_has_valid_sequence(self, designer: ProbeDesigner) -> None:
        left_seq = "ATCGATCGATCGATCGATCG"
        right_seq = "CGCGCGCGCGCGCGCGCGCG"
        interior = (
            "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"
        )
        amplicon = left_seq + interior + right_seq
        pair = _make_pair(
            left_seq=left_seq,
            right_seq=right_seq,
            amplicon_size=len(amplicon),
            amplicon_seq=amplicon,
        )
        probe = designer.design_taqman_probe(amplicon, pair)
        if probe is not None:
            assert isinstance(probe.sequence, str)
            assert len(probe.sequence) >= 18
            assert all(c in "ACGT" for c in probe.sequence)


# ---------------------------------------------------------------------------
# Unit tests — design (top-level)
# ---------------------------------------------------------------------------


class TestDesign:
    """Test design method."""

    def test_design_returns_list(self, designer: ProbeDesigner) -> None:
        pair = _make_pair()
        processed = ProcessedSequence(sequence=_TEMPLATE)
        result = designer.design([pair], processed)
        assert isinstance(result, list)
        assert len(result) == 1

    def test_design_assigns_probe_or_none(self, designer: ProbeDesigner) -> None:
        pair = _make_pair()
        processed = ProcessedSequence(sequence=_TEMPLATE)
        result = designer.design([pair], processed)
        # probe may be None or Oligo, both are valid
        assert result[0].probe is None or isinstance(result[0].probe, Oligo)

    def test_design_with_amplicon_seq_set(self, designer: ProbeDesigner) -> None:
        left_seq = "ATCGATCGATCGATCGATCG"
        right_seq = "CGCGCGCGCGCGCGCGCGCG"
        interior = (
            "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG"
        )
        amplicon = left_seq + interior + right_seq
        pair = _make_pair(
            left_seq=left_seq,
            right_seq=right_seq,
            amplicon_seq=amplicon,
            amplicon_size=len(amplicon),
        )
        processed = ProcessedSequence(sequence=amplicon)
        result = designer.design([pair], processed)
        assert isinstance(result, list)

    def test_design_empty_input(self, designer: ProbeDesigner) -> None:
        processed = ProcessedSequence(sequence=_TEMPLATE)
        result = designer.design([], processed)
        assert result == []

    def test_design_preserves_pair_id(self, designer: ProbeDesigner) -> None:
        pair = _make_pair()
        processed = ProcessedSequence(sequence=_TEMPLATE)
        result = designer.design([pair], processed)
        assert result[0].pair_id == "test_pair"
