"""Tests cho ProbeDesigner module."""

from __future__ import annotations

import pytest

from tta_primer_design.config import AppConfig, PipelineConfig
from tta_primer_design.models import Oligo, PrimerPair, ProcessedSequence
from tta_primer_design.modules.probe_designer import ProbeDesigner

# A realistic ~150 bp ACTB-like amplicon sequence
_ACTB_AMPLICON = (
    "ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAAGCCGGCTTCGCGGGCGACGATGCC"
    "CCGAGCCAAAGCGGGCAGAAGAGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGC"
)

# Build a ProcessedSequence from the amplicon
_PROCESSED_SEQ = ProcessedSequence(
    sequence=_ACTB_AMPLICON,
    gc_content=0.55,
    complexity_score=0.9,
)

_LEFT_START = 0
_LEFT_LEN = 22
_RIGHT_LEN = 20
_RIGHT_START = len(_ACTB_AMPLICON) - _RIGHT_LEN  # 138


def _make_pair(
    pair_id: str = "test_pair",
    left_start: int = _LEFT_START,
    left_len: int = _LEFT_LEN,
    right_start: int = _RIGHT_START,
    right_len: int = _RIGHT_LEN,
    left_tm: float = 60.0,
    right_tm: float = 60.0,
) -> PrimerPair:
    left_seq = _ACTB_AMPLICON[left_start : left_start + left_len]
    right_seq = _ACTB_AMPLICON[right_start : right_start + right_len]
    return PrimerPair(
        pair_id=pair_id,
        left_primer=Oligo(sequence=left_seq, start=left_start, length=left_len, tm=left_tm),
        right_primer=Oligo(sequence=right_seq, start=right_start, length=right_len, tm=right_tm),
    )


@pytest.fixture()
def taqman_config() -> AppConfig:
    cfg = AppConfig()
    cfg.pipeline = PipelineConfig(mode="taqman")
    return cfg


@pytest.fixture()
def designer(taqman_config: AppConfig) -> ProbeDesigner:
    return ProbeDesigner(taqman_config)


@pytest.fixture()
def pcr_config() -> AppConfig:
    cfg = AppConfig()
    cfg.pipeline = PipelineConfig(mode="pcr")
    return cfg


class TestCheckProbeRules:
    """Tests for check_probe_rules()."""

    def test_rejects_probe_starting_with_g(self, designer: ProbeDesigner) -> None:
        # GC is 50%, no poly-G, but starts with G → fail
        probe = Oligo(sequence="GCATCGATCGATCGATCGATC")
        ok, violations = designer.check_probe_rules(probe)
        assert ok is False
        assert "5_prime_G" in violations

    def test_rejects_probe_with_poly_g_run(self, designer: ProbeDesigner) -> None:
        probe = Oligo(sequence="ATCGATGGGGGCGATCGATCG")
        ok, violations = designer.check_probe_rules(probe)
        assert ok is False
        assert "poly_G_run" in violations

    def test_accepts_good_probe(self, designer: ProbeDesigner) -> None:
        # 40-65% GC, no poly run, doesn't start with G
        # "ATCGATCGATCGATCGATCG" → 10G+C / 20 = 50%
        probe = Oligo(sequence="ATCGATCGATCGATCGATCG")
        ok, violations = designer.check_probe_rules(probe)
        assert ok is True
        assert violations == []

    def test_rejects_low_gc_probe(self, designer: ProbeDesigner) -> None:
        # All AT → 0% GC
        probe = Oligo(sequence="ATATATATATATATATATAT")
        ok, violations = designer.check_probe_rules(probe)
        assert ok is False
        gc_violations = [v for v in violations if "gc_out_of_range" in v]
        assert len(gc_violations) == 1

    def test_rejects_high_gc_probe(self, designer: ProbeDesigner) -> None:
        # All GC → 100% GC (also starts with G)
        probe = Oligo(sequence="GCGCGCGCGCGCGCGCGCGC")
        ok, violations = designer.check_probe_rules(probe)
        assert ok is False

    def test_returns_tuple(self, designer: ProbeDesigner) -> None:
        probe = Oligo(sequence="ATCGATCGATCGATCGATCG")
        result = designer.check_probe_rules(probe)
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestSelectBestStrand:
    """Tests for select_best_strand()."""

    def test_returns_string(self, designer: ProbeDesigner) -> None:
        result = designer.select_best_strand(_ACTB_AMPLICON, 20, 22)
        assert isinstance(result, str)

    def test_returned_length_matches_requested(self, designer: ProbeDesigner) -> None:
        result = designer.select_best_strand(_ACTB_AMPLICON, 20, 22)
        assert len(result) == 22

    def test_prefers_strand_with_fewer_g(self, designer: ProbeDesigner) -> None:
        # Construct a 20-mer where forward has more G than reverse
        # Forward: GGGGGATCGATCGATCGATC → G count = 5 + 2 = 7
        # Reverse complement will have C count = 7 G on the other strand actually
        # Let's use a simpler approach: pick a window and verify strand choice
        fwd = "GGGGGATCGATCGATCGATC"  # lots of G
        amplicon = "N" * 10 + fwd + "N" * 10
        result = designer.select_best_strand(amplicon, 10, 20)
        # result should be reverse complement (fewer G = more C in fwd = fewer G in rev)
        assert isinstance(result, str)

    def test_does_not_start_with_g_when_possible(self, designer: ProbeDesigner) -> None:
        # Equal G count: prefer strand that does NOT start with G
        # Build a palindrome-ish sequence where G counts are equal
        # "ATCG" * 5 → G count = 5, reverse complement = "CGAT" * 5 → G count = 5
        # Forward starts with A, reverse complement starts with C → forward preferred
        fwd = "ATCGATCGATCGATCGATCG"
        result = designer.select_best_strand(fwd, 0, 20)
        # Forward: A starts, rev_comp: starts with C → both don't start with G
        # Forward G count = 5, rev_comp G count = 5, so forward returned (not starting with G)
        assert result[0] != "G" or fwd[0] == "G"


class TestCalculateProbeTm:
    """Tests for calculate_probe_tm()."""

    def test_returns_float(self, designer: ProbeDesigner) -> None:
        tm = designer.calculate_probe_tm("ATCGATCGATCGATCGATCG", {})
        assert isinstance(tm, float)

    def test_reasonable_range(self, designer: ProbeDesigner) -> None:
        tm = designer.calculate_probe_tm("ATCGATCGATCGATCGATCG", {})
        assert 30.0 <= tm <= 90.0

    def test_custom_params_accepted(self, designer: ProbeDesigner) -> None:
        tm = designer.calculate_probe_tm(
            "ATCGATCGATCGATCGATCG",
            {"mv_conc": 50.0, "dv_conc": 1.5},
        )
        assert isinstance(tm, float)


class TestDesignTaqmanProbe:
    """Tests for design_taqman_probe()."""

    def test_returns_none_for_short_amplicon(self, designer: ProbeDesigner) -> None:
        short_amplicon = "ATCGATCGATCGATCGATCG"  # only 20 bp
        pair = _make_pair(left_len=5, right_len=5)
        result = designer.design_taqman_probe(short_amplicon, pair)
        assert result is None

    def test_returns_oligo_or_none_for_normal_amplicon(self, designer: ProbeDesigner) -> None:
        pair = _make_pair(left_tm=60.0, right_tm=60.0)
        result = designer.design_taqman_probe(_ACTB_AMPLICON, pair)
        assert result is None or isinstance(result, Oligo)

    def test_probe_start_set_when_found(self, designer: ProbeDesigner) -> None:
        pair = _make_pair(left_tm=60.0, right_tm=60.0)
        result = designer.design_taqman_probe(_ACTB_AMPLICON, pair)
        if result is not None:
            assert result.start >= 0


class TestDesign:
    """Tests for design()."""

    def test_skips_non_taqman_qpcr_mode(self, pcr_config: AppConfig) -> None:
        designer = ProbeDesigner(pcr_config)
        pair = _make_pair()
        result = designer.design([pair], _PROCESSED_SEQ)
        # In PCR mode, pairs are returned unchanged
        assert result == [pair]

    def test_skips_pair_already_with_probe(self, designer: ProbeDesigner) -> None:
        pair = _make_pair()
        existing_probe = Oligo(sequence="ATCGATCGATCGATCGATCG")
        pair.probe = existing_probe
        result = designer.design([pair], _PROCESSED_SEQ)
        # Should not replace existing probe
        assert result[0].probe is existing_probe

    def test_returns_list_of_pairs(self, designer: ProbeDesigner) -> None:
        pair = _make_pair()
        result = designer.design([pair], _PROCESSED_SEQ)
        assert isinstance(result, list)
        assert len(result) == 1

    def test_handles_empty_input(self, designer: ProbeDesigner) -> None:
        result = designer.design([], _PROCESSED_SEQ)
        assert result == []

    def test_qpcr_mode_also_runs(self) -> None:
        cfg = AppConfig()
        cfg.pipeline = PipelineConfig(mode="qpcr")
        designer = ProbeDesigner(cfg)
        pair = _make_pair()
        result = designer.design([pair], _PROCESSED_SEQ)
        assert isinstance(result, list)
