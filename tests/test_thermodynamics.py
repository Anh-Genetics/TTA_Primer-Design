"""Tests for Thermodynamics module."""

from __future__ import annotations

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.thermodynamics import (
    Thermodynamics,
    ThermoProfile,
    _calc_end_stability,
)

# Typical qPCR primer sequences
PRIMER_LEFT = "ATCGATCGATCGATCGATCG"
PRIMER_RIGHT = "GCTAGCTAGCTAGCTAGCTA"
PROBE_SEQ = "GCATGCATGCATGCATGCATGCATGC"

# A primer with known issues
REPEAT_SEQ = "AAAAAAAAAAATCGATCG"  # 12 A's — repeat run
LOW_GC_SEQ = "ATATATATATATATATAT"  # all AT, no GC clamp


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def thermo(config: AppConfig) -> Thermodynamics:
    return Thermodynamics(config)


@pytest.fixture
def left_oligo() -> Oligo:
    return Oligo(sequence=PRIMER_LEFT)


@pytest.fixture
def right_oligo() -> Oligo:
    return Oligo(sequence=PRIMER_RIGHT)


@pytest.fixture
def primer_pair(left_oligo: Oligo, right_oligo: Oligo) -> PrimerPair:
    return PrimerPair(pair_id="pair_0", left_primer=left_oligo, right_primer=right_oligo)


# ---------------------------------------------------------------------------
# calculate_tm
# ---------------------------------------------------------------------------


class TestCalculateTm:
    def test_returns_float(self, thermo: Thermodynamics) -> None:
        tm = thermo.calculate_tm(PRIMER_LEFT)
        assert isinstance(tm, float)

    def test_reasonable_range(self, thermo: Thermodynamics) -> None:
        tm = thermo.calculate_tm(PRIMER_LEFT)
        assert 40.0 <= tm <= 80.0

    def test_longer_seq_higher_tm(self, thermo: Thermodynamics) -> None:
        short = "ATCGATCG"
        long = "ATCGATCGATCGATCGATCG"
        assert thermo.calculate_tm(long) > thermo.calculate_tm(short)

    def test_gc_rich_higher_tm(self, thermo: Thermodynamics) -> None:
        at_rich = "ATATATATATATAT"
        gc_rich = "GCGCGCGCGCGCGC"
        assert thermo.calculate_tm(gc_rich) > thermo.calculate_tm(at_rich)

    def test_lowercase_accepted(self, thermo: Thermodynamics) -> None:
        tm = thermo.calculate_tm(PRIMER_LEFT.lower())
        assert isinstance(tm, float)

    def test_custom_params_accepted(self, thermo: Thermodynamics) -> None:
        tm = thermo.calculate_tm(PRIMER_LEFT, params={"mv_conc": 100.0})
        assert isinstance(tm, float)


# ---------------------------------------------------------------------------
# calculate_hairpin_dg
# ---------------------------------------------------------------------------


class TestCalculateHairpinDg:
    def test_returns_float(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_hairpin_dg(PRIMER_LEFT)
        assert isinstance(dg, float)

    def test_units_kcal_per_mol(self, thermo: Thermodynamics) -> None:
        # primer3 returns cal/mol → we convert to kcal/mol; expected range -20 to +5
        dg = thermo.calculate_hairpin_dg(PRIMER_LEFT)
        assert -50.0 < dg < 10.0

    def test_self_complementary_negative(self, thermo: Thermodynamics) -> None:
        # A self-complementary sequence should have strong hairpin
        seq = "GCGCAAAAAGCGC"  # short stem-loop
        dg = thermo.calculate_hairpin_dg(seq)
        assert dg < 0  # stable hairpin → negative dG


# ---------------------------------------------------------------------------
# calculate_dimer_dg
# ---------------------------------------------------------------------------


class TestCalculateDimerDg:
    def test_returns_float(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_dimer_dg(PRIMER_LEFT, PRIMER_RIGHT)
        assert isinstance(dg, float)

    def test_self_dimer_same_as_homodimer(self, thermo: Thermodynamics) -> None:
        dg_self = thermo.calculate_dimer_dg(PRIMER_LEFT, PRIMER_LEFT)
        assert isinstance(dg_self, float)

    def test_complementary_stronger_dimer(self, thermo: Thermodynamics) -> None:
        # Perfectly complementary pair should have more negative ΔG
        seq = "ATCGATCGATCG"
        complement = "CGATCGATCGAT"
        dg_comp = thermo.calculate_dimer_dg(seq, complement)
        dg_non = thermo.calculate_dimer_dg(seq, "TTTTTTTTTTTT")
        assert dg_comp < dg_non


# ---------------------------------------------------------------------------
# check_gc_clamp
# ---------------------------------------------------------------------------


class TestCheckGcClamp:
    def test_good_clamp(self, thermo: Thermodynamics) -> None:
        # Ends with GC — should pass
        assert thermo.check_gc_clamp("ATCGATCGATCG") is True

    def test_no_gc_at_end_fails(self, thermo: Thermodynamics) -> None:
        assert thermo.check_gc_clamp("GCGCGCTATAT") is False

    def test_too_many_gc_fails(self, thermo: Thermodynamics) -> None:
        # 4 GC in last 5 → fails max_gc=3
        assert thermo.check_gc_clamp("ATCGGGGGG", window=5, min_gc=1, max_gc=3) is False

    def test_empty_seq_fails(self, thermo: Thermodynamics) -> None:
        assert thermo.check_gc_clamp("") is False

    def test_exactly_min_gc(self, thermo: Thermodynamics) -> None:
        # Last 5 bases "ATCGX" where last base is C — exactly 1 GC
        assert thermo.check_gc_clamp("ATATATATC", window=5, min_gc=1, max_gc=3) is True


# ---------------------------------------------------------------------------
# check_repeat_runs
# ---------------------------------------------------------------------------


class TestCheckRepeatRuns:
    def test_no_repeats(self, thermo: Thermodynamics) -> None:
        assert thermo.check_repeat_runs("ATCGATCGATCG") is True

    def test_long_mono_repeat_fails(self, thermo: Thermodynamics) -> None:
        # 5 identical bases → exceeds max_run=4
        assert thermo.check_repeat_runs("ATCGAAAAATCG") is False

    def test_exactly_max_run_passes(self, thermo: Thermodynamics) -> None:
        # Exactly 4 identical bases — OK
        assert thermo.check_repeat_runs("ATCGAAAATCG", max_run=4) is True

    def test_dinucleotide_repeat_fails(self, thermo: Thermodynamics) -> None:
        # ATATATATAT — 5 AT repeats; dinucleotide check with max_run=4
        assert thermo.check_repeat_runs("ATATATATATAT", max_run=4) is False

    def test_custom_max_run(self, thermo: Thermodynamics) -> None:
        assert thermo.check_repeat_runs("AAAAAAA", max_run=8) is True
        assert thermo.check_repeat_runs("AAAAAAA", max_run=6) is False


# ---------------------------------------------------------------------------
# full_thermodynamic_profile
# ---------------------------------------------------------------------------


class TestFullThermodynamicProfile:
    def test_returns_thermo_profile(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile, ThermoProfile)

    def test_sequence_preserved(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert profile.sequence == PRIMER_LEFT

    def test_tm_reasonable(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert 40.0 <= profile.tm <= 80.0

    def test_gc_percent_range(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert 0.0 <= profile.gc_percent <= 100.0

    def test_hairpin_dg_is_float(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile.hairpin_dg, float)

    def test_pass_all_good_primer(self, thermo: Thermodynamics) -> None:
        good = Oligo(sequence="ATCGATCGATCGATCGATCG")
        profile = thermo.full_thermodynamic_profile(good)
        assert isinstance(profile.pass_all, bool)

    def test_repeat_seq_fails_repeat_ok(self, thermo: Thermodynamics) -> None:
        repeat_oligo = Oligo(sequence=REPEAT_SEQ)
        profile = thermo.full_thermodynamic_profile(repeat_oligo)
        assert profile.repeat_ok is False
        assert profile.pass_all is False

    def test_no_gc_clamp_fails(self, thermo: Thermodynamics) -> None:
        oligo = Oligo(sequence=LOW_GC_SEQ)
        profile = thermo.full_thermodynamic_profile(oligo)
        assert profile.gc_clamp_ok is False


# ---------------------------------------------------------------------------
# validate_all
# ---------------------------------------------------------------------------


class TestValidateAll:
    def test_returns_list(self, thermo: Thermodynamics, primer_pair: PrimerPair) -> None:
        result = thermo.validate_all([primer_pair])
        assert isinstance(result, list)
        assert len(result) == 1

    def test_oligo_tm_updated(self, thermo: Thermodynamics, primer_pair: PrimerPair) -> None:
        thermo.validate_all([primer_pair])
        assert primer_pair.left_primer.tm > 0

    def test_oligo_gc_percent_updated(
        self, thermo: Thermodynamics, primer_pair: PrimerPair
    ) -> None:
        thermo.validate_all([primer_pair])
        assert primer_pair.left_primer.gc_percent > 0

    def test_empty_list(self, thermo: Thermodynamics) -> None:
        result = thermo.validate_all([])
        assert result == []

    def test_probe_annotated(self, thermo: Thermodynamics) -> None:
        left = Oligo(sequence=PRIMER_LEFT)
        right = Oligo(sequence=PRIMER_RIGHT)
        probe = Oligo(sequence=PROBE_SEQ)
        pair = PrimerPair(pair_id="p0", left_primer=left, right_primer=right, probe=probe)
        thermo.validate_all([pair])
        assert pair.probe is not None
        assert pair.probe.tm > 0


# ---------------------------------------------------------------------------
# _calc_end_stability helper
# ---------------------------------------------------------------------------


class TestCalcEndStability:
    def test_returns_float(self) -> None:
        result = _calc_end_stability(PRIMER_LEFT)
        assert isinstance(result, float)
