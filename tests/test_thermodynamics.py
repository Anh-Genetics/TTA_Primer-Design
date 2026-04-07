"""Tests cho Thermodynamics module."""

from __future__ import annotations

import pytest

from tta_primer_design.config import AppConfig
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.thermodynamics import Thermodynamics, ThermoProfile

_LEFT_SEQ = "GCAAGGAATGGTTTCAGAAATCCA"
_RIGHT_SEQ = "CAGGACTCCATGTCGTCCA"
_REPEAT_FAIL = "AAAAAATCGATCG"
_NO_GC_CLAMP = "GCATATATAT"  # last 5 = "TATAT" → 0 GC → out of [1,3]
_GC_CLAMP_OK = "ACGATCGATCGATCG"  # last 5 = "GATCG" → 3 GC → in [1,3]


@pytest.fixture()
def thermo() -> Thermodynamics:
    return Thermodynamics(AppConfig())


@pytest.fixture()
def left_oligo() -> Oligo:
    return Oligo(sequence=_LEFT_SEQ)


@pytest.fixture()
def right_oligo() -> Oligo:
    return Oligo(sequence=_RIGHT_SEQ)


@pytest.fixture()
def primer_pair(left_oligo: Oligo, right_oligo: Oligo) -> PrimerPair:
    return PrimerPair(pair_id="test_pair", left_primer=left_oligo, right_primer=right_oligo)


class TestCalculateTm:
    """Tests for calculate_tm()."""

    def test_returns_float(self, thermo: Thermodynamics) -> None:
        tm = thermo.calculate_tm(_LEFT_SEQ)
        assert isinstance(tm, float)

    def test_reasonable_range_left_primer(self, thermo: Thermodynamics) -> None:
        tm = thermo.calculate_tm(_LEFT_SEQ)
        assert 50.0 <= tm <= 80.0

    def test_reasonable_range_right_primer(self, thermo: Thermodynamics) -> None:
        tm = thermo.calculate_tm(_RIGHT_SEQ)
        assert 50.0 <= tm <= 80.0

    def test_longer_sequence_higher_tm(self, thermo: Thermodynamics) -> None:
        short = "ATCGATCGATCGATCG"
        long_seq = "ATCGATCGATCGATCGATCGATCGATCGATCG"
        tm_short = thermo.calculate_tm(short)
        tm_long = thermo.calculate_tm(long_seq)
        assert tm_long > tm_short

    def test_lowercase_input(self, thermo: Thermodynamics) -> None:
        tm_upper = thermo.calculate_tm(_LEFT_SEQ.upper())
        tm_lower = thermo.calculate_tm(_LEFT_SEQ.lower())
        assert abs(tm_upper - tm_lower) < 0.01


class TestCalculateHairpinDg:
    """Tests for calculate_hairpin_dg()."""

    def test_returns_float(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_hairpin_dg(_LEFT_SEQ)
        assert isinstance(dg, float)

    def test_kcal_units_range(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_hairpin_dg(_LEFT_SEQ)
        # Reasonable ΔG range for typical oligos (kcal/mol)
        assert -50.0 < dg < 10.0

    def test_right_primer_returns_float(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_hairpin_dg(_RIGHT_SEQ)
        assert isinstance(dg, float)


class TestCalculateDimerDg:
    """Tests for calculate_dimer_dg()."""

    def test_returns_float(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_dimer_dg(_LEFT_SEQ, _RIGHT_SEQ)
        assert isinstance(dg, float)

    def test_symmetric_inputs_returns_float(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_dimer_dg(_LEFT_SEQ, _LEFT_SEQ)
        assert isinstance(dg, float)

    def test_kcal_units_range(self, thermo: Thermodynamics) -> None:
        dg = thermo.calculate_dimer_dg(_LEFT_SEQ, _RIGHT_SEQ)
        assert -100.0 < dg < 10.0


class TestCheckGcClamp:
    """Tests for check_gc_clamp()."""

    def test_returns_true_for_gc_ending(self, thermo: Thermodynamics) -> None:
        # last 5 of "ACGATCGATCGATCG" = "GATCG" → G, A, T, C, G = 3 GC → True
        assert thermo.check_gc_clamp(_GC_CLAMP_OK) is True

    def test_returns_false_for_no_gc_in_tail(self, thermo: Thermodynamics) -> None:
        # last 5 of "GCATATATAT" = "TATAT" → 0 GC → False (below min_gc=1)
        assert thermo.check_gc_clamp(_NO_GC_CLAMP) is False

    def test_returns_false_for_too_many_gc(self, thermo: Thermodynamics) -> None:
        # "ATCGCGCGCG" → last 5 = "GCGCG" → 5 GC → False (exceeds max_gc=3)
        assert thermo.check_gc_clamp("ATCGCGCGCG") is False

    def test_one_gc_at_end_passes(self, thermo: Thermodynamics) -> None:
        # "ATATATATATG" → last 5 = "ATATG" → 1 GC → True
        assert thermo.check_gc_clamp("ATATATATATG") is True

    def test_custom_window(self, thermo: Thermodynamics) -> None:
        result = thermo.check_gc_clamp("ATCGATCGATCG", window=3, min_gc=1, max_gc=2)
        assert isinstance(result, bool)


class TestCheckRepeatRuns:
    """Tests for check_repeat_runs()."""

    def test_returns_false_for_long_poly_a(self, thermo: Thermodynamics) -> None:
        # 6 A's → triggers pattern (.)\1{4,} → returns False
        assert thermo.check_repeat_runs(_REPEAT_FAIL) is False

    def test_returns_true_for_no_repeat(self, thermo: Thermodynamics) -> None:
        assert thermo.check_repeat_runs("ATCGATCG") is True

    def test_exactly_four_in_a_row_passes(self, thermo: Thermodynamics) -> None:
        # 4 A's = "AAAA" is borderline: pattern is {4,} meaning 5+ total (.)\1{4,}
        # "AAAA" → A followed by A{3} which is (.)\1{3} not {4} → passes
        assert thermo.check_repeat_runs("ATCGAAAATCG") is True

    def test_five_in_a_row_fails(self, thermo: Thermodynamics) -> None:
        assert thermo.check_repeat_runs("ATCGAAAAATCG") is False

    def test_poly_c_run_fails(self, thermo: Thermodynamics) -> None:
        assert thermo.check_repeat_runs("ATGCCCCCATG") is False

    def test_left_primer_passes(self, thermo: Thermodynamics) -> None:
        assert thermo.check_repeat_runs(_LEFT_SEQ) is True


class TestFullThermodynamicProfile:
    """Tests for full_thermodynamic_profile()."""

    def test_returns_thermo_profile(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile, ThermoProfile)

    def test_sequence_preserved(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert profile.sequence == _LEFT_SEQ

    def test_tm_populated(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert profile.tm > 0.0

    def test_gc_percent_populated(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert 0.0 < profile.gc_percent <= 100.0

    def test_hairpin_dg_is_float(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile.hairpin_dg, float)

    def test_self_dimer_dg_is_float(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile.self_dimer_dg, float)

    def test_end_stability_dg_is_float(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile.end_stability_dg, float)

    def test_gc_clamp_ok_is_bool(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile.gc_clamp_ok, bool)

    def test_repeat_ok_is_bool(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile.repeat_ok, bool)

    def test_pass_all_is_bool(self, thermo: Thermodynamics, left_oligo: Oligo) -> None:
        profile = thermo.full_thermodynamic_profile(left_oligo)
        assert isinstance(profile.pass_all, bool)

    def test_repeat_fail_sequence(self, thermo: Thermodynamics) -> None:
        oligo = Oligo(sequence=_REPEAT_FAIL)
        profile = thermo.full_thermodynamic_profile(oligo)
        assert profile.repeat_ok is False


class TestValidateAll:
    """Tests for validate_all()."""

    def test_returns_list(self, thermo: Thermodynamics, primer_pair: PrimerPair) -> None:
        result = thermo.validate_all([primer_pair])
        assert isinstance(result, list)
        assert len(result) == 1

    def test_updates_tm_on_left_primer(
        self, thermo: Thermodynamics, primer_pair: PrimerPair
    ) -> None:
        assert primer_pair.left_primer.tm == 0.0
        thermo.validate_all([primer_pair])
        assert primer_pair.left_primer.tm > 0.0

    def test_updates_tm_on_right_primer(
        self, thermo: Thermodynamics, primer_pair: PrimerPair
    ) -> None:
        thermo.validate_all([primer_pair])
        assert primer_pair.right_primer.tm > 0.0

    def test_updates_gc_percent(self, thermo: Thermodynamics, primer_pair: PrimerPair) -> None:
        thermo.validate_all([primer_pair])
        assert primer_pair.left_primer.gc_percent > 0.0

    def test_handles_empty_list(self, thermo: Thermodynamics) -> None:
        result = thermo.validate_all([])
        assert result == []

    def test_updates_probe_when_present(
        self, thermo: Thermodynamics, primer_pair: PrimerPair
    ) -> None:
        primer_pair.probe = Oligo(sequence="TGCAGCCACACTTTCTACAATGAGC")
        thermo.validate_all([primer_pair])
        assert primer_pair.probe.tm > 0.0

    def test_multiple_pairs(
        self, thermo: Thermodynamics, left_oligo: Oligo, right_oligo: Oligo
    ) -> None:
        pair1 = PrimerPair(
            pair_id="p1",
            left_primer=Oligo(sequence=_LEFT_SEQ),
            right_primer=Oligo(sequence=_RIGHT_SEQ),
        )
        pair2 = PrimerPair(
            pair_id="p2",
            left_primer=Oligo(sequence=_LEFT_SEQ),
            right_primer=Oligo(sequence=_RIGHT_SEQ),
        )
        result = thermo.validate_all([pair1, pair2])
        assert len(result) == 2
        for p in result:
            assert p.left_primer.tm > 0.0
