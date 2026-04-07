"""Tests cho PrimerEvaluator — kiểm thử nhiệt động học và BLAST specificity."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.modules.blast_specificity import BlastHit, OffTargetAmplicon
from tta_primer_design.modules.primer_evaluator import (
    EvaluationReport,
    PrimerEvaluator,
    ThermoProfile,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Primer pair thực tế cho ACTB (beta-actin) — qPCR primers đã được kiểm tra
LEFT_SEQ = "GCACTGACCTCCCACTTCAA"
RIGHT_SEQ = "TTGCTGATCCACATCTGCTG"
# Primer với GC% thấp, Tm thấp — dùng để test FAIL case
LOW_GC_LEFT = "ATATATATAT"
LOW_GC_RIGHT = "TATATATATA"
# Primer có repeat run dài
REPEAT_SEQ = "AAAAATCGATCGATCGATCG"


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def evaluator(config: AppConfig) -> PrimerEvaluator:
    return PrimerEvaluator(config)


# ---------------------------------------------------------------------------
# Helpers: mock BLAST để tránh gọi mạng trong unit test
# ---------------------------------------------------------------------------


def _make_blast_hit(
    subject_id: str,
    sbjct_start: int,
    sbjct_end: int,
    identity: float = 100.0,
    evalue: float = 1e-5,
) -> BlastHit:
    """Tạo một BlastHit giả cho kiểm thử."""
    return BlastHit(
        subject_id=subject_id,
        subject_title=f"Mock hit {subject_id}",
        identity=identity,
        alignment_length=20,
        mismatches=int(20 * (1 - identity / 100)),
        gaps=0,
        query_start=1,
        query_end=20,
        subject_start=sbjct_start,
        subject_end=sbjct_end,
        evalue=evalue,
        bit_score=40.0,
        mismatches_3prime=0,
    )


# ---------------------------------------------------------------------------
# Test: _validate_sequence
# ---------------------------------------------------------------------------


class TestValidateSequence:
    """Kiểm tra hàm validate và chuẩn hoá chuỗi nucleotide."""

    def test_valid_uppercase(self, evaluator: PrimerEvaluator) -> None:
        result = evaluator._validate_sequence("ATCGATCG", "left")
        assert result == "ATCGATCG"

    def test_lowercase_converted(self, evaluator: PrimerEvaluator) -> None:
        result = evaluator._validate_sequence("atcgatcg", "left")
        assert result == "ATCGATCG"

    def test_strip_whitespace(self, evaluator: PrimerEvaluator) -> None:
        result = evaluator._validate_sequence("  ATCG  ", "left")
        assert result == "ATCG"

    def test_empty_raises(self, evaluator: PrimerEvaluator) -> None:
        with pytest.raises(ValueError, match="rỗng"):
            evaluator._validate_sequence("", "left")

    def test_invalid_chars_raises(self, evaluator: PrimerEvaluator) -> None:
        with pytest.raises(ValueError, match="không hợp lệ"):
            evaluator._validate_sequence("ATCG123", "left")

    def test_iupac_codes_accepted(self, evaluator: PrimerEvaluator) -> None:
        # IUPAC ambiguity codes phải được chấp nhận
        result = evaluator._validate_sequence("ATCGRYN", "left")
        assert result == "ATCGRYN"


# ---------------------------------------------------------------------------
# Test: _calc_gc_percent
# ---------------------------------------------------------------------------


class TestCalcGCPercent:
    """Kiểm tra tính tỷ lệ GC."""

    def test_pure_gc(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._calc_gc_percent("GCGCGC") == pytest.approx(100.0)

    def test_pure_at(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._calc_gc_percent("ATATAT") == pytest.approx(0.0)

    def test_mixed(self, evaluator: PrimerEvaluator) -> None:
        # GCATGCAT → 4 GC / 8 total = 50%
        assert evaluator._calc_gc_percent("GCATGCAT") == pytest.approx(50.0)

    def test_empty_returns_zero(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._calc_gc_percent("") == pytest.approx(0.0)

    def test_typical_primer(self, evaluator: PrimerEvaluator) -> None:
        # LEFT_SEQ = "GCACTGACCTCCCACTTCAA"
        # G,C,A,C,T,G,A,C,C,T,C,C,C,A,C,T,T,C,A,A → 10 GC / 20 = 50%
        gc = evaluator._calc_gc_percent(LEFT_SEQ)
        assert 40.0 <= gc <= 70.0


# ---------------------------------------------------------------------------
# Test: _check_gc_clamp
# ---------------------------------------------------------------------------


class TestCheckGCClamp:
    """Kiểm tra GC clamp ở 3' end."""

    def test_good_clamp(self, evaluator: PrimerEvaluator) -> None:
        # Kết thúc bằng "GCATG" — 2 GC trong 5 bp cuối
        assert evaluator._check_gc_clamp("ATCGATCGCATG") is True

    def test_no_gc_at_end(self, evaluator: PrimerEvaluator) -> None:
        # Kết thúc bằng "AAAAT" — 0 GC → FAIL
        assert evaluator._check_gc_clamp("GCGCGCAAAAT") is False

    def test_too_many_gc_at_end(self, evaluator: PrimerEvaluator) -> None:
        # Kết thúc bằng "GCGCG" — 5 GC → FAIL (max=3)
        assert evaluator._check_gc_clamp("ATCGCGCGCGCGCGCG") is False

    def test_exactly_min_gc(self, evaluator: PrimerEvaluator) -> None:
        # 1 GC ở cuối 5 bp → PASS
        assert evaluator._check_gc_clamp("ATCGATCGATCG") is True

    def test_empty_sequence(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._check_gc_clamp("") is False


# ---------------------------------------------------------------------------
# Test: _check_repeat_runs
# ---------------------------------------------------------------------------


class TestCheckRepeatRuns:
    """Kiểm tra repeat run."""

    def test_no_repeat(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._check_repeat_runs("ATCGATCG") is True

    def test_run_of_4_ok(self, evaluator: PrimerEvaluator) -> None:
        # Chính xác 4 → OK (max_run=4 nghĩa là phải > 4 mới fail)
        assert evaluator._check_repeat_runs("AAAATCG") is True

    def test_run_of_5_fail(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._check_repeat_runs("AAAAATCG") is False

    def test_long_run_fail(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._check_repeat_runs(REPEAT_SEQ) is False

    def test_empty_sequence(self, evaluator: PrimerEvaluator) -> None:
        assert evaluator._check_repeat_runs("") is True


# ---------------------------------------------------------------------------
# Test: _calc_thermodynamics (mocked primer3)
# ---------------------------------------------------------------------------


class TestCalcThermodynamics:
    """Kiểm tra tính nhiệt động học dùng primer3-py (mock để tách biệt)."""

    def _mock_thermo_result(self, dg: float = -2000.0, tm: float = 40.0) -> MagicMock:
        """Tạo ThermoResult mock với .dg và .tm."""
        mock = MagicMock()
        mock.dg = dg  # cal/mol
        mock.tm = tm
        return mock

    def test_returns_thermo_profile(self, evaluator: PrimerEvaluator) -> None:
        mock_result = self._mock_thermo_result()
        with patch("primer3.calc_tm", return_value=60.0), patch(
            "primer3.calc_hairpin", return_value=mock_result
        ), patch("primer3.calc_homodimer", return_value=mock_result):
            profile = evaluator._calc_thermodynamics(LEFT_SEQ)
        assert isinstance(profile, ThermoProfile)

    def test_tm_assigned(self, evaluator: PrimerEvaluator) -> None:
        mock_result = self._mock_thermo_result()
        with patch("primer3.calc_tm", return_value=60.0), patch(
            "primer3.calc_hairpin", return_value=mock_result
        ), patch("primer3.calc_homodimer", return_value=mock_result):
            profile = evaluator._calc_thermodynamics(LEFT_SEQ)
        assert profile.tm == pytest.approx(60.0)

    def test_dg_converted_to_kcal(self, evaluator: PrimerEvaluator) -> None:
        # dg = -5000 cal/mol → -5.0 kcal/mol
        mock_result = self._mock_thermo_result(dg=-5000.0)
        with patch("primer3.calc_tm", return_value=60.0), patch(
            "primer3.calc_hairpin", return_value=mock_result
        ), patch("primer3.calc_homodimer", return_value=mock_result):
            profile = evaluator._calc_thermodynamics(LEFT_SEQ)
        assert profile.hairpin_dg == pytest.approx(-5.0)

    def test_pass_all_true_when_in_range(self, evaluator: PrimerEvaluator) -> None:
        mock_result = self._mock_thermo_result(dg=-2000.0)  # -2 kcal/mol → OK
        with patch("primer3.calc_tm", return_value=60.0), patch(
            "primer3.calc_hairpin", return_value=mock_result
        ), patch("primer3.calc_homodimer", return_value=mock_result):
            profile = evaluator._calc_thermodynamics("GCACTGACCTCCCACTTCAA")
        assert profile.pass_all is True

    def test_fail_when_tm_out_of_range(self, evaluator: PrimerEvaluator) -> None:
        mock_result = self._mock_thermo_result(dg=-2000.0)
        with patch("primer3.calc_tm", return_value=45.0), patch(  # Tm quá thấp
            "primer3.calc_hairpin", return_value=mock_result
        ), patch("primer3.calc_homodimer", return_value=mock_result):
            profile = evaluator._calc_thermodynamics(LEFT_SEQ)
        assert profile.pass_all is False
        assert any("Tm" in w for w in profile.warnings)

    def test_fail_when_hairpin_too_strong(self, evaluator: PrimerEvaluator) -> None:
        mock_result = self._mock_thermo_result(dg=-15000.0)  # -15 kcal/mol → FAIL
        with patch("primer3.calc_tm", return_value=60.0), patch(
            "primer3.calc_hairpin", return_value=mock_result
        ), patch("primer3.calc_homodimer", return_value=self._mock_thermo_result()):
            profile = evaluator._calc_thermodynamics(LEFT_SEQ)
        assert profile.pass_all is False
        assert any("Hairpin" in w for w in profile.warnings)

    def test_probe_uses_different_tm_threshold(self, evaluator: PrimerEvaluator) -> None:
        mock_result = self._mock_thermo_result(dg=-2000.0)
        # Tm=60°C → OK cho primer nhưng FAIL cho probe (cần 65–72°C)
        with patch("primer3.calc_tm", return_value=60.0), patch(
            "primer3.calc_hairpin", return_value=mock_result
        ), patch("primer3.calc_homodimer", return_value=mock_result):
            profile_probe = evaluator._calc_thermodynamics(LEFT_SEQ, is_probe=True)
        assert profile_probe.pass_all is False
        assert any("probe" in w.lower() for w in profile_probe.warnings)


# ---------------------------------------------------------------------------
# Test: _predict_offtarget_amplicons
# ---------------------------------------------------------------------------


class TestPredictOfftargetAmplicons:
    """Kiểm tra dự đoán off-target amplicon từ BLAST hits."""

    def test_no_hits_returns_empty(self, evaluator: PrimerEvaluator) -> None:
        result = evaluator._predict_offtarget_amplicons([], [])
        assert result == []

    def test_left_plus_right_minus_same_subject_small_size(
        self, evaluator: PrimerEvaluator
    ) -> None:
        # Left: + strand tại 100–119, Right: - strand tại 300–281 (size ≈ 201)
        lh = _make_blast_hit("NM_001101", sbjct_start=100, sbjct_end=119)
        rh = _make_blast_hit("NM_001101", sbjct_start=300, sbjct_end=281)
        result = evaluator._predict_offtarget_amplicons([lh], [rh])
        assert len(result) == 1
        assert result[0].subject_id == "NM_001101"
        assert 0 < result[0].amplicon_size <= evaluator.max_amplicon_size

    def test_different_subjects_no_amplicon(self, evaluator: PrimerEvaluator) -> None:
        lh = _make_blast_hit("NM_001101", sbjct_start=100, sbjct_end=119)
        rh = _make_blast_hit("NM_999999", sbjct_start=300, sbjct_end=281)
        result = evaluator._predict_offtarget_amplicons([lh], [rh])
        assert result == []

    def test_amplicon_too_large_excluded(self, evaluator: PrimerEvaluator) -> None:
        # Khoảng cách lớn hơn max_amplicon_size (default 4000)
        lh = _make_blast_hit("NM_001101", sbjct_start=100, sbjct_end=119)
        rh = _make_blast_hit("NM_001101", sbjct_start=10000, sbjct_end=9981)
        result = evaluator._predict_offtarget_amplicons([lh], [rh])
        assert result == []

    def test_both_plus_strand_no_amplicon(self, evaluator: PrimerEvaluator) -> None:
        # Cả hai đều + strand → không tạo amplicon
        lh = _make_blast_hit("NM_001101", sbjct_start=100, sbjct_end=119)
        rh = _make_blast_hit("NM_001101", sbjct_start=200, sbjct_end=219)
        result = evaluator._predict_offtarget_amplicons([lh], [rh])
        assert result == []

    def test_overlap_no_amplicon(self, evaluator: PrimerEvaluator) -> None:
        # Right primer ở trước left primer → overlap → không tạo amplicon
        lh = _make_blast_hit("NM_001101", sbjct_start=200, sbjct_end=219)
        rh = _make_blast_hit("NM_001101", sbjct_start=150, sbjct_end=131)
        result = evaluator._predict_offtarget_amplicons([lh], [rh])
        assert result == []

    def test_multiple_off_targets(self, evaluator: PrimerEvaluator) -> None:
        lh = _make_blast_hit("NM_001", sbjct_start=100, sbjct_end=119)
        rh1 = _make_blast_hit("NM_001", sbjct_start=300, sbjct_end=281)
        rh2 = _make_blast_hit("NM_001", sbjct_start=500, sbjct_end=481)
        result = evaluator._predict_offtarget_amplicons([lh], [rh1, rh2])
        assert len(result) == 2


# ---------------------------------------------------------------------------
# Test: evaluate() — mock cả primer3 và BLAST
# ---------------------------------------------------------------------------


class TestEvaluate:
    """Kiểm tra luồng đánh giá tổng thể với mock."""

    def _setup_mocks(self, tm: float = 60.0, dg: float = -2000.0) -> dict:
        """Tạo mock context managers cho primer3 và BLAST."""
        mock_thermo = MagicMock()
        mock_thermo.dg = dg

        patches = {
            "primer3.calc_tm": patch("primer3.calc_tm", return_value=tm),
            "primer3.calc_hairpin": patch("primer3.calc_hairpin", return_value=mock_thermo),
            "primer3.calc_homodimer": patch("primer3.calc_homodimer", return_value=mock_thermo),
            "primer3.calc_heterodimer": patch(
                "primer3.calc_heterodimer", return_value=mock_thermo
            ),
        }
        return patches

    def test_evaluate_returns_evaluation_report(self, evaluator: PrimerEvaluator) -> None:
        patches = self._setup_mocks()
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq=LEFT_SEQ,
                right_seq=RIGHT_SEQ,
                run_blast=False,
            )
        assert isinstance(report, EvaluationReport)

    def test_evaluate_no_blast(self, evaluator: PrimerEvaluator) -> None:
        patches = self._setup_mocks()
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq=LEFT_SEQ,
                right_seq=RIGHT_SEQ,
                run_blast=False,
            )
        assert report.blast_performed is False
        assert report.specificity is None

    def test_evaluate_pair_id_preserved(self, evaluator: PrimerEvaluator) -> None:
        patches = self._setup_mocks()
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            report = evaluator.evaluate(
                pair_id="MY_PAIR_123",
                left_seq=LEFT_SEQ,
                right_seq=RIGHT_SEQ,
                run_blast=False,
            )
        assert report.pair_id == "MY_PAIR_123"

    def test_evaluate_sequences_normalized(self, evaluator: PrimerEvaluator) -> None:
        patches = self._setup_mocks()
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq=LEFT_SEQ.lower(),
                right_seq=RIGHT_SEQ.lower(),
                run_blast=False,
            )
        assert report.left_sequence == LEFT_SEQ.upper()
        assert report.right_sequence == RIGHT_SEQ.upper()

    def test_evaluate_with_probe(self, evaluator: PrimerEvaluator) -> None:
        patches = self._setup_mocks(tm=68.0)
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq=LEFT_SEQ,
                right_seq=RIGHT_SEQ,
                probe_seq="AATTTCCCGGGAAATTTCCC",
                run_blast=False,
            )
        assert report.probe_thermo is not None
        assert report.probe_sequence == "AATTTCCCGGGAAATTTCCC"

    def test_evaluate_pass_recommendation_good_primers(
        self, evaluator: PrimerEvaluator
    ) -> None:
        # Tm=60°C, ΔG=-2 kcal/mol → tất cả đều trong ngưỡng
        patches = self._setup_mocks(tm=60.0, dg=-2000.0)
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq="GCACTGACCTCCCACTTCAA",  # 50% GC, GC clamp OK
                right_seq="TTGCTGATCCACATCTGCTG",  # 50% GC, GC clamp OK
                run_blast=False,
            )
        assert report.overall_recommendation == "PASS"

    def test_evaluate_fail_recommendation_bad_primers(
        self, evaluator: PrimerEvaluator
    ) -> None:
        # Tm=45°C → dưới ngưỡng → FAIL
        patches = self._setup_mocks(tm=45.0, dg=-2000.0)
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq=LEFT_SEQ,
                right_seq=RIGHT_SEQ,
                run_blast=False,
            )
        assert report.overall_recommendation == "FAIL"
        assert report.pass_thermodynamics is False

    def test_evaluate_invalid_sequence_raises(self, evaluator: PrimerEvaluator) -> None:
        with pytest.raises(ValueError):
            evaluator.evaluate(
                pair_id="bad",
                left_seq="ATCG123XYZ",
                right_seq=RIGHT_SEQ,
                run_blast=False,
            )

    def test_evaluate_blast_failure_handled_gracefully(
        self, evaluator: PrimerEvaluator
    ) -> None:
        """Nếu BLAST thất bại, báo cáo vẫn trả về (không crash)."""
        patches = self._setup_mocks(tm=60.0, dg=-2000.0)
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
            patch.object(
                evaluator,
                "_blast_pair",
                side_effect=RuntimeError("Network error"),
            ),
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq="GCACTGACCTCCCACTTCAA",
                right_seq="TTGCTGATCCACATCTGCTG",
                run_blast=True,
            )
        assert report.blast_performed is False
        assert report.specificity is None

    def test_evaluate_blast_with_offtargets_is_fail(
        self, evaluator: PrimerEvaluator
    ) -> None:
        """Nếu có off-target amplicon → FAIL."""
        mock_specificity = MagicMock()
        mock_specificity.off_target_amplicons = [
            OffTargetAmplicon(
                subject_id="NM_999",
                amplicon_size=300,
                left_hit=_make_blast_hit("NM_999", 100, 119),
                right_hit=_make_blast_hit("NM_999", 300, 281),
            )
        ]
        mock_specificity.specificity_score = 80.0
        mock_specificity.is_specific = False

        patches = self._setup_mocks(tm=60.0, dg=-2000.0)
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
            patch.object(evaluator, "_blast_pair", return_value=mock_specificity),
        ):
            report = evaluator.evaluate(
                pair_id="test",
                left_seq="GCACTGACCTCCCACTTCAA",
                right_seq="TTGCTGATCCACATCTGCTG",
                run_blast=True,
            )
        assert report.overall_recommendation == "FAIL"
        assert report.blast_performed is True


# ---------------------------------------------------------------------------
# Test: CLI evaluate command
# ---------------------------------------------------------------------------


class TestCLIEvaluateCommand:
    """Kiểm tra CLI command ``evaluate``."""

    def _setup_mocks(self, tm: float = 60.0, dg: float = -2000.0) -> dict:
        mock_thermo = MagicMock()
        mock_thermo.dg = dg
        return {
            "primer3.calc_tm": patch("primer3.calc_tm", return_value=tm),
            "primer3.calc_hairpin": patch("primer3.calc_hairpin", return_value=mock_thermo),
            "primer3.calc_homodimer": patch("primer3.calc_homodimer", return_value=mock_thermo),
            "primer3.calc_heterodimer": patch(
                "primer3.calc_heterodimer", return_value=mock_thermo
            ),
        }

    def test_evaluate_help(self) -> None:
        from click.testing import CliRunner

        from tta_primer_design.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["evaluate", "--help"])
        assert result.exit_code == 0
        assert "--left-primer" in result.output
        assert "--right-primer" in result.output
        assert "--no-blast" in result.output

    def test_evaluate_missing_left_primer(self) -> None:
        from click.testing import CliRunner

        from tta_primer_design.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["evaluate", "--right-primer", RIGHT_SEQ])
        assert result.exit_code != 0

    def test_evaluate_no_blast_exits_0_for_good_primers(self) -> None:
        from click.testing import CliRunner

        from tta_primer_design.cli import main

        runner = CliRunner()
        patches = self._setup_mocks(tm=60.0, dg=-2000.0)
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            result = runner.invoke(
                main,
                [
                    "evaluate",
                    "--left-primer",
                    "GCACTGACCTCCCACTTCAA",
                    "--right-primer",
                    "TTGCTGATCCACATCTGCTG",
                    "--no-blast",
                ],
            )
        assert result.exit_code == 0
        assert "PASS" in result.output
        assert "FAIL" not in result.output

    def test_evaluate_bad_primers_exits_1(self) -> None:
        from click.testing import CliRunner

        from tta_primer_design.cli import main

        runner = CliRunner()
        patches = self._setup_mocks(tm=40.0, dg=-2000.0)  # Tm quá thấp
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            result = runner.invoke(
                main,
                [
                    "evaluate",
                    "--left-primer",
                    LEFT_SEQ,
                    "--right-primer",
                    RIGHT_SEQ,
                    "--no-blast",
                ],
            )
        assert result.exit_code == 1
        assert "FAIL" in result.output

    def test_evaluate_saves_json_to_output(self, tmp_path) -> None:
        import json

        from click.testing import CliRunner

        from tta_primer_design.cli import main

        runner = CliRunner()
        patches = self._setup_mocks(tm=60.0, dg=-2000.0)
        with (
            patches["primer3.calc_tm"],
            patches["primer3.calc_hairpin"],
            patches["primer3.calc_homodimer"],
            patches["primer3.calc_heterodimer"],
        ):
            result = runner.invoke(
                main,
                [
                    "evaluate",
                    "--left-primer",
                    "GCACTGACCTCCCACTTCAA",
                    "--right-primer",
                    "TTGCTGATCCACATCTGCTG",
                    "--pair-id",
                    "test_pair",
                    "--no-blast",
                    "--output",
                    str(tmp_path),
                ],
            )
        assert result.exit_code == 0
        out_file = tmp_path / "evaluation_test_pair.json"
        assert out_file.exists()
        data = json.loads(out_file.read_text())
        assert data["pair_id"] == "test_pair"

    def test_evaluate_invalid_sequence_exits_2(self) -> None:
        from click.testing import CliRunner

        from tta_primer_design.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "evaluate",
                "--left-primer",
                "ATCG123XYZ",
                "--right-primer",
                RIGHT_SEQ,
                "--no-blast",
            ],
        )
        assert result.exit_code == 2
