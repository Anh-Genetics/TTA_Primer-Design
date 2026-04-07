"""Tests cho CLI — tta-primer-design --help và run command."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from tta_primer_design.cli import main
from tta_primer_design.modules.blast_specificity import (
    BlastHit,
    SpecificityResult,
)


class TestCLIHelp:
    """Test --help output."""

    def test_main_help(self) -> None:
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "tta-primer-design" in result.output.lower() or "primer" in result.output.lower()

    def test_main_version(self) -> None:
        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.output

    def test_run_help(self) -> None:
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "--config" in result.output
        assert "--input" in result.output
        assert "--output" in result.output

    def test_run_help_shows_log_level(self) -> None:
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--help"])
        assert result.exit_code == 0
        assert "--log-level" in result.output


class TestCLIRunCommand:
    """Test lệnh ``run`` với file thực tế."""

    def test_run_missing_input_fails(self) -> None:
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--output", "/tmp/out"])
        assert result.exit_code != 0

    def test_run_missing_output_fails(self) -> None:
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--input", "dummy.csv"])
        assert result.exit_code != 0

    def test_run_nonexistent_input_exits_2(self, tmp_path) -> None:
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "run",
                "--input",
                str(tmp_path / "nonexistent.csv"),
                "--output",
                str(tmp_path / "out"),
            ],
        )
        # Click validates --input path exists; exits with error
        assert result.exit_code != 0


_LEFT = "GCAAGGAATGGTTTCAGAAATCCA"
_RIGHT = "CAGGACTCCATGTCGTCCA"
_PROBE = "TGCAGCCACACTTTCTACAATGAGC"


def _make_spec_result(pair_id: str = "pair_001") -> SpecificityResult:
    left_hit = BlastHit(subject_id="NM_001101", mismatches=0, evalue=1e-5)
    right_hit = BlastHit(subject_id="NM_001101", mismatches=0, evalue=1e-5)
    return SpecificityResult(
        primer_pair_id=pair_id,
        blast_hits_left=[left_hit],
        blast_hits_right=[right_hit],
        off_target_amplicons=[],
        specificity_score=100.0,
        is_specific=True,
    )


class TestCLIEvaluateCommand:
    """Tests for the ``evaluate`` CLI command."""

    def test_evaluate_help(self) -> None:
        runner = CliRunner()
        result = runner.invoke(main, ["evaluate", "--help"])
        assert result.exit_code == 0
        assert "--left" in result.output
        assert "--right" in result.output

    def test_evaluate_skip_blast_basic(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["evaluate", "--left", _LEFT, "--right", _RIGHT, "--skip-blast"],
        )
        assert result.exit_code == 0
        assert "Tm" in result.output or "tm" in result.output.lower()

    def test_evaluate_invalid_left_sequence(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["evaluate", "--left", "ATCGX123", "--right", _RIGHT, "--skip-blast"],
        )
        assert result.exit_code == 2

    def test_evaluate_invalid_right_sequence(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["evaluate", "--left", _LEFT, "--right", "ATCG!@#", "--skip-blast"],
        )
        assert result.exit_code == 2

    def test_evaluate_invalid_probe_sequence(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "evaluate",
                "--left",
                _LEFT,
                "--right",
                _RIGHT,
                "--probe",
                "ATCGXYZ",
                "--skip-blast",
            ],
        )
        assert result.exit_code == 2

    def test_evaluate_saves_json(self, tmp_path: Path) -> None:
        out = tmp_path / "eval_result.json"
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "evaluate",
                "--left",
                _LEFT,
                "--right",
                _RIGHT,
                "--skip-blast",
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0
        assert out.exists()
        data = json.loads(out.read_text())
        assert "left_primer" in data
        assert "right_primer" in data
        assert "score" in data

    def test_evaluate_skip_blast_with_probe(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "evaluate",
                "--left",
                _LEFT,
                "--right",
                _RIGHT,
                "--probe",
                _PROBE,
                "--skip-blast",
            ],
        )
        assert result.exit_code == 0

    def test_evaluate_plot_ascii(self) -> None:
        runner = CliRunner()
        spec = _make_spec_result()
        with patch(
            "tta_primer_design.modules.blast_specificity.BlastSpecificity.check_pair",
            return_value=spec,
        ):
            result = runner.invoke(
                main,
                [
                    "evaluate",
                    "--left",
                    _LEFT,
                    "--right",
                    _RIGHT,
                    "--plot",
                ],
            )
        # Should not crash; exit code may vary based on blast mock
        assert result.exit_code in (0, 1)
        # ASCII chart was triggered (contains "BLAST" or hit info)
        assert "BLAST" in result.output or "Hit" in result.output or "hit" in result.output

    def test_evaluate_plot_output_file(self, tmp_path: Path) -> None:
        plot_out = tmp_path / "blast_plot.png"
        runner = CliRunner()
        spec = _make_spec_result()
        with patch(
            "tta_primer_design.modules.blast_specificity.BlastSpecificity.check_pair",
            return_value=spec,
        ):
            result = runner.invoke(
                main,
                [
                    "evaluate",
                    "--left",
                    _LEFT,
                    "--right",
                    _RIGHT,
                    "--plot-output",
                    str(plot_out),
                ],
            )
        assert result.exit_code in (0, 1)
        assert plot_out.exists()
