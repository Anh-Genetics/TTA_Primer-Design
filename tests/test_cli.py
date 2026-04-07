"""Tests cho CLI — tta-primer-design --help và run command."""

from __future__ import annotations

from click.testing import CliRunner

from tta_primer_design.cli import main


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
