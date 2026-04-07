"""Tests cho ReportGenerator module."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tta_primer_design.config import AppConfig, OutputConfig
from tta_primer_design.models import DesignResult, DesignTarget, Oligo, PrimerPair
from tta_primer_design.modules.report_generator import ReportGenerator

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_oligo(seq: str, tm: float = 60.0, gc: float = 50.0) -> Oligo:
    return Oligo(sequence=seq, tm=tm, gc_percent=gc)


def _make_pair(pair_id: str = "pair_01", include_probe: bool = False) -> PrimerPair:
    pair = PrimerPair(
        pair_id=pair_id,
        left_primer=_make_oligo("GCAAGGAATGGTTTCAGAAATCCA", tm=61.5, gc=41.7),
        right_primer=_make_oligo("CAGGACTCCATGTCGTCCA", tm=59.8, gc=57.9),
        amplicon_size=120,
        pair_penalty=0.5,
        score=78.0,
    )
    if include_probe:
        pair.probe = _make_oligo("TGCAGCCACACTTTCTACAATGAGC", tm=68.0, gc=52.0)
    return pair


def _make_target(target_id: str = "ACTB") -> DesignTarget:
    return DesignTarget(target_id=target_id, input_type="sequence")


def _make_result(
    target_id: str = "ACTB",
    n_pairs: int = 2,
    include_probe: bool = False,
    status: str = "success",
) -> DesignResult:
    target = _make_target(target_id)
    pairs = [_make_pair(f"pair_{i:02d}", include_probe=include_probe) for i in range(n_pairs)]
    return DesignResult(target=target, primer_pairs=pairs, status=status)


def _make_empty_result(target_id: str = "GENE_FAIL") -> DesignResult:
    return DesignResult(
        target=_make_target(target_id),
        primer_pairs=[],
        status="no_primers",
    )


@pytest.fixture()
def results() -> list[DesignResult]:
    return [_make_result("ACTB"), _make_result("GAPDH", include_probe=True)]


@pytest.fixture()
def empty_results() -> list[DesignResult]:
    return [_make_empty_result("GENE_FAIL")]


def _make_reporter(formats: list[str]) -> ReportGenerator:
    cfg = AppConfig()
    cfg.output = OutputConfig(formats=formats)
    return ReportGenerator(cfg)


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------


class TestGenerateCsv:
    """Tests for generate_csv()."""

    def test_creates_file(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["csv"])
        out = tmp_path / "summary.csv"
        reporter.generate_csv(results, out)
        assert out.exists()

    def test_csv_contains_target_id_column(
        self, tmp_path: Path, results: list[DesignResult]
    ) -> None:
        reporter = _make_reporter(["csv"])
        out = tmp_path / "summary.csv"
        reporter.generate_csv(results, out)
        content = out.read_text()
        assert "target_id" in content

    def test_csv_contains_left_primer_column(
        self, tmp_path: Path, results: list[DesignResult]
    ) -> None:
        reporter = _make_reporter(["csv"])
        out = tmp_path / "summary.csv"
        reporter.generate_csv(results, out)
        content = out.read_text()
        assert "left_primer" in content

    def test_csv_contains_data_rows(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["csv"])
        out = tmp_path / "summary.csv"
        reporter.generate_csv(results, out)
        lines = out.read_text().strip().split("\n")
        # Header + at least one data row per result pair
        assert len(lines) > 1

    def test_csv_empty_result_writes_row(
        self, tmp_path: Path, empty_results: list[DesignResult]
    ) -> None:
        reporter = _make_reporter(["csv"])
        out = tmp_path / "empty.csv"
        reporter.generate_csv(empty_results, out)
        content = out.read_text()
        assert "no_primers" in content


# ---------------------------------------------------------------------------
# JSON
# ---------------------------------------------------------------------------


class TestGenerateJson:
    """Tests for generate_json()."""

    def test_creates_file(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["json"])
        out = tmp_path / "results.json"
        reporter.generate_json(results, out)
        assert out.exists()

    def test_json_is_valid(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["json"])
        out = tmp_path / "results.json"
        reporter.generate_json(results, out)
        data = json.loads(out.read_text())
        assert isinstance(data, list)

    def test_json_contains_target_id(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["json"])
        out = tmp_path / "results.json"
        reporter.generate_json(results, out)
        data = json.loads(out.read_text())
        target_ids = {item.get("target", {}).get("target_id") for item in data}
        assert "ACTB" in target_ids

    def test_json_empty_result(self, tmp_path: Path, empty_results: list[DesignResult]) -> None:
        reporter = _make_reporter(["json"])
        out = tmp_path / "empty.json"
        reporter.generate_json(empty_results, out)
        data = json.loads(out.read_text())
        assert len(data) == 1


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------


class TestGenerateFasta:
    """Tests for generate_fasta()."""

    def test_creates_file(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["fasta"])
        out = tmp_path / "primers.fasta"
        reporter.generate_fasta(results, out)
        assert out.exists()

    def test_fasta_has_left_header(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["fasta"])
        out = tmp_path / "primers.fasta"
        reporter.generate_fasta(results, out)
        content = out.read_text()
        assert "_LEFT" in content

    def test_fasta_has_right_header(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["fasta"])
        out = tmp_path / "primers.fasta"
        reporter.generate_fasta(results, out)
        content = out.read_text()
        assert "_RIGHT" in content

    def test_fasta_has_probe_header_when_probe_present(self, tmp_path: Path) -> None:
        reporter = _make_reporter(["fasta"])
        results = [_make_result("ACTB", include_probe=True)]
        out = tmp_path / "primers.fasta"
        reporter.generate_fasta(results, out)
        content = out.read_text()
        assert "_PROBE" in content

    def test_fasta_starts_with_gt(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["fasta"])
        out = tmp_path / "primers.fasta"
        reporter.generate_fasta(results, out)
        lines = out.read_text().strip().split("\n")
        headers = [line for line in lines if line.startswith(">")]
        assert len(headers) > 0

    def test_empty_result_creates_empty_fasta(
        self, tmp_path: Path, empty_results: list[DesignResult]
    ) -> None:
        reporter = _make_reporter(["fasta"])
        out = tmp_path / "empty.fasta"
        reporter.generate_fasta(empty_results, out)
        # File created but no sequences
        content = out.read_text()
        assert ">" not in content


# ---------------------------------------------------------------------------
# Excel
# ---------------------------------------------------------------------------


class TestGenerateExcel:
    """Tests for generate_excel()."""

    def test_creates_file(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["xlsx"])
        out = tmp_path / "results.xlsx"
        reporter.generate_excel(results, out)
        assert out.exists()

    def test_excel_has_summary_sheet(self, tmp_path: Path, results: list[DesignResult]) -> None:
        from openpyxl import load_workbook

        reporter = _make_reporter(["xlsx"])
        out = tmp_path / "results.xlsx"
        reporter.generate_excel(results, out)
        wb = load_workbook(out)
        assert "Summary" in wb.sheetnames

    def test_excel_has_primers_sheet(self, tmp_path: Path, results: list[DesignResult]) -> None:
        from openpyxl import load_workbook

        reporter = _make_reporter(["xlsx"])
        out = tmp_path / "results.xlsx"
        reporter.generate_excel(results, out)
        wb = load_workbook(out)
        assert "Primers" in wb.sheetnames

    def test_excel_summary_has_data(self, tmp_path: Path, results: list[DesignResult]) -> None:
        from openpyxl import load_workbook

        reporter = _make_reporter(["xlsx"])
        out = tmp_path / "results.xlsx"
        reporter.generate_excel(results, out)
        wb = load_workbook(out)
        ws = wb["Summary"]
        # Row 1 = header, Row 2+ = data
        assert ws.max_row >= 2


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------


class TestGenerateHtml:
    """Tests for generate_html()."""

    def test_creates_file(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["html"])
        out = tmp_path / "report.html"
        reporter.generate_html(results, out)
        assert out.exists()

    def test_html_contains_doctype(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["html"])
        out = tmp_path / "report.html"
        reporter.generate_html(results, out)
        content = out.read_text()
        assert "<!DOCTYPE html>" in content or "<!doctype html>" in content.lower()

    def test_html_contains_target_id(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["html"])
        out = tmp_path / "report.html"
        reporter.generate_html(results, out)
        content = out.read_text()
        assert "ACTB" in content

    def test_html_no_primers_shows_message(
        self, tmp_path: Path, empty_results: list[DesignResult]
    ) -> None:
        reporter = _make_reporter(["html"])
        out = tmp_path / "empty.html"
        reporter.generate_html(empty_results, out)
        content = out.read_text()
        assert "No primer" in content or "no primer" in content.lower()


# ---------------------------------------------------------------------------
# generate() orchestrator
# ---------------------------------------------------------------------------


class TestGenerate:
    """Tests for the generate() orchestrator."""

    def test_generates_requested_formats(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["csv", "json"])
        generated = reporter.generate(results, tmp_path)
        names = {p.name for p in generated}
        assert "results_summary.csv" in names
        assert "results.json" in names

    def test_skips_unsupported_format_without_raising(
        self, tmp_path: Path, results: list[DesignResult]
    ) -> None:
        reporter = _make_reporter(["csv", "xml"])  # xml is unsupported
        # Should not raise
        generated = reporter.generate(results, tmp_path)
        # csv should still be generated
        assert any(p.suffix == ".csv" for p in generated)

    def test_creates_output_dir(self, tmp_path: Path, results: list[DesignResult]) -> None:
        out_dir = tmp_path / "nested" / "output"
        reporter = _make_reporter(["json"])
        reporter.generate(results, out_dir)
        assert out_dir.exists()

    def test_returns_list_of_paths(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["csv"])
        generated = reporter.generate(results, tmp_path)
        assert isinstance(generated, list)
        for p in generated:
            assert isinstance(p, Path)

    def test_empty_results_no_error(self, tmp_path: Path) -> None:
        reporter = _make_reporter(["csv", "json"])
        results = [_make_empty_result()]
        generated = reporter.generate(results, tmp_path)
        assert isinstance(generated, list)

    def test_generates_fasta_format(self, tmp_path: Path, results: list[DesignResult]) -> None:
        reporter = _make_reporter(["fasta"])
        generated = reporter.generate(results, tmp_path)
        assert any(p.suffix == ".fasta" for p in generated)
