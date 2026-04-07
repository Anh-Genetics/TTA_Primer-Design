"""Tests cho BlastVisualizer module."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from tta_primer_design.modules.blast_specificity import (
    BlastHit,
    OffTargetAmplicon,
    SpecificityResult,
)
from tta_primer_design.modules.blast_visualizer import BlastVisualizer

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_hit(subject_id: str, mismatches: int, evalue: float = 0.01) -> BlastHit:
    return BlastHit(
        subject_id=subject_id,
        mismatches=mismatches,
        evalue=evalue,
        identity=100.0 - mismatches * 5.0,
    )


def _make_spec_result(
    pair_id: str = "ACTB_pair1",
    left_mismatches: list[int] | None = None,
    right_mismatches: list[int] | None = None,
    off_target_sizes: list[int] | None = None,
) -> SpecificityResult:
    if left_mismatches is None:
        left_mismatches = [0, 0, 1, 2, 3]
    if right_mismatches is None:
        right_mismatches = [0, 1, 2]
    if off_target_sizes is None:
        off_target_sizes = [250]

    left_hits = [_make_hit(f"NM_{i:04d}", mm) for i, mm in enumerate(left_mismatches)]
    right_hits = [_make_hit(f"NM_{i:04d}", mm) for i, mm in enumerate(right_mismatches)]

    off_targets = [
        OffTargetAmplicon(
            subject_id=f"off_{i}",
            amplicon_size=size,
            left_hit=left_hits[0],
            right_hit=right_hits[0],
        )
        for i, size in enumerate(off_target_sizes)
    ]

    return SpecificityResult(
        primer_pair_id=pair_id,
        blast_hits_left=left_hits,
        blast_hits_right=right_hits,
        off_target_amplicons=off_targets,
        specificity_score=100.0 - len(off_targets) * 20.0,
        is_specific=len(off_targets) == 0,
    )


@pytest.fixture()
def spec_result() -> SpecificityResult:
    return _make_spec_result()


@pytest.fixture()
def empty_spec_result() -> SpecificityResult:
    return SpecificityResult(
        primer_pair_id="empty_pair",
        blast_hits_left=[],
        blast_hits_right=[],
        off_target_amplicons=[],
        specificity_score=100.0,
        is_specific=True,
    )


# ---------------------------------------------------------------------------
# print_ascii_chart
# ---------------------------------------------------------------------------


class TestPrintAsciiChart:
    """Tests for BlastVisualizer.print_ascii_chart()."""

    def test_runs_without_error(
        self, spec_result: SpecificityResult, capsys: pytest.CaptureFixture
    ) -> None:
        viz = BlastVisualizer(spec_result, "ACTB_pair1")
        viz.print_ascii_chart()  # should not raise

    def test_outputs_to_stdout(
        self, spec_result: SpecificityResult, capsys: pytest.CaptureFixture
    ) -> None:
        viz = BlastVisualizer(spec_result, "ACTB_pair1")
        viz.print_ascii_chart()
        captured = capsys.readouterr()
        assert len(captured.out) > 0

    def test_shows_pair_name_in_output(
        self, spec_result: SpecificityResult, capsys: pytest.CaptureFixture
    ) -> None:
        viz = BlastVisualizer(spec_result, "ACTB_pair1")
        viz.print_ascii_chart()
        captured = capsys.readouterr()
        assert "ACTB_pair1" in captured.out

    def test_shows_no_hits_for_empty_result(
        self, empty_spec_result: SpecificityResult, capsys: pytest.CaptureFixture
    ) -> None:
        viz = BlastVisualizer(empty_spec_result, "empty_pair")
        viz.print_ascii_chart()
        captured = capsys.readouterr()
        assert "(no hits)" in captured.out

    def test_shows_off_target_amplicons_section(
        self, spec_result: SpecificityResult, capsys: pytest.CaptureFixture
    ) -> None:
        viz = BlastVisualizer(spec_result)
        viz.print_ascii_chart()
        captured = capsys.readouterr()
        # Should mention off-target amplicons
        assert "amplicon" in captured.out.lower() or "off-target" in captured.out.lower()

    def test_shows_no_off_target_when_empty(
        self, empty_spec_result: SpecificityResult, capsys: pytest.CaptureFixture
    ) -> None:
        viz = BlastVisualizer(empty_spec_result)
        viz.print_ascii_chart()
        captured = capsys.readouterr()
        # Should show the "no off-target" message
        assert "off-target" in captured.out.lower() or "amplicon" in captured.out.lower()

    def test_pair_name_defaults_to_result_id(self, spec_result: SpecificityResult) -> None:
        viz = BlastVisualizer(spec_result)
        assert viz.pair_name == spec_result.primer_pair_id

    def test_custom_pair_name_used(self, spec_result: SpecificityResult) -> None:
        viz = BlastVisualizer(spec_result, pair_name="CUSTOM_NAME")
        assert viz.pair_name == "CUSTOM_NAME"

    def test_shows_mismatch_distribution_header(
        self, spec_result: SpecificityResult, capsys: pytest.CaptureFixture
    ) -> None:
        viz = BlastVisualizer(spec_result)
        viz.print_ascii_chart()
        captured = capsys.readouterr()
        assert "Mismatch" in captured.out or "mismatch" in captured.out.lower()


# ---------------------------------------------------------------------------
# save_plot
# ---------------------------------------------------------------------------


class TestSavePlot:
    """Tests for BlastVisualizer.save_plot()."""

    def test_saves_png_file(self, spec_result: SpecificityResult, tmp_path: Path) -> None:
        out = tmp_path / "blast_plot.png"
        viz = BlastVisualizer(spec_result)
        result_path = viz.save_plot(out)
        assert out.exists()
        assert result_path == out

    def test_creates_parent_dirs(self, spec_result: SpecificityResult, tmp_path: Path) -> None:
        out = tmp_path / "nested" / "subdir" / "plot.png"
        viz = BlastVisualizer(spec_result)
        viz.save_plot(out)
        assert out.exists()

    def test_raises_import_error_when_matplotlib_missing(
        self, spec_result: SpecificityResult, tmp_path: Path
    ) -> None:
        out = tmp_path / "plot.png"
        viz = BlastVisualizer(spec_result)
        # Temporarily hide matplotlib
        with patch.dict(sys.modules, {"matplotlib": None, "matplotlib.pyplot": None}):
            with pytest.raises(ImportError, match="matplotlib"):
                viz.save_plot(out)

    def test_save_plot_with_no_off_targets(
        self, empty_spec_result: SpecificityResult, tmp_path: Path
    ) -> None:
        out = tmp_path / "empty_plot.png"
        viz = BlastVisualizer(empty_spec_result)
        viz.save_plot(out)
        assert out.exists()

    def test_save_plot_with_off_targets(self, tmp_path: Path) -> None:
        # Multiple off-target amplicons → triggers subplot 3
        result = _make_spec_result(off_target_sizes=[250, 300, 450])
        out = tmp_path / "with_off_targets.png"
        viz = BlastVisualizer(result)
        viz.save_plot(out)
        assert out.exists()

    def test_png_file_is_non_empty(self, spec_result: SpecificityResult, tmp_path: Path) -> None:
        out = tmp_path / "plot.png"
        viz = BlastVisualizer(spec_result)
        viz.save_plot(out)
        assert out.stat().st_size > 0
