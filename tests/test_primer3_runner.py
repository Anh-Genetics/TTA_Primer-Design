"""Tests for Primer3Runner module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import DesignTarget, PrimerPair, ProcessedSequence
from tta_primer_design.modules.primer3_runner import (
    Primer3Runner,
    _apply_relax_step,
    _parse_oligo,
)

CONFIG_DIR = Path(__file__).parent.parent / "config"

# Sequences long enough for Primer3 to design primers
SEQ_200 = (
    "GATTCGATCGATCGAATTCGATCGATCGAATTCGATCGATCGAATTCGATCGATCGAATTCGATCGATCGAAT"
    "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    "ATCGATCGATCGATCGATCGATCG"
)


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def runner(config: AppConfig) -> Primer3Runner:
    return Primer3Runner(config, config_dir=CONFIG_DIR)


@pytest.fixture
def processed_seq() -> ProcessedSequence:
    return ProcessedSequence(sequence=SEQ_200)


@pytest.fixture
def target() -> DesignTarget:
    return DesignTarget(target_id="TEST_GENE")


# ---------------------------------------------------------------------------
# load_primer3_params
# ---------------------------------------------------------------------------


class TestLoadPrimer3Params:
    def test_loads_qpcr_params(self, runner: Primer3Runner) -> None:
        params = runner.load_primer3_params()
        assert isinstance(params, dict)
        assert len(params) > 0

    def test_contains_primer_opt_tm(self, runner: Primer3Runner) -> None:
        params = runner.load_primer3_params()
        assert "PRIMER_OPT_TM" in params

    def test_missing_file_returns_empty(self, config: AppConfig, tmp_path: Path) -> None:
        r = Primer3Runner(config, config_dir=tmp_path)
        params = r.load_primer3_params()
        assert params == {}

    def test_pcr_mode_loads_pcr_config(self, tmp_path: Path) -> None:
        from tta_primer_design.config import PipelineConfig

        cfg = AppConfig()
        cfg.pipeline = PipelineConfig(mode="pcr")
        pcr_file = tmp_path / "primer3_pcr_params.yaml"
        pcr_file.write_text("PRIMER_OPT_TM: 55.0\n")
        r = Primer3Runner(cfg, config_dir=tmp_path)
        params = r.load_primer3_params()
        assert params.get("PRIMER_OPT_TM") == 55.0


# ---------------------------------------------------------------------------
# _build_seq_args
# ---------------------------------------------------------------------------


class TestBuildSeqArgs:
    def test_contains_template(
        self, runner: Primer3Runner, processed_seq: ProcessedSequence, target: DesignTarget
    ) -> None:
        args = runner._build_seq_args(processed_seq, target)
        assert args["SEQUENCE_TEMPLATE"] == SEQ_200

    def test_contains_id(
        self, runner: Primer3Runner, processed_seq: ProcessedSequence, target: DesignTarget
    ) -> None:
        args = runner._build_seq_args(processed_seq, target)
        assert args["SEQUENCE_ID"] == "TEST_GENE"

    def test_excluded_regions_converted(self, runner: Primer3Runner) -> None:
        ps = ProcessedSequence(sequence=SEQ_200, excluded_regions=[(10, 30), (50, 70)])
        t = DesignTarget(target_id="T1")
        args = runner._build_seq_args(ps, t)
        assert "SEQUENCE_EXCLUDED_REGION" in args
        assert args["SEQUENCE_EXCLUDED_REGION"] == [[10, 20], [50, 20]]

    def test_included_region_present(self, runner: Primer3Runner) -> None:
        ps = ProcessedSequence(sequence=SEQ_200, included_region=(5, 150))
        t = DesignTarget(target_id="T1")
        args = runner._build_seq_args(ps, t)
        assert args["SEQUENCE_INCLUDED_REGION"] == [5, 150]

    def test_no_excluded_when_empty(
        self, runner: Primer3Runner, processed_seq: ProcessedSequence, target: DesignTarget
    ) -> None:
        args = runner._build_seq_args(processed_seq, target)
        assert "SEQUENCE_EXCLUDED_REGION" not in args

    def test_target_regions_converted(self, runner: Primer3Runner) -> None:
        ps = ProcessedSequence(sequence=SEQ_200, target_regions=[(20, 80)])
        t = DesignTarget(target_id="T1")
        args = runner._build_seq_args(ps, t)
        assert "SEQUENCE_TARGET" in args
        assert args["SEQUENCE_TARGET"] == [[20, 60]]


# ---------------------------------------------------------------------------
# _parse_primer3_output
# ---------------------------------------------------------------------------


class TestParsePrimer3Output:
    def _make_result(self, n: int = 2) -> dict[str, Any]:
        result: dict[str, Any] = {"PRIMER_PAIR_NUM_RETURNED": n}
        for i in range(n):
            result[f"PRIMER_LEFT_{i}"] = [10 + i * 5, 20]
            result[f"PRIMER_LEFT_{i}_SEQUENCE"] = f"ATCGATCGATCGATCG{i:04d}"
            result[f"PRIMER_LEFT_{i}_TM"] = 60.0 + i
            result[f"PRIMER_LEFT_{i}_GC_PERCENT"] = 50.0
            result[f"PRIMER_LEFT_{i}_SELF_ANY_TH"] = 10.0
            result[f"PRIMER_LEFT_{i}_SELF_END_TH"] = 5.0
            result[f"PRIMER_LEFT_{i}_HAIRPIN_TH"] = 8.0
            result[f"PRIMER_LEFT_{i}_END_STABILITY"] = 3.0
            result[f"PRIMER_LEFT_{i}_PENALTY"] = 0.5
            result[f"PRIMER_RIGHT_{i}"] = [120 + i * 5, 20]
            result[f"PRIMER_RIGHT_{i}_SEQUENCE"] = f"GCTAGCTAGCTAGCTA{i:04d}"
            result[f"PRIMER_RIGHT_{i}_TM"] = 60.0 + i
            result[f"PRIMER_RIGHT_{i}_GC_PERCENT"] = 50.0
            result[f"PRIMER_RIGHT_{i}_SELF_ANY_TH"] = 10.0
            result[f"PRIMER_RIGHT_{i}_SELF_END_TH"] = 5.0
            result[f"PRIMER_RIGHT_{i}_HAIRPIN_TH"] = 8.0
            result[f"PRIMER_RIGHT_{i}_END_STABILITY"] = 3.0
            result[f"PRIMER_RIGHT_{i}_PENALTY"] = 0.5
            result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"] = 100 + i * 10
            result[f"PRIMER_PAIR_{i}_PENALTY"] = 1.0
        return result

    def test_returns_correct_count(self, runner: Primer3Runner) -> None:
        pairs = runner._parse_primer3_output(self._make_result(2), "GENE1")
        assert len(pairs) == 2

    def test_pair_ids(self, runner: Primer3Runner) -> None:
        pairs = runner._parse_primer3_output(self._make_result(1), "GENE1")
        assert pairs[0].pair_id == "GENE1_pair_0"

    def test_left_primer_sequence(self, runner: Primer3Runner) -> None:
        pairs = runner._parse_primer3_output(self._make_result(1), "GENE1")
        assert "ATCG" in pairs[0].left_primer.sequence

    def test_right_primer_tm(self, runner: Primer3Runner) -> None:
        pairs = runner._parse_primer3_output(self._make_result(1), "GENE1")
        assert pairs[0].right_primer.tm == pytest.approx(60.0)

    def test_amplicon_size_set(self, runner: Primer3Runner) -> None:
        pairs = runner._parse_primer3_output(self._make_result(1), "GENE1")
        assert pairs[0].amplicon_size == 100

    def test_zero_returned_gives_empty(self, runner: Primer3Runner) -> None:
        pairs = runner._parse_primer3_output({"PRIMER_PAIR_NUM_RETURNED": 0}, "GENE1")
        assert pairs == []

    def test_probe_parsed_when_present(self, runner: Primer3Runner) -> None:
        result = self._make_result(1)
        result["PRIMER_INTERNAL_0"] = [60, 22]
        result["PRIMER_INTERNAL_0_SEQUENCE"] = "GCATGCATGCATGCATGCATGC"
        result["PRIMER_INTERNAL_0_TM"] = 68.0
        result["PRIMER_INTERNAL_0_GC_PERCENT"] = 50.0
        result["PRIMER_INTERNAL_0_SELF_ANY_TH"] = 5.0
        result["PRIMER_INTERNAL_0_SELF_END_TH"] = 2.0
        result["PRIMER_INTERNAL_0_HAIRPIN_TH"] = 3.0
        result["PRIMER_INTERNAL_0_END_STABILITY"] = 1.0
        result["PRIMER_INTERNAL_0_PENALTY"] = 0.2
        pairs = runner._parse_primer3_output(result, "GENE1")
        assert pairs[0].probe is not None
        assert pairs[0].probe.tm == pytest.approx(68.0)


# ---------------------------------------------------------------------------
# run() — integration with primer3-py (needs suitable sequence)
# ---------------------------------------------------------------------------


class TestRun:
    def _relaxed_global_args(self) -> dict[str, Any]:
        """Minimal global_args that work for SEQ_200."""
        return {
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 15,
            "PRIMER_MAX_SIZE": 25,
            "PRIMER_NUM_RETURN": 3,
            "PRIMER_PRODUCT_SIZE_RANGE": "50-170",
            "PRIMER_OPT_TM": 55.0,
            "PRIMER_MIN_TM": 45.0,
            "PRIMER_MAX_TM": 65.0,
            "PRIMER_MIN_GC": 30.0,
            "PRIMER_MAX_GC": 70.0,
        }

    def test_run_returns_list(
        self, runner: Primer3Runner, processed_seq: ProcessedSequence, target: DesignTarget
    ) -> None:
        # Override global args via custom_params for a relaxed design
        t = DesignTarget(target_id="TEST_GENE", custom_params=self._relaxed_global_args())
        result = runner.run(processed_seq, t)
        assert isinstance(result, list)

    def test_run_finds_primers_with_relaxed_constraints(
        self, runner: Primer3Runner, processed_seq: ProcessedSequence
    ) -> None:
        t = DesignTarget(target_id="TEST_GENE", custom_params=self._relaxed_global_args())
        pairs = runner.run(processed_seq, t)
        assert len(pairs) > 0

    def test_run_primer_pair_type(
        self, runner: Primer3Runner, processed_seq: ProcessedSequence
    ) -> None:
        t = DesignTarget(target_id="TEST_GENE", custom_params=self._relaxed_global_args())
        pairs = runner.run(processed_seq, t)
        for pair in pairs:
            assert isinstance(pair, PrimerPair)

    def test_run_empty_result_when_impossible(self, runner: Primer3Runner) -> None:
        # Very short sequence → Primer3 cannot design primers
        ps = ProcessedSequence(sequence="ATCGATCGATCG")
        t = DesignTarget(target_id="TINY")
        pairs = runner.run(ps, t)
        assert isinstance(pairs, list)

    def test_run_custom_params_merged(
        self, runner: Primer3Runner, processed_seq: ProcessedSequence
    ) -> None:
        custom = dict(self._relaxed_global_args())
        custom["PRIMER_NUM_RETURN"] = 1
        t = DesignTarget(target_id="T1", custom_params=custom)
        pairs = runner.run(processed_seq, t)
        assert len(pairs) <= 1


# ---------------------------------------------------------------------------
# _parse_oligo helper
# ---------------------------------------------------------------------------


class TestParseOligo:
    def test_returns_oligo(self) -> None:
        result: dict[str, Any] = {
            "PRIMER_LEFT_0": [5, 20],
            "PRIMER_LEFT_0_SEQUENCE": "ATCGATCGATCGATCGATCG",
            "PRIMER_LEFT_0_TM": 60.0,
            "PRIMER_LEFT_0_GC_PERCENT": 50.0,
            "PRIMER_LEFT_0_SELF_ANY_TH": 10.0,
            "PRIMER_LEFT_0_SELF_END_TH": 5.0,
            "PRIMER_LEFT_0_HAIRPIN_TH": 8.0,
            "PRIMER_LEFT_0_END_STABILITY": 3.0,
            "PRIMER_LEFT_0_PENALTY": 0.5,
        }
        oligo = _parse_oligo(result, "LEFT", 0)
        assert oligo is not None
        assert oligo.sequence == "ATCGATCGATCGATCGATCG"
        assert oligo.tm == pytest.approx(60.0)
        assert oligo.start == 5

    def test_returns_none_when_no_sequence(self) -> None:
        assert _parse_oligo({}, "LEFT", 0) is None


# ---------------------------------------------------------------------------
# _apply_relax_step helper
# ---------------------------------------------------------------------------


class TestApplyRelaxStep:
    def test_min_tm_decreased(self) -> None:
        params = {"PRIMER_MIN_TM": 58.0}
        relaxed = _apply_relax_step(params, {"PRIMER_MIN_TM": 2.0})
        assert relaxed["PRIMER_MIN_TM"] == pytest.approx(56.0)

    def test_max_tm_increased(self) -> None:
        params = {"PRIMER_MAX_TM": 62.0}
        relaxed = _apply_relax_step(params, {"PRIMER_MAX_TM": 2.0})
        assert relaxed["PRIMER_MAX_TM"] == pytest.approx(64.0)

    def test_missing_key_not_added(self) -> None:
        params: dict[str, Any] = {}
        relaxed = _apply_relax_step(params, {"PRIMER_MIN_TM": 2.0})
        assert "PRIMER_MIN_TM" not in relaxed

    def test_original_not_mutated(self) -> None:
        params = {"PRIMER_MIN_TM": 58.0}
        _apply_relax_step(params, {"PRIMER_MIN_TM": 2.0})
        assert params["PRIMER_MIN_TM"] == pytest.approx(58.0)
