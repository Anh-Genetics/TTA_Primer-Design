"""Tests cho FilterRanker — filter + rank primer pairs."""

from __future__ import annotations

import pytest

from tta_primer_design.config import AppConfig, FiltersConfig, PipelineConfig, load_config
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.filter_ranker import FilterRanker

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def ranker(config: AppConfig) -> FilterRanker:
    return FilterRanker(config)


def _make_oligo(
    seq: str = "ATCGATCGATCGATCGATCG",
    tm: float = 60.0,
    gc: float = 50.0,
    self_any_th: float = 0.0,
    hairpin_th: float = 0.0,
    penalty: float = 0.0,
) -> Oligo:
    return Oligo(
        sequence=seq,
        tm=tm,
        gc_percent=gc,
        self_any_th=self_any_th,
        hairpin_th=hairpin_th,
        penalty=penalty,
    )


def _make_pair(
    pair_id: str = "pair_0",
    left_tm: float = 60.0,
    right_tm: float = 60.0,
    amplicon_size: int = 120,
    snp_flags: list[str] | None = None,
    specificity_result=None,
    pair_penalty: float = 1.0,
    probe: Oligo | None = None,
) -> PrimerPair:
    return PrimerPair(
        pair_id=pair_id,
        left_primer=_make_oligo(tm=left_tm),
        right_primer=_make_oligo(seq="GCTAGCTAGCTAGCTAGCTA", tm=right_tm),
        amplicon_size=amplicon_size,
        snp_flags=snp_flags or [],
        specificity_result=specificity_result,
        pair_penalty=pair_penalty,
        probe=probe,
    )


# ---------------------------------------------------------------------------
# Hard-filter helpers
# ---------------------------------------------------------------------------


class _MockSpecificity:
    def __init__(self, off_targets: int = 0, specificity_score: float = 100.0) -> None:

        self.off_target_amplicons = [object()] * off_targets
        self.specificity_score = specificity_score
        self.is_specific = off_targets == 0


# ---------------------------------------------------------------------------
# Tests: apply_hard_filters
# ---------------------------------------------------------------------------


class TestApplyHardFilters:
    def test_clean_pair_passes(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair()]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 1

    def test_off_target_removed(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(specificity_result=_MockSpecificity(off_targets=1))]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 0

    def test_zero_off_targets_passes(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(specificity_result=_MockSpecificity(off_targets=0))]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 1

    def test_critical_snp_fail_removed(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(snp_flags=["FAIL:snp_3prime"])]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 0

    def test_warning_snp_not_removed(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(snp_flags=["WARNING:snp_mid"])]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 1

    def test_tm_too_low_removed(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(left_tm=30.0)]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 0

    def test_tm_too_high_removed(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(right_tm=85.0)]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 0

    def test_tm_zero_not_removed(self, ranker: FilterRanker) -> None:
        """Tm=0 berarti belum dihitung — tidak di-filter."""
        pairs = [_make_pair(left_tm=0.0, right_tm=0.0)]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 1

    def test_multiple_pairs_mixed(self, ranker: FilterRanker) -> None:
        pairs = [
            _make_pair("good"),
            _make_pair("bad_off_target", specificity_result=_MockSpecificity(off_targets=2)),
            _make_pair("bad_snp", snp_flags=["FAIL:probe_snp"]),
            _make_pair("good2"),
        ]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 2
        ids = [p.pair_id for p in result]
        assert "good" in ids
        assert "good2" in ids

    def test_empty_input(self, ranker: FilterRanker) -> None:
        assert ranker.apply_hard_filters([]) == []

    def test_probe_tm_out_of_range_removed(self, ranker: FilterRanker) -> None:
        probe = _make_oligo(tm=90.0)
        pairs = [_make_pair(probe=probe)]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 0

    def test_config_max_off_targets_1(self) -> None:
        cfg = AppConfig()
        cfg.filters = FiltersConfig(max_off_targets=1)
        ranker = FilterRanker(cfg)
        pairs = [_make_pair(specificity_result=_MockSpecificity(off_targets=1))]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 1  # 1 == max_off_targets → còn pass

    def test_config_max_off_targets_exceeded(self) -> None:
        cfg = AppConfig()
        cfg.filters = FiltersConfig(max_off_targets=1)
        ranker = FilterRanker(cfg)
        pairs = [_make_pair(specificity_result=_MockSpecificity(off_targets=2))]
        result = ranker.apply_hard_filters(pairs)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# Tests: apply_soft_filters
# ---------------------------------------------------------------------------


class TestApplySoftFilters:
    def test_no_flags_unchanged(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair()]
        result = ranker.apply_soft_filters(pairs)
        assert result[0].snp_flags == []

    def test_dg_warning_added(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        pair.left_primer.self_any_th = -15.0  # below threshold
        result = ranker.apply_soft_filters([pair])
        assert any("dg_high" in f for f in result[0].snp_flags)

    def test_hairpin_dg_warning_added(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        pair.right_primer.hairpin_th = -12.0
        result = ranker.apply_soft_filters([pair])
        assert any("dg_high" in f for f in result[0].snp_flags)

    def test_snp_warning_flag_added(self, ranker: FilterRanker) -> None:
        pair = _make_pair(snp_flags=["WARNING:mid_position"])
        result = ranker.apply_soft_filters([pair])
        assert any("snp_warning" in f for f in result[0].snp_flags)

    def test_soft_filter_does_not_remove_pairs(self, ranker: FilterRanker) -> None:
        pairs = [
            _make_pair("p1", snp_flags=["WARNING:mid"]),
            _make_pair("p2"),
        ]
        result = ranker.apply_soft_filters(pairs)
        assert len(result) == 2

    def test_no_duplicate_dg_flag(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        pair.left_primer.self_any_th = -15.0
        pairs = [pair]
        ranker.apply_soft_filters(pairs)
        ranker.apply_soft_filters(pairs)  # apply twice
        dg_flags = [f for f in pair.snp_flags if "dg_high" in f]
        assert len(dg_flags) == 1

    def test_empty_input(self, ranker: FilterRanker) -> None:
        assert ranker.apply_soft_filters([]) == []


# ---------------------------------------------------------------------------
# Tests: calculate_score
# ---------------------------------------------------------------------------


class TestCalculateScore:
    def test_score_in_range(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        score = ranker.calculate_score(pair)
        assert 0.0 <= score <= 100.0

    def test_perfect_pair_high_score(self, ranker: FilterRanker) -> None:
        pair = _make_pair(
            specificity_result=_MockSpecificity(off_targets=0, specificity_score=100.0),
            amplicon_size=120,
            pair_penalty=0.0,
        )
        score = ranker.calculate_score(pair)
        assert score > 60.0

    def test_bad_pair_lower_score(self, ranker: FilterRanker) -> None:
        good = _make_pair(
            specificity_result=_MockSpecificity(specificity_score=100.0),
            pair_penalty=0.0,
        )
        bad = _make_pair(
            snp_flags=["WARNING:a", "WARNING:b"],
            pair_penalty=15.0,
        )
        good_score = ranker.calculate_score(good)
        bad_score = ranker.calculate_score(bad)
        assert good_score > bad_score

    def test_no_specificity_result_uses_default(self, ranker: FilterRanker) -> None:
        pair = _make_pair(specificity_result=None)
        score = ranker.calculate_score(pair)
        assert 0.0 <= score <= 100.0

    def test_score_with_probe(self, ranker: FilterRanker) -> None:
        probe = _make_oligo(tm=67.0)
        pair = _make_pair(probe=probe)
        score = ranker.calculate_score(pair)
        assert 0.0 <= score <= 100.0

    def test_tm_in_optimal_range_gets_bonus(self, ranker: FilterRanker) -> None:
        good = _make_pair(left_tm=60.0, right_tm=60.0, pair_penalty=0.0)
        bad = _make_pair(left_tm=50.0, right_tm=50.0, pair_penalty=0.0)
        good.score = 0.0
        bad.score = 0.0
        assert ranker.calculate_score(good) >= ranker.calculate_score(bad)

    @pytest.mark.parametrize(
        "size,mode,expect_high",
        [
            (120, "qpcr", True),
            (50, "pcr", False),  # below PCR min (100)
            (500, "pcr", True),
            (5000, "pcr", False),
        ],
    )
    def test_amplicon_size_score(self, size: int, mode: str, expect_high: bool) -> None:
        cfg = AppConfig()
        cfg.pipeline = PipelineConfig(mode=mode)
        ranker = FilterRanker(cfg)
        pair = _make_pair(amplicon_size=size)
        amplicon_component = ranker._amplicon_size_score(pair)
        if expect_high:
            assert amplicon_component >= 0.7
        else:
            assert amplicon_component < 1.0


# ---------------------------------------------------------------------------
# Tests: rank_pairs & get_top_n
# ---------------------------------------------------------------------------


class TestRankAndTopN:
    def test_rank_by_score_desc(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair("p1"), _make_pair("p2"), _make_pair("p3")]
        pairs[0].score = 50.0
        pairs[1].score = 90.0
        pairs[2].score = 70.0
        ranked = ranker.rank_pairs(pairs)
        assert ranked[0].pair_id == "p2"
        assert ranked[1].pair_id == "p3"
        assert ranked[2].pair_id == "p1"

    def test_get_top_n_default(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(10)]
        for i, p in enumerate(pairs):
            p.score = float(i)
        ranked = ranker.rank_pairs(pairs)
        top = ranker.get_top_n(ranked)
        assert len(top) == ranker.config.pipeline.top_n_pairs

    def test_get_top_n_explicit(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(8)]
        result = ranker.get_top_n(pairs, n=3)
        assert len(result) == 3

    def test_get_top_n_fewer_than_n(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair("p0")]
        result = ranker.get_top_n(pairs, n=5)
        assert len(result) == 1

    def test_rank_empty(self, ranker: FilterRanker) -> None:
        assert ranker.rank_pairs([]) == []


# ---------------------------------------------------------------------------
# Tests: process (end-to-end)
# ---------------------------------------------------------------------------


class TestProcess:
    def test_process_returns_list(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair("p0"), _make_pair("p1")]
        result = ranker.process(pairs)
        assert isinstance(result, list)

    def test_process_scores_assigned(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair("p0"), _make_pair("p1")]
        result = ranker.process(pairs)
        for p in result:
            assert p.score > 0.0

    def test_process_respects_top_n(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(20)]
        result = ranker.process(pairs)
        assert len(result) <= ranker.config.pipeline.top_n_pairs

    def test_process_removes_hard_filter_pairs(self, ranker: FilterRanker) -> None:
        good = _make_pair("good")
        bad = _make_pair("bad", snp_flags=["FAIL:critical"])
        result = ranker.process([good, bad])
        ids = [p.pair_id for p in result]
        assert "bad" not in ids

    def test_process_empty_input(self, ranker: FilterRanker) -> None:
        result = ranker.process([])
        assert result == []

    def test_process_all_filtered_out(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(snp_flags=["FAIL:all"]) for _ in range(5)]
        result = ranker.process(pairs)
        assert result == []

    def test_process_sorted_desc(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(10)]
        result = ranker.process(pairs)
        scores = [p.score for p in result]
        assert scores == sorted(scores, reverse=True)
