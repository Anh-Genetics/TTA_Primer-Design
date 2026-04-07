"""Tests cho FilterRanker module."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from tta_primer_design.config import AppConfig, FiltersConfig
from tta_primer_design.models import Oligo, PrimerPair
from tta_primer_design.modules.filter_ranker import FilterRanker


def _make_oligo(seq: str = "ATCGATCGATCGATCGATCG", tm: float = 60.0) -> Oligo:
    return Oligo(sequence=seq, tm=tm)


def _make_pair(
    pair_id: str = "pair_01",
    left_tm: float = 60.0,
    right_tm: float = 60.0,
    snp_flags: list[str] | None = None,
    off_targets: int = 0,
    specificity_score: float = 100.0,
    amplicon_size: int = 120,
    pair_penalty: float = 0.0,
) -> PrimerPair:
    left = _make_oligo(tm=left_tm)
    right = _make_oligo(tm=right_tm)
    pair = PrimerPair(
        pair_id=pair_id,
        left_primer=left,
        right_primer=right,
        amplicon_size=amplicon_size,
        pair_penalty=pair_penalty,
        snp_flags=snp_flags or [],
    )
    if off_targets >= 0:
        mock_spec = MagicMock()
        mock_spec.off_target_amplicons = [MagicMock()] * off_targets
        mock_spec.specificity_score = specificity_score
        pair.specificity_result = mock_spec
    return pair


@pytest.fixture()
def config_zero_off_targets() -> AppConfig:
    cfg = AppConfig()
    cfg.filters = FiltersConfig(max_off_targets=0)
    return cfg


@pytest.fixture()
def ranker(config_zero_off_targets: AppConfig) -> FilterRanker:
    return FilterRanker(config_zero_off_targets)


class TestApplyHardFilters:
    """Tests for apply_hard_filters()."""

    def test_removes_pair_with_off_target_amplicons(self, ranker: FilterRanker) -> None:
        good = _make_pair("good", off_targets=0)
        bad = _make_pair("bad", off_targets=1)
        result = ranker.apply_hard_filters([good, bad])
        ids = [p.pair_id for p in result]
        assert "good" in ids
        assert "bad" not in ids

    def test_keeps_pair_with_no_off_targets(self, ranker: FilterRanker) -> None:
        pair = _make_pair("clean", off_targets=0)
        result = ranker.apply_hard_filters([pair])
        assert len(result) == 1
        assert result[0].pair_id == "clean"

    def test_removes_pair_with_fail_snp_flag(self, ranker: FilterRanker) -> None:
        pair = _make_pair("snp_fail", snp_flags=["FAIL_SNP_3PRIME"])
        result = ranker.apply_hard_filters([pair])
        assert len(result) == 0

    def test_keeps_pair_with_no_snp_flags(self, ranker: FilterRanker) -> None:
        pair = _make_pair("clean_snp")
        result = ranker.apply_hard_filters([pair])
        assert len(result) == 1

    def test_keeps_pair_with_warning_snp_flag(self, ranker: FilterRanker) -> None:
        pair = _make_pair("warn", snp_flags=["SNP_WARNING"])
        result = ranker.apply_hard_filters([pair])
        assert len(result) == 1

    def test_empty_input(self, ranker: FilterRanker) -> None:
        assert ranker.apply_hard_filters([]) == []

    def test_pair_without_specificity_result_kept(self, ranker: FilterRanker) -> None:
        pair = PrimerPair(
            pair_id="no_spec",
            left_primer=_make_oligo(),
            right_primer=_make_oligo(),
        )
        pair.specificity_result = None
        result = ranker.apply_hard_filters([pair])
        assert len(result) == 1


class TestApplySoftFilters:
    """Tests for apply_soft_filters()."""

    def test_adds_snp_warning_for_snp_prefix_flag(self, ranker: FilterRanker) -> None:
        pair = _make_pair("warn", snp_flags=["SNP_MIDDLE"])
        ranker.apply_soft_filters([pair])
        assert "SNP_WARNING" in pair.snp_flags

    def test_no_duplicate_snp_warning(self, ranker: FilterRanker) -> None:
        pair = _make_pair("warn2", snp_flags=["SNP_LOW_MAF", "SNP_WARNING"])
        ranker.apply_soft_filters([pair])
        assert pair.snp_flags.count("SNP_WARNING") == 1

    def test_no_snp_warning_without_snp_prefix(self, ranker: FilterRanker) -> None:
        pair = _make_pair("clean")
        ranker.apply_soft_filters([pair])
        assert "SNP_WARNING" not in pair.snp_flags

    def test_returns_all_pairs(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(3)]
        result = ranker.apply_soft_filters(pairs)
        assert len(result) == 3

    def test_hairpin_warning_added_for_high_hairpin(self, ranker: FilterRanker) -> None:
        pair = _make_pair("hp")
        pair.left_primer.hairpin_th = 25.0
        ranker.apply_soft_filters([pair])
        hairpin_flags = [f for f in pair.snp_flags if "HAIRPIN_WARNING" in f]
        assert len(hairpin_flags) == 1


class TestCalculateScore:
    """Tests for calculate_score()."""

    def test_returns_float(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        score = ranker.calculate_score(pair)
        assert isinstance(score, float)

    def test_score_in_0_to_100_range(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        score = ranker.calculate_score(pair)
        assert 0.0 <= score <= 100.0

    def test_good_tm_scores_higher_than_bad_tm(self, ranker: FilterRanker) -> None:
        good = _make_pair("good_tm", left_tm=60.0, right_tm=60.0)
        bad = _make_pair("bad_tm", left_tm=40.0, right_tm=40.0)
        assert ranker.calculate_score(good) > ranker.calculate_score(bad)

    def test_no_snp_scores_higher_than_fail_snp(self, ranker: FilterRanker) -> None:
        clean = _make_pair("clean_snp")
        fail = _make_pair("fail_snp", snp_flags=["FAIL_SOMETHING"])
        assert ranker.calculate_score(clean) > ranker.calculate_score(fail)

    def test_no_specificity_result_uses_default(self, ranker: FilterRanker) -> None:
        pair = PrimerPair(
            pair_id="no_spec",
            left_primer=_make_oligo(tm=60.0),
            right_primer=_make_oligo(tm=60.0),
        )
        pair.specificity_result = None
        score = ranker.calculate_score(pair)
        assert 0.0 <= score <= 100.0

    def test_high_specificity_scores_higher(self, ranker: FilterRanker) -> None:
        high = _make_pair("high_spec", specificity_score=100.0)
        low = _make_pair("low_spec", specificity_score=0.0)
        assert ranker.calculate_score(high) > ranker.calculate_score(low)

    def test_optimal_amplicon_size_scores_better(self, ranker: FilterRanker) -> None:
        optimal = _make_pair("opt", amplicon_size=125)
        extreme = _make_pair("ext", amplicon_size=500)
        assert ranker.calculate_score(optimal) > ranker.calculate_score(extreme)


class TestRankPairs:
    """Tests for rank_pairs()."""

    def test_sorted_by_score_descending(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(3)]
        pairs[0].score = 50.0
        pairs[1].score = 80.0
        pairs[2].score = 60.0
        ranked = ranker.rank_pairs(pairs)
        scores = [p.score for p in ranked]
        assert scores == sorted(scores, reverse=True)

    def test_empty_list(self, ranker: FilterRanker) -> None:
        assert ranker.rank_pairs([]) == []

    def test_returns_all_pairs(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(5)]
        for i, p in enumerate(pairs):
            p.score = float(i * 10)
        assert len(ranker.rank_pairs(pairs)) == 5


class TestGetTopN:
    """Tests for get_top_n()."""

    def test_limits_to_n(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(10)]
        result = ranker.get_top_n(pairs, n=3)
        assert len(result) == 3

    def test_returns_first_n(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(10)]
        result = ranker.get_top_n(pairs, n=2)
        assert result == pairs[:2]

    def test_uses_config_top_n_when_n_not_given(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(10)]
        # Default config has top_n_pairs=5
        result = ranker.get_top_n(pairs)
        assert len(result) == 5

    def test_fewer_pairs_than_n(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}") for i in range(2)]
        result = ranker.get_top_n(pairs, n=10)
        assert len(result) == 2


class TestProcess:
    """Tests for process()."""

    def test_returns_list(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        result = ranker.process([pair])
        assert isinstance(result, list)

    def test_empty_input(self, ranker: FilterRanker) -> None:
        assert ranker.process([]) == []

    def test_removes_off_target_pairs(self, ranker: FilterRanker) -> None:
        good = _make_pair("good", off_targets=0)
        bad = _make_pair("bad", off_targets=2)
        result = ranker.process([good, bad])
        assert all(p.pair_id == "good" for p in result)

    def test_sets_score_on_pairs(self, ranker: FilterRanker) -> None:
        pair = _make_pair()
        result = ranker.process([pair])
        for p in result:
            assert p.score > 0.0

    def test_result_sorted_by_score(self, ranker: FilterRanker) -> None:
        pairs = [_make_pair(f"p{i}", left_tm=50.0 + i * 2) for i in range(5)]
        result = ranker.process(pairs)
        scores = [p.score for p in result]
        assert scores == sorted(scores, reverse=True)
