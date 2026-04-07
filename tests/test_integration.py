"""Integration tests — kiểm tra luồng filter → rank → report end-to-end."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from tta_primer_design.config import OutputConfig, PipelineConfig, load_config
from tta_primer_design.models import DesignResult, DesignTarget, Oligo, PrimerPair
from tta_primer_design.modules.filter_ranker import FilterRanker
from tta_primer_design.modules.report_generator import ReportGenerator

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _make_target(tid: str = "ACTB") -> DesignTarget:
    return DesignTarget(target_id=tid, accession="NM_001101", design_mode="qpcr")


def _make_oligo(
    seq: str = "ATCGATCGATCGATCGATCG",
    tm: float = 60.0,
    gc: float = 50.0,
    self_any_th: float = 0.0,
    hairpin_th: float = 0.0,
) -> Oligo:
    return Oligo(sequence=seq, tm=tm, gc_percent=gc, self_any_th=self_any_th, hairpin_th=hairpin_th)


def _make_pair(
    pair_id: str = "pair_0",
    left_tm: float = 60.0,
    right_tm: float = 60.0,
    amplicon_size: int = 120,
    snp_flags: list[str] | None = None,
    pair_penalty: float = 1.0,
    probe: Oligo | None = None,
) -> PrimerPair:
    return PrimerPair(
        pair_id=pair_id,
        left_primer=_make_oligo(tm=left_tm),
        right_primer=_make_oligo(seq="GCTAGCTAGCTAGCTAGCTA", tm=right_tm),
        amplicon_size=amplicon_size,
        snp_flags=snp_flags or [],
        pair_penalty=pair_penalty,
        probe=probe,
    )


# ---------------------------------------------------------------------------
# Tests: FilterRanker → ReportGenerator pipeline
# ---------------------------------------------------------------------------


class TestFilterToReport:
    """Kiểm tra toàn bộ luồng: filter + rank → tạo báo cáo."""

    def test_full_pipeline_csv(self, tmp_path: Path) -> None:
        """FilterRanker xử lý pairs → ReportGenerator tạo CSV."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["csv"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [
            _make_pair("p0", left_tm=60.0, right_tm=60.0, amplicon_size=120),
            _make_pair("p1", left_tm=62.0, right_tm=61.0, amplicon_size=110),
            _make_pair("bad", snp_flags=["FAIL:critical"]),
        ]
        ranked = ranker.process(pairs)

        target = _make_target()
        results = [DesignResult(target=target, primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)

        assert len(paths) == 1
        csv_path = paths[0]
        assert csv_path.exists()
        with csv_path.open() as fh:
            rows = list(csv.DictReader(fh))
        # "bad" pair filtered out
        pair_ids = [r["pair_id"] for r in rows]
        assert "bad" not in pair_ids

    def test_full_pipeline_json(self, tmp_path: Path) -> None:
        """FilterRanker → ReportGenerator tạo JSON hợp lệ."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["json"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [_make_pair("p0"), _make_pair("p1")]
        ranked = ranker.process(pairs)

        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)

        data = json.loads(paths[0].read_text())
        assert len(data) == 1
        assert data[0]["status"] == "success"

    def test_scores_assigned_before_report(self, tmp_path: Path) -> None:
        """Score phải được gán bởi FilterRanker trước khi report."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["csv"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [_make_pair("p0"), _make_pair("p1")]
        ranked = ranker.process(pairs)
        for p in ranked:
            assert p.score > 0.0

        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)
        with paths[0].open() as fh:
            rows = list(csv.DictReader(fh))
        for row in rows:
            assert float(row["score"]) > 0.0

    def test_ranked_order_preserved_in_csv(self, tmp_path: Path) -> None:
        """Thứ tự xếp hạng (score cao → thấp) phải được giữ trong CSV."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["csv"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [_make_pair(f"p{i}", pair_penalty=float(i)) for i in range(5)]
        ranked = ranker.process(pairs)

        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)
        with paths[0].open() as fh:
            rows = list(csv.DictReader(fh))
        scores = [float(r["score"]) for r in rows]
        assert scores == sorted(scores, reverse=True)

    def test_hard_filtered_pairs_absent_from_all_reports(self, tmp_path: Path) -> None:
        """Pairs bị hard filter phải vắng mặt trong tất cả các loại báo cáo."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["csv", "json", "fasta"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [
            _make_pair("keep", left_tm=60.0),
            _make_pair("drop", snp_flags=["FAIL:3prime"]),
        ]
        ranked = ranker.process(pairs)

        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)
        path_map = {p.name: p for p in paths}

        # CSV
        with path_map["results_summary.csv"].open() as fh:
            csv_ids = [r["pair_id"] for r in csv.DictReader(fh)]
        assert "drop" not in csv_ids

        # JSON
        data = json.loads(path_map["results.json"].read_text())
        json_pair_ids = [p["pair_id"] for result in data for p in result.get("primer_pairs", [])]
        assert "drop" not in json_pair_ids

        # FASTA
        fasta_content = path_map["primers.fasta"].read_text()
        assert "drop" not in fasta_content

    def test_multiple_targets_in_report(self, tmp_path: Path) -> None:
        """Nhiều targets phải xuất hiện đầy đủ trong báo cáo."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["csv"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        results = []
        for tid in ["ACTB", "GAPDH", "TP53"]:
            pairs = [_make_pair(f"{tid}_p0"), _make_pair(f"{tid}_p1")]
            ranked = ranker.process(pairs)
            results.append(
                DesignResult(target=_make_target(tid), primer_pairs=ranked, status="success")
            )

        paths = reporter.generate(results, tmp_path)
        with paths[0].open() as fh:
            rows = list(csv.DictReader(fh))
        target_ids = {r["target_id"] for r in rows}
        assert {"ACTB", "GAPDH", "TP53"} == target_ids

    def test_failed_result_in_report(self, tmp_path: Path) -> None:
        """DesignResult với status=failed phải xuất hiện đúng trong báo cáo JSON."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["json"])
        reporter = ReportGenerator(cfg)

        results = [
            DesignResult(target=_make_target("ACTB"), primer_pairs=[], status="success"),
            DesignResult(
                target=_make_target("BAD"), primer_pairs=[], status="failed", error="Network error"
            ),
        ]
        paths = reporter.generate(results, tmp_path)
        data = json.loads(paths[0].read_text())
        statuses = {d["target"]["target_id"]: d["status"] for d in data}
        assert statuses["ACTB"] == "success"
        assert statuses["BAD"] == "failed"

    def test_top_n_respected_in_report(self, tmp_path: Path) -> None:
        """Top-N từ config phải được phản ánh trong báo cáo."""
        cfg = load_config(None)
        cfg.pipeline = PipelineConfig(top_n_pairs=2)
        cfg.output = OutputConfig(formats=["csv"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [_make_pair(f"p{i}") for i in range(10)]
        ranked = ranker.process(pairs)
        assert len(ranked) == 2

        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)
        with paths[0].open() as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 2

    def test_all_formats_generated(self, tmp_path: Path) -> None:
        """Tất cả 5 định dạng phải được tạo trong một lần gọi generate()."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["csv", "json", "xlsx", "fasta", "html"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [_make_pair("p0")]
        ranked = ranker.process(pairs)
        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]

        paths = reporter.generate(results, tmp_path)
        assert len(paths) == 5

        expected_names = {
            "results_summary.csv",
            "results.json",
            "results_detailed.xlsx",
            "primers.fasta",
            "report.html",
        }
        actual_names = {p.name for p in paths}
        assert expected_names == actual_names

    def test_fasta_sequences_match_primer_sequences(self, tmp_path: Path) -> None:
        """Sequences trong FASTA phải khớp chính xác với primer sequences."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["fasta"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        left_seq = "AAACCCGGGTTTAAACCCGGG"
        right_seq = "TTTGGGCCCAAATTTGGGCCC"
        pair = PrimerPair(
            pair_id="p_test",
            left_primer=Oligo(sequence=left_seq, tm=60.0),
            right_primer=Oligo(sequence=right_seq, tm=60.0),
            amplicon_size=150,
        )
        ranked = ranker.process([pair])
        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)
        content = paths[0].read_text()
        assert left_seq in content
        assert right_seq in content

    def test_html_contains_score(self, tmp_path: Path) -> None:
        """Score phải xuất hiện trong báo cáo HTML."""
        cfg = load_config(None)
        cfg.output = OutputConfig(formats=["html"])
        ranker = FilterRanker(cfg)
        reporter = ReportGenerator(cfg)

        pairs = [_make_pair("p0")]
        ranked = ranker.process(pairs)
        results = [DesignResult(target=_make_target(), primer_pairs=ranked, status="success")]
        paths = reporter.generate(results, tmp_path)
        content = paths[0].read_text(encoding="utf-8")
        # Score should appear in the HTML (formatted as float)
        score_str = f"{ranked[0].score:.1f}" if ranked else None
        if score_str:
            assert score_str in content


class TestFilterRankerIntegration:
    """Kiểm tra tích hợp chi tiết hơn cho FilterRanker."""

    def test_soft_flag_does_not_filter_pair(self) -> None:
        """Pair có soft flag (WARNING) không bị loại bỏ."""
        cfg = load_config(None)
        ranker = FilterRanker(cfg)
        pair = _make_pair("p0", snp_flags=["WARNING:mid_position"])
        result = ranker.process([pair])
        assert len(result) == 1

    def test_dg_warning_pair_still_ranked(self) -> None:
        """Pair có ΔG warning vẫn được xếp hạng (score thấp hơn pair sạch)."""
        cfg = load_config(None)
        ranker = FilterRanker(cfg)

        clean = _make_pair("clean")
        noisy = _make_pair("noisy")
        noisy.left_primer.self_any_th = -15.0

        result = ranker.process([clean, noisy])
        assert len(result) == 2
        assert result[0].score >= result[1].score

    def test_process_preserves_probe(self) -> None:
        """Probe phải được giữ nguyên sau khi process."""
        cfg = load_config(None)
        ranker = FilterRanker(cfg)
        probe = _make_oligo(seq="TTTTGGGGCCCCAAAATTTT", tm=67.0)
        pair = _make_pair("p0", probe=probe)
        result = ranker.process([pair])
        assert result[0].probe is not None
        assert result[0].probe.sequence == "TTTTGGGGCCCCAAAATTTT"

    def test_all_hard_filtered_returns_empty(self) -> None:
        cfg = load_config(None)
        ranker = FilterRanker(cfg)
        pairs = [_make_pair(f"p{i}", snp_flags=["FAIL:critical"]) for i in range(5)]
        result = ranker.process(pairs)
        assert result == []

    def test_score_increases_with_better_specificity(self) -> None:
        """Score phải tăng khi specificity_score cao hơn."""
        cfg = load_config(None)
        ranker = FilterRanker(cfg)

        class MockSpec:
            def __init__(self, s: float) -> None:
                self.specificity_score = s
                self.off_target_amplicons: list = []
                self.is_specific = True

        p_high = _make_pair("high", pair_penalty=1.0)
        p_low = _make_pair("low", pair_penalty=1.0)
        p_high.specificity_result = MockSpec(100.0)
        p_low.specificity_result = MockSpec(20.0)

        p_high.score = ranker.calculate_score(p_high)
        p_low.score = ranker.calculate_score(p_low)
        assert p_high.score > p_low.score
