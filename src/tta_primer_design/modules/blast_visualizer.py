"""Module — BlastVisualizer: trực quan hoá kết quả BLAST specificity.

Cung cấp hai cách xuất:
    - :meth:`BlastVisualizer.print_ascii_chart` — in biểu đồ ASCII ra stdout.
    - :meth:`BlastVisualizer.save_plot` — lưu hình PNG bằng matplotlib (tuỳ chọn).

matplotlib là dependency *tuỳ chọn*; nếu chưa cài thì ``save_plot`` sẽ raise
``ImportError`` với hướng dẫn cài đặt.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tta_primer_design.modules.blast_specificity import BlastHit, SpecificityResult

logger = logging.getLogger("tta_primer_design.modules.blast_visualizer")

# Độ rộng cột mặc định cho ASCII chart
_BAR_WIDTH = 40
_TABLE_TOP_N = 10


def _mismatch_category(mismatches: int) -> str:
    """Phân loại số mismatch thành nhãn hiển thị."""
    if mismatches == 0:
        return "0"
    if mismatches == 1:
        return "1"
    if mismatches == 2:
        return "2"
    return "3+"


def _count_by_mismatch(hits: list[BlastHit]) -> dict[str, int]:
    """Đếm số hit theo nhóm mismatch (0, 1, 2, 3+)."""
    counts: dict[str, int] = {"0": 0, "1": 0, "2": 0, "3+": 0}
    for hit in hits:
        counts[_mismatch_category(hit.mismatches)] += 1
    return counts


class BlastVisualizer:
    """Trực quan hoá kết quả BLAST specificity cho một cặp primer.

    Args:
        result: :class:`~tta_primer_design.modules.blast_specificity.SpecificityResult`
            đã chạy BLAST.
        pair_name: Tên/ID cặp mồi (mặc định lấy từ ``result.primer_pair_id``).

    Example::

        viz = BlastVisualizer(spec_result, pair_name="ACTB_pair1")
        viz.print_ascii_chart()
        viz.save_plot(Path("results/blast_plot.png"))
    """

    def __init__(self, result: SpecificityResult, pair_name: str | None = None) -> None:
        self.result = result
        self.pair_name = pair_name or result.primer_pair_id

    # ------------------------------------------------------------------
    # ASCII chart
    # ------------------------------------------------------------------

    def print_ascii_chart(self) -> None:
        """In biểu đồ ASCII trực quan hoá kết quả BLAST ra stdout.

        Bao gồm:
        - Thanh ngang (horizontal bar) thể hiện số hit mỗi primer.
        - Bảng top BLAST hits (subject_id, identity%, mismatches, e-value).
        - Off-target amplicons nếu có.
        """
        sep = "─" * 60
        print()
        print(f"  BLAST Visualization — {self.pair_name}")
        print(sep)

        left_hits = self.result.blast_hits_left
        right_hits = self.result.blast_hits_right

        # ── Hit count bars ──────────────────────────────────────────
        print("  Hit counts per primer:")
        self._print_bar("Left ", len(left_hits))
        self._print_bar("Right", len(right_hits))
        print()

        # ── Mismatch distribution ────────────────────────────────────
        print("  Mismatch distribution:")
        left_mm = _count_by_mismatch(left_hits)
        right_mm = _count_by_mismatch(right_hits)
        print(f"  {'Cat':>3}  {'Left':>6}  {'Right':>6}")
        print(f"  {'───':>3}  {'──────':>6}  {'──────':>6}")
        for cat in ("0", "1", "2", "3+"):
            print(f"  {cat:>3}  {left_mm[cat]:>6}  {right_mm[cat]:>6}")
        print()

        # ── Top hits table ────────────────────────────────────────────
        for label, hits in (("Left primer", left_hits), ("Right primer", right_hits)):
            print(f"  Top BLAST hits — {label}:")
            if not hits:
                print("    (no hits)")
            else:
                header = (
                    f"  {'#':>3}  {'Subject ID':<20}  {'Identity%':>9}  {'MM':>4}  {'E-value':>12}"
                )
                print(header)
                print(
                    f"  {'───':>3}  {'──────────────────────':20}  {'─────────':>9}  {'────':>4}  {'────────────':>12}"
                )
                for i, hit in enumerate(hits[:_TABLE_TOP_N], 1):
                    subj = hit.subject_id[:20]
                    eval_str = f"{hit.evalue:.2e}" if hit.evalue > 0 else "0.0"
                    print(
                        f"  {i:>3}  {subj:<20}  {hit.identity:>8.1f}%"
                        f"  {hit.mismatches:>4}  {eval_str:>12}"
                    )
                if len(hits) > _TABLE_TOP_N:
                    print(f"    … và {len(hits) - _TABLE_TOP_N} hit khác")
            print()

        # ── Off-target amplicons ──────────────────────────────────────
        off = self.result.off_target_amplicons
        if off:
            print(f"  ⚠️  Off-target amplicons ({len(off)} tổng cộng):")
            print(f"  {'#':>3}  {'Subject ID':<20}  {'Size (bp)':>10}")
            print(f"  {'───':>3}  {'──────────────────────':20}  {'──────────':>10}")
            for i, amp in enumerate(off[:10], 1):
                print(f"  {i:>3}  {amp.subject_id[:20]:<20}  {amp.amplicon_size:>10}")
            if len(off) > 10:
                print(f"    … và {len(off) - 10} amplicon khác")
        else:
            print("  ✅ Không có off-target amplicon.")

        print(sep)
        print()

    @staticmethod
    def _print_bar(label: str, count: int, max_count: int = 50) -> None:
        """In một thanh ASCII tương ứng với ``count``."""
        filled = min(int(count / max(max_count, 1) * _BAR_WIDTH), _BAR_WIDTH)
        bar = "█" * filled + "░" * (_BAR_WIDTH - filled)
        print(f"  {label}: |{bar}| {count}")

    # ------------------------------------------------------------------
    # matplotlib plot
    # ------------------------------------------------------------------

    def save_plot(self, output_path: Path) -> Path:
        """Tạo biểu đồ matplotlib và lưu ra file PNG.

        Các subplot:
        - Subplot 1: Bar chart số hit theo nhóm mismatch (0/1/2/3+) cho mỗi primer.
        - Subplot 2: Scatter/lollipop e-value theo chỉ số hit cho cả hai primer.
        - Subplot 3 (chỉ khi có off-target amplicon): Bar chart kích thước amplicon.

        Args:
            output_path: Đường dẫn file PNG đầu ra.

        Returns:
            ``output_path`` sau khi lưu thành công.

        Raises:
            ImportError: Nếu matplotlib chưa được cài đặt.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "matplotlib is required for save_plot(). "
                "Run: pip install 'tta-primer-design[viz]'"
            ) from exc

        left_hits = self.result.blast_hits_left
        right_hits = self.result.blast_hits_right
        off = self.result.off_target_amplicons
        has_off_target = bool(off)

        n_subplots = 3 if has_off_target else 2
        fig, axes = plt.subplots(1, n_subplots, figsize=(6 * n_subplots, 5))
        if n_subplots == 1:
            axes = [axes]

        fig.suptitle(f"BLAST Specificity — {self.pair_name}", fontsize=13, fontweight="bold")

        # ── Subplot 1: mismatch category bar chart ────────────────────
        ax1 = axes[0]
        categories = ["0", "1", "2", "3+"]
        left_counts = [_count_by_mismatch(left_hits)[c] for c in categories]
        right_counts = [_count_by_mismatch(right_hits)[c] for c in categories]

        x = range(len(categories))
        width = 0.35
        ax1.bar([xi - width / 2 for xi in x], left_counts, width, label="Left", color="#4C72B0")
        ax1.bar([xi + width / 2 for xi in x], right_counts, width, label="Right", color="#DD8452")
        ax1.set_xticks(list(x))
        ax1.set_xticklabels(categories)
        ax1.set_xlabel("Mismatches")
        ax1.set_ylabel("Hit count")
        ax1.set_title("Hits by mismatch category")
        ax1.legend()

        # ── Subplot 2: e-value lollipop ───────────────────────────────
        ax2 = axes[1]
        for hits, color, label in (
            (left_hits, "#4C72B0", "Left"),
            (right_hits, "#DD8452", "Right"),
        ):
            if not hits:
                continue
            xs = list(range(len(hits)))
            # Chuyển e-value sang -log10; e-value = 0 → dùng giá trị rất nhỏ
            ys = [-math.log10(max(h.evalue, 1e-200)) for h in hits]
            ax2.scatter(xs, ys, color=color, s=20, zorder=3, label=label)
            for xi, yi in zip(xs, ys, strict=True):
                ax2.plot([xi, xi], [0, yi], color=color, linewidth=0.5, alpha=0.6)

        ax2.set_xlabel("Hit index")
        ax2.set_ylabel("-log₁₀(E-value)")
        ax2.set_title("BLAST hit E-values")
        ax2.legend()

        # ── Subplot 3: off-target amplicon sizes (nếu có) ─────────────
        if has_off_target:
            ax3 = axes[2]
            sizes = [amp.amplicon_size for amp in off]
            ax3.bar(range(len(sizes)), sizes, color="#55A868")
            ax3.set_xlabel("Amplicon index")
            ax3.set_ylabel("Size (bp)")
            ax3.set_title(f"Off-target amplicons (n={len(sizes)})")

        plt.tight_layout()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

        logger.info("Saved BLAST plot to %s", output_path)
        return output_path
