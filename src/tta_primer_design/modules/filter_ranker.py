"""Module 10 — FilterRanker: lọc và xếp hạng primer pairs.

Scoring Algorithm (tổng 100 điểm):
    - Specificity score:         30 điểm
    - Thermodynamic score:       25 điểm (Tm, ΔG, GC clamp, no dimers)
    - SNP-free score:            20 điểm
    - Primer3 penalty (inverse): 15 điểm
    - Amplicon size score:       10 điểm

Filter Pipeline (theo thứ tự):
    1. Hard filter: off-target amplicon → LOẠI
    2. Hard filter: SNP critical position → LOẠI
    3. Hard filter: Tm out of range → LOẠI
    4. Soft filter: SNP warning → FLAG
    5. Soft filter: ΔG warning → FLAG
    6. Score & rank các pairs còn lại
    7. Trả về top-N pairs

TODO (Sprint 4):
    - Implement apply_hard_filters()
    - Implement apply_soft_filters()
    - Implement calculate_score()
    - Implement rank_pairs() và get_top_n()
"""

from __future__ import annotations

import logging

from tta_primer_design.config import AppConfig
from tta_primer_design.models import PrimerPair

logger = logging.getLogger("tta_primer_design.modules.filter_ranker")

# Trọng số điểm cho scoring
_WEIGHTS = {
    "specificity": 30.0,
    "thermodynamic": 25.0,
    "snp_free": 20.0,
    "primer3_penalty": 15.0,
    "amplicon_size": 10.0,
}


class FilterRanker:
    """Lọc và xếp hạng primer pairs theo chất lượng tổng hợp.

    Args:
        config: AppConfig đã load.

    Example::

        ranker = FilterRanker(config)
        final_pairs = ranker.process(primer_pairs)
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def process(self, primer_pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Chạy toàn bộ filter + rank pipeline.

        Args:
            primer_pairs: Danh sách PrimerPair đầu vào.

        Returns:
            Top-N PrimerPair đã được lọc và xếp hạng.
        """
        if not primer_pairs:
            return []
        pairs = self.apply_hard_filters(primer_pairs)
        logger.debug("After hard filters: %d/%d pairs remain", len(pairs), len(primer_pairs))
        pairs = self.apply_soft_filters(pairs)
        for pair in pairs:
            pair.score = self.calculate_score(pair)
        pairs = self.rank_pairs(pairs)
        return self.get_top_n(pairs)

    def apply_hard_filters(self, pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Áp dụng hard filters — loại bỏ các pair không đạt.

        Hard filters:
            - Off-target amplicon tồn tại
            - SNP ở vị trí critical (3' end primer, trong probe)
            - Tm ngoài khoảng cho phép

        Args:
            pairs: Danh sách PrimerPair.

        Returns:
            Danh sách sau khi loại bỏ các pair không đạt.
        """
        filtered = []
        for pair in pairs:
            if (
                pair.specificity_result is not None
                and hasattr(pair.specificity_result, "off_target_amplicons")
                and len(pair.specificity_result.off_target_amplicons)
                > self.config.filters.max_off_targets
            ):
                continue
            if any("FAIL" in flag for flag in pair.snp_flags):
                continue
            filtered.append(pair)
        return filtered

    def apply_soft_filters(self, pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Áp dụng soft filters — gán flag cảnh báo.

        Soft filters:
            - SNP warning (MAF thấp hoặc ở giữa primer)
            - ΔG warning

        Args:
            pairs: Danh sách PrimerPair.

        Returns:
            Danh sách với flags đã được gán.
        """
        for pair in pairs:
            if any("WARNING" in flag and "HAIRPIN" not in flag for flag in pair.snp_flags):
                if "SNP_WARNING" not in pair.snp_flags:
                    pair.snp_flags.append("SNP_WARNING")
            for oligo in (pair.left_primer, pair.right_primer):
                if oligo.hairpin_th > 20.0:
                    flag = f"HAIRPIN_WARNING:{oligo.hairpin_th:.1f}C"
                    if flag not in pair.snp_flags:
                        pair.snp_flags.append(flag)
        return pairs

    def calculate_score(self, pair: PrimerPair) -> float:
        """Tính điểm tổng hợp cho một PrimerPair (0–100).

        Args:
            pair: PrimerPair đã có đầy đủ thông tin.

        Returns:
            Điểm tổng (0–100).
        """
        score = 0.0

        # Specificity (30 pts)
        if pair.specificity_result is not None and hasattr(
            pair.specificity_result, "specificity_score"
        ):
            score += _WEIGHTS["specificity"] * (pair.specificity_result.specificity_score / 100.0)
        else:
            score += _WEIGHTS["specificity"] * 0.5

        # Thermodynamic (25 pts)
        tm_diff = (abs(pair.left_primer.tm - 60.0) + abs(pair.right_primer.tm - 60.0)) / 2.0
        tm_score = max(0.0, 1.0 - tm_diff / 5.0)
        max_hairpin = max(pair.left_primer.hairpin_th, pair.right_primer.hairpin_th)
        hairpin_score = max(0.0, 1.0 - max_hairpin / 45.0)
        score += _WEIGHTS["thermodynamic"] * ((tm_score + hairpin_score) / 2.0)

        # SNP-free (20 pts)
        if not pair.snp_flags:
            snp_score = 1.0
        elif all("WARNING" in f for f in pair.snp_flags):
            snp_score = 0.5
        else:
            snp_score = 0.0
        score += _WEIGHTS["snp_free"] * snp_score

        # Primer3 penalty inverse (15 pts)
        penalty_score = max(0.0, 1.0 - pair.pair_penalty / 5.0)
        score += _WEIGHTS["primer3_penalty"] * penalty_score

        # Amplicon size (10 pts)
        if pair.amplicon_size > 0:
            if 70 <= pair.amplicon_size <= 200:
                dist = abs(pair.amplicon_size - 125)
                size_score = max(0.0, 1.0 - dist / 75.0)
            else:
                size_score = 0.0
        else:
            size_score = 0.5
        score += _WEIGHTS["amplicon_size"] * size_score

        return round(min(100.0, max(0.0, score)), 2)

    def rank_pairs(self, pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Xếp hạng primer pairs theo score giảm dần.

        Args:
            pairs: Danh sách PrimerPair đã có score.

        Returns:
            Danh sách đã được sắp xếp (score cao nhất đứng đầu).
        """
        return sorted(pairs, key=lambda p: p.score, reverse=True)

    def get_top_n(self, pairs: list[PrimerPair], n: int | None = None) -> list[PrimerPair]:
        """Lấy top-N primer pairs.

        Args:
            pairs: Danh sách đã được rank.
            n: Số pair muốn lấy. None = dùng config.pipeline.top_n_pairs.

        Returns:
            Danh sách top-N pair.
        """
        top_n = n if n is not None else self.config.pipeline.top_n_pairs
        return pairs[:top_n]
