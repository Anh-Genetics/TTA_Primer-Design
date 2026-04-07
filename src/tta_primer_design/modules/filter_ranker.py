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

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("process chưa được implement (Sprint 4)")

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

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("apply_hard_filters chưa được implement (Sprint 4)")

    def apply_soft_filters(self, pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Áp dụng soft filters — gán flag cảnh báo.

        Soft filters:
            - SNP warning (MAF thấp hoặc ở giữa primer)
            - ΔG warning

        Args:
            pairs: Danh sách PrimerPair.

        Returns:
            Danh sách với flags đã được gán.

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("apply_soft_filters chưa được implement (Sprint 4)")

    def calculate_score(self, pair: PrimerPair) -> float:
        """Tính điểm tổng hợp cho một PrimerPair (0–100).

        Args:
            pair: PrimerPair đã có đầy đủ thông tin.

        Returns:
            Điểm tổng (0–100).

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("calculate_score chưa được implement (Sprint 4)")

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
