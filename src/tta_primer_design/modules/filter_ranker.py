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

# Ngưỡng nhiệt động học
_TM_MIN_HARD = 45.0  # Tm tối thiểu tuyệt đối (°C)
_TM_MAX_HARD = 80.0  # Tm tối đa tuyệt đối (°C)
_TM_OPTIMAL_MIN = 55.0  # Tm tối ưu tối thiểu (°C)
_TM_OPTIMAL_MAX = 68.0  # Tm tối ưu tối đa (°C)
_DG_WARNING_THRESHOLD = -9.0  # ΔG cảnh báo (kcal/mol)

# Ngưỡng amplicon size hợp lý cho qPCR
_AMPLICON_QPCR_MIN = 80
_AMPLICON_QPCR_MAX = 150
_AMPLICON_PCR_MIN = 100
_AMPLICON_PCR_MAX = 2000

# Penalty tối đa Primer3 để chuẩn hoá
_MAX_PRIMER3_PENALTY = 20.0


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
        logger.info("FilterRanker: processing %d primer pairs", len(primer_pairs))

        pairs = self.apply_hard_filters(primer_pairs)
        logger.info("After hard filters: %d pairs remain", len(pairs))

        pairs = self.apply_soft_filters(pairs)

        for pair in pairs:
            pair.score = self.calculate_score(pair)

        pairs = self.rank_pairs(pairs)
        pairs = self.get_top_n(pairs)

        logger.info("FilterRanker: returning top %d pairs", len(pairs))
        return pairs

    def apply_hard_filters(self, pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Áp dụng hard filters — loại bỏ các pair không đạt.

        Hard filters:
            - Off-target amplicon tồn tại (vượt quá max_off_targets)
            - SNP ở vị trí critical (3' end primer, trong probe) — flag "FAIL"
            - Tm ngoài khoảng cho phép tuyệt đối

        Args:
            pairs: Danh sách PrimerPair.

        Returns:
            Danh sách sau khi loại bỏ các pair không đạt.
        """
        passed: list[PrimerPair] = []
        for pair in pairs:
            reason = self._hard_filter_reason(pair)
            if reason:
                logger.debug("Hard-filtered pair '%s': %s", pair.pair_id, reason)
            else:
                passed.append(pair)
        return passed

    def _hard_filter_reason(self, pair: PrimerPair) -> str | None:
        """Trả về lý do loại bỏ nếu pair không đạt hard filter, ngược lại None."""
        # 1. Off-target amplicons
        if pair.specificity_result is not None:
            off_targets = getattr(pair.specificity_result, "off_target_amplicons", [])
            if len(off_targets) > self.config.filters.max_off_targets:
                return (
                    f"off_target amplicons={len(off_targets)} > "
                    f"max={self.config.filters.max_off_targets}"
                )

        # 2. Critical SNP flags
        for flag in pair.snp_flags:
            if "FAIL" in flag.upper():
                return f"critical SNP flag: {flag}"

        # 3. Tm out of hard range
        for label, oligo in [
            ("left_primer", pair.left_primer),
            ("right_primer", pair.right_primer),
        ]:
            if oligo.tm > 0 and not (_TM_MIN_HARD <= oligo.tm <= _TM_MAX_HARD):
                return (
                    f"{label} Tm={oligo.tm:.1f}°C out of hard range [{_TM_MIN_HARD},{_TM_MAX_HARD}]"
                )

        if pair.probe is not None and pair.probe.tm > 0:
            if not (_TM_MIN_HARD <= pair.probe.tm <= _TM_MAX_HARD):
                return f"probe Tm={pair.probe.tm:.1f}°C out of hard range"

        return None

    def apply_soft_filters(self, pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Áp dụng soft filters — gán flag cảnh báo.

        Soft filters:
            - SNP warning (MAF thấp hoặc ở giữa primer)
            - ΔG warning (self-dimer/hairpin quá mạnh)

        Args:
            pairs: Danh sách PrimerPair.

        Returns:
            Danh sách với flags đã được gán (không loại bỏ pair nào).
        """
        for pair in pairs:
            # SNP WARNING flags (không phải FAIL — đã xử lý ở hard filter)
            for flag in pair.snp_flags:
                if "WARNING" in flag.upper() and "WARNING:snp_warning" not in pair.snp_flags:
                    pair.snp_flags.append("WARNING:snp_warning")
                    break

            # ΔG warning
            dg_warn = False
            for oligo in [pair.left_primer, pair.right_primer]:
                if oligo.self_any_th < _DG_WARNING_THRESHOLD:
                    dg_warn = True
                    break
                if oligo.hairpin_th < _DG_WARNING_THRESHOLD:
                    dg_warn = True
                    break
            if pair.probe is not None:
                if pair.probe.self_any_th < _DG_WARNING_THRESHOLD:
                    dg_warn = True
                if pair.probe.hairpin_th < _DG_WARNING_THRESHOLD:
                    dg_warn = True

            if dg_warn and "WARNING:dg_high" not in pair.snp_flags:
                pair.snp_flags.append("WARNING:dg_high")

        return pairs

    def calculate_score(self, pair: PrimerPair) -> float:
        """Tính điểm tổng hợp cho một PrimerPair (0–100).

        Args:
            pair: PrimerPair đã có đầy đủ thông tin.

        Returns:
            Điểm tổng (0–100).
        """
        score = 0.0

        # --- Specificity score (30 điểm) ---
        specificity_raw = 50.0  # mặc định khi chưa có kết quả BLAST
        if pair.specificity_result is not None:
            specificity_raw = float(getattr(pair.specificity_result, "specificity_score", 50.0))
        score += _WEIGHTS["specificity"] * (specificity_raw / 100.0)

        # --- Thermodynamic score (25 điểm) ---
        thermo_score = self._thermodynamic_score(pair)
        score += _WEIGHTS["thermodynamic"] * thermo_score

        # --- SNP-free score (20 điểm) ---
        snp_score = self._snp_score(pair)
        score += _WEIGHTS["snp_free"] * snp_score

        # --- Primer3 penalty inverse (15 điểm) ---
        # Penalty thấp → điểm cao
        penalty = min(pair.pair_penalty, _MAX_PRIMER3_PENALTY)
        penalty_score = 1.0 - (penalty / _MAX_PRIMER3_PENALTY)
        score += _WEIGHTS["primer3_penalty"] * penalty_score

        # --- Amplicon size score (10 điểm) ---
        amplicon_score = self._amplicon_size_score(pair)
        score += _WEIGHTS["amplicon_size"] * amplicon_score

        return round(min(max(score, 0.0), 100.0), 4)

    def _thermodynamic_score(self, pair: PrimerPair) -> float:
        """Trả về thành phần nhiệt động học (0.0–1.0)."""
        points = 0.0
        total = 0.0

        oligos = [pair.left_primer, pair.right_primer]
        if pair.probe is not None:
            oligos.append(pair.probe)

        for oligo in oligos:
            total += 1.0
            if oligo.tm > 0:
                if _TM_OPTIMAL_MIN <= oligo.tm <= _TM_OPTIMAL_MAX:
                    points += 1.0
                else:
                    # Partial credit khi gần khoảng tối ưu
                    dist = min(
                        abs(oligo.tm - _TM_OPTIMAL_MIN),
                        abs(oligo.tm - _TM_OPTIMAL_MAX),
                    )
                    points += max(0.0, 1.0 - dist / 10.0)
            else:
                # Tm chưa tính — cho điểm trung bình
                points += 0.5

        if total == 0:
            return 0.5
        return points / total

    def _snp_score(self, pair: PrimerPair) -> float:
        """Trả về thành phần SNP-free (0.0–1.0)."""
        if not pair.snp_flags:
            return 1.0
        warnings = sum(1 for f in pair.snp_flags if "WARNING" in f.upper())
        # Mỗi warning giảm 0.3; không thể âm
        return max(0.0, 1.0 - warnings * 0.3)

    def _amplicon_size_score(self, pair: PrimerPair) -> float:
        """Trả về thành phần amplicon size (0.0–1.0)."""
        if pair.amplicon_size == 0:
            return 0.5  # Chưa có thông tin — cho điểm trung bình

        mode = self.config.pipeline.mode
        if mode in ("qpcr", "taqman", "sybr"):
            opt_min, opt_max = _AMPLICON_QPCR_MIN, _AMPLICON_QPCR_MAX
        else:
            opt_min, opt_max = _AMPLICON_PCR_MIN, _AMPLICON_PCR_MAX

        if opt_min <= pair.amplicon_size <= opt_max:
            return 1.0
        # Partial credit khi gần khoảng tối ưu
        dist = min(
            abs(pair.amplicon_size - opt_min),
            abs(pair.amplicon_size - opt_max),
        )
        return max(0.0, 1.0 - dist / opt_max)

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
