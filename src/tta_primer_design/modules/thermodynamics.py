"""Module 9 — Thermodynamics: tính toán nhiệt động học cho primer/probe.

Công cụ sử dụng:
    - primer3-py: Tm, ΔG self-dimer, hairpin (ưu tiên)
    - ViennaRNA (tùy chọn): RNA/DNA secondary structure chi tiết

Thông số tính toán:
    1. Tm (SantaLucia 1998 — nearest-neighbor model)
    2. ΔG hairpin (secondary structure)
    3. ΔG self-dimer
    4. ΔG hetero-dimer (left vs right, primer vs probe)
    5. ΔG 3'-end stability
    6. GC clamp check (2–3 G/C ở 3' cuối)
    7. Repeat runs (≤4 identical bases liên tiếp)
    8. Dinucleotide repeats (ví dụ: ATATAT...)

Ngưỡng tham chiếu:
    Tm primer:      58–62°C
    Tm probe:       65–72°C
    ΔG hairpin:     > -9 kcal/mol
    ΔG self-dimer:  > -9 kcal/mol
    GC clamp:       1–3 G/C trong 5 bases cuối 3'

TODO (Sprint 2):
    - Implement calculate_tm() dùng primer3-py hoặc công thức Wallace/SantaLucia
    - Implement calculate_hairpin_dg() và calculate_dimer_dg()
    - Implement check_gc_clamp() và check_repeat_runs()
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from tta_primer_design.config import AppConfig
from tta_primer_design.models import Oligo, PrimerPair

logger = logging.getLogger("tta_primer_design.modules.thermodynamics")


@dataclass
class ThermoProfile:
    """Profile nhiệt động học đầy đủ cho một oligo.

    Attributes:
        sequence: Chuỗi oligo.
        tm: Nhiệt độ nóng chảy (°C).
        gc_percent: Tỷ lệ GC (%).
        hairpin_dg: ΔG hairpin (kcal/mol).
        self_dimer_dg: ΔG self-dimer (kcal/mol).
        end_stability_dg: ΔG 3'-end stability (kcal/mol).
        gc_clamp_ok: True nếu có GC clamp hợp lệ.
        repeat_ok: True nếu không có repeat runs > max_run.
        pass_all: True nếu tất cả thông số đều hợp lệ.
    """

    sequence: str
    tm: float = 0.0
    gc_percent: float = 0.0
    hairpin_dg: float = 0.0
    self_dimer_dg: float = 0.0
    end_stability_dg: float = 0.0
    gc_clamp_ok: bool = True
    repeat_ok: bool = True
    pass_all: bool = True


class Thermodynamics:
    """Tính toán và validate thông số nhiệt động học.

    Args:
        config: AppConfig đã load.

    Example::

        thermo = Thermodynamics(config)
        pairs = thermo.validate_all(primer_pairs)
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def validate_all(self, primer_pairs: list[PrimerPair]) -> list[PrimerPair]:
        """Validate nhiệt động học cho tất cả primer pairs.

        Args:
            primer_pairs: Danh sách PrimerPair.

        Returns:
            PrimerPair đã được cập nhật thông tin nhiệt động học.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("validate_all chưa được implement (Sprint 2)")

    def calculate_tm(self, sequence: str, params: dict | None = None) -> float:
        """Tính Tm (SantaLucia 1998).

        Args:
            sequence: Chuỗi oligo (5'→3').
            params: Dict tham số (salt, DNA conc, v.v.)

        Returns:
            Tm (°C).

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("calculate_tm chưa được implement (Sprint 2)")

    def calculate_hairpin_dg(self, sequence: str) -> float:
        """Tính ΔG cấu trúc hairpin.

        Args:
            sequence: Chuỗi oligo.

        Returns:
            ΔG hairpin (kcal/mol). Âm = có hairpin.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("calculate_hairpin_dg chưa được implement (Sprint 2)")

    def calculate_dimer_dg(self, seq1: str, seq2: str) -> float:
        """Tính ΔG hetero-dimer giữa 2 oligo.

        Args:
            seq1: Chuỗi oligo thứ nhất.
            seq2: Chuỗi oligo thứ hai.

        Returns:
            ΔG dimer (kcal/mol).

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("calculate_dimer_dg chưa được implement (Sprint 2)")

    def check_gc_clamp(
        self, sequence: str, window: int = 5, min_gc: int = 1, max_gc: int = 3
    ) -> bool:
        """Kiểm tra GC clamp ở 3' end.

        Args:
            sequence: Chuỗi oligo.
            window: Số base ở 3' end cần kiểm tra.
            min_gc: Số G/C tối thiểu.
            max_gc: Số G/C tối đa.

        Returns:
            True nếu GC clamp hợp lệ.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("check_gc_clamp chưa được implement (Sprint 2)")

    def check_repeat_runs(self, sequence: str, max_run: int = 4) -> bool:
        """Kiểm tra không có repeat runs quá dài.

        Args:
            sequence: Chuỗi oligo.
            max_run: Số base giống nhau liên tiếp tối đa cho phép.

        Returns:
            True nếu không có repeat runs > max_run.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("check_repeat_runs chưa được implement (Sprint 2)")

    def full_thermodynamic_profile(self, oligo: Oligo) -> ThermoProfile:
        """Tính toán đầy đủ profile nhiệt động học cho một oligo.

        Args:
            oligo: Oligo cần tính toán.

        Returns:
            ThermoProfile.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("full_thermodynamic_profile chưa được implement (Sprint 2)")
