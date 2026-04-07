"""Module 7 — ProbeDesigner: thiết kế probe TaqMan / hydrolysis probe.

Quy tắc thiết kế TaqMan probe:
    1. Probe nằm TRONG amplicon, không overlap với primer
    2. Tm probe > Tm primer + 8–10°C
    3. Tránh G ở đầu 5' (dập tắt fluorophore reporter)
    4. Nếu G > C trong probe: dùng reverse complement strand
    5. Không có polyG/polyC quá 4 nucleotide liên tiếp
    6. ΔG hairpin < -2 kcal/mol (tránh secondary structure)
    7. Không overlap với vị trí primer

Chế độ hỗ trợ (theo kế hoạch):
    - TaqMan (5'-nuclease): FAM-probe-BHQ1/MGB
    - SYBR Green: chỉ cần primer (không probe)
    - Molecular Beacon: stem-loop (stub)
    - Scorpion probe (stub)
    - LNA-modified probe (stub)

TODO (Sprint 3):
    - Implement ``design_taqman_probe()``
    - Implement ``check_probe_rules()``
    - Integrate với thermodynamics.py để tính Tm probe
"""

from __future__ import annotations

import logging

from tta_primer_design.config import AppConfig
from tta_primer_design.models import Oligo, PrimerPair, ProcessedSequence

logger = logging.getLogger("tta_primer_design.modules.probe_designer")


class ProbeDesigner:
    """Thiết kế probe TaqMan cho các primer pair.

    Args:
        config: AppConfig đã load.

    Example::

        designer = ProbeDesigner(config)
        pairs_with_probe = designer.design(primer_pairs, processed_seq)
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def design(
        self,
        primer_pairs: list[PrimerPair],
        processed_seq: ProcessedSequence,
    ) -> list[PrimerPair]:
        """Thiết kế probe cho tất cả primer pairs.

        Args:
            primer_pairs: Danh sách PrimerPair (đã có left/right primer).
            processed_seq: ProcessedSequence chứa sequence gốc.

        Returns:
            Danh sách PrimerPair đã được gán probe.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("design chưa được implement (Sprint 3)")

    def design_taqman_probe(
        self,
        amplicon: str,
        pair: PrimerPair,
    ) -> Oligo | None:
        """Thiết kế TaqMan probe cho một amplicon.

        Args:
            amplicon: Chuỗi amplicon (giữa 2 primer).
            pair: PrimerPair tương ứng.

        Returns:
            Oligo probe hoặc None nếu không thiết kế được.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("design_taqman_probe chưa được implement (Sprint 3)")

    def check_probe_rules(self, probe: Oligo) -> tuple[bool, list[str]]:
        """Kiểm tra probe theo các quy tắc TaqMan.

        Args:
            probe: Oligo probe cần kiểm tra.

        Returns:
            (pass_flag, list_of_violations)

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("check_probe_rules chưa được implement (Sprint 3)")

    def select_best_strand(self, amplicon: str, pos: int, length: int) -> str:
        """Chọn strand tốt hơn (ít G hơn) cho probe.

        Args:
            amplicon: Chuỗi amplicon.
            pos: Vị trí bắt đầu probe trên amplicon (0-based).
            length: Độ dài probe.

        Returns:
            Chuỗi probe trên strand tốt hơn.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("select_best_strand chưa được implement (Sprint 3)")

    def calculate_probe_tm(self, sequence: str, params: dict) -> float:
        """Tính Tm của probe.

        Args:
            sequence: Chuỗi probe.
            params: Dict tham số nhiệt động học.

        Returns:
            Tm (°C).

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("calculate_probe_tm chưa được implement (Sprint 3)")
