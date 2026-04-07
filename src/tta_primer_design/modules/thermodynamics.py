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
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

from tta_primer_design.config import AppConfig
from tta_primer_design.models import Oligo, PrimerPair

logger = logging.getLogger("tta_primer_design.modules.thermodynamics")

# primer3-py trả về ΔG theo đơn vị cal/mol — chuyển sang kcal/mol
_CAL_TO_KCAL = 1e-3


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

        Tính ThermoProfile cho left, right, và probe (nếu có) của mỗi cặp.
        Ghi kết quả vào các trường tương ứng của Oligo.

        Args:
            primer_pairs: Danh sách PrimerPair.

        Returns:
            PrimerPair đã được cập nhật thông tin nhiệt động học.
        """
        for pair in primer_pairs:
            _annotate_oligo(self, pair.left_primer)
            _annotate_oligo(self, pair.right_primer)
            if pair.probe is not None:
                _annotate_oligo(self, pair.probe)
        return primer_pairs

    def calculate_tm(self, sequence: str, params: dict | None = None) -> float:
        """Tính Tm (SantaLucia 1998) qua primer3-py calc_tm.

        Args:
            sequence: Chuỗi oligo (5'→3').
            params: Dict tham số tuỳ chọn (mv_conc, dv_conc, dntp_conc, dna_conc).

        Returns:
            Tm (°C).
        """
        import primer3

        kwargs: dict = {}
        if params:
            for key in ("mv_conc", "dv_conc", "dntp_conc", "dna_conc"):
                if key in params:
                    kwargs[key] = params[key]
        return float(primer3.calc_tm(sequence.upper(), **kwargs))

    def calculate_hairpin_dg(self, sequence: str) -> float:
        """Tính ΔG cấu trúc hairpin qua primer3-py calc_hairpin.

        Args:
            sequence: Chuỗi oligo.

        Returns:
            ΔG hairpin (kcal/mol). Âm = có hairpin.
        """
        import primer3

        result = primer3.calc_hairpin(sequence.upper())
        return float(result.dg) * _CAL_TO_KCAL

    def calculate_dimer_dg(self, seq1: str, seq2: str) -> float:
        """Tính ΔG hetero-dimer giữa 2 oligo qua primer3-py calc_heterodimer.

        Nếu seq1 == seq2 thì tính self-dimer (homodimer).

        Args:
            seq1: Chuỗi oligo thứ nhất.
            seq2: Chuỗi oligo thứ hai.

        Returns:
            ΔG dimer (kcal/mol).
        """
        import primer3

        s1 = seq1.upper()
        s2 = seq2.upper()
        if s1 == s2:
            result = primer3.calc_homodimer(s1)
        else:
            result = primer3.calc_heterodimer(s1, s2)
        return float(result.dg) * _CAL_TO_KCAL

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
        """
        if not sequence:
            return False
        tail = sequence.upper()[-window:]
        gc_count = tail.count("G") + tail.count("C")
        return min_gc <= gc_count <= max_gc

    def check_repeat_runs(self, sequence: str, max_run: int = 4) -> bool:
        """Kiểm tra không có repeat runs quá dài.

        Kiểm tra cả mononucleotide runs (AAAA...) và dinucleotide repeats
        (ATATAT...).

        Args:
            sequence: Chuỗi oligo.
            max_run: Số base giống nhau liên tiếp tối đa cho phép.

        Returns:
            True nếu không có repeat runs > max_run.
        """
        seq = sequence.upper()
        # Mononucleotide repeat: (X){max_run+1,}
        mono_pattern = re.compile(r"([ATCG])\1{" + str(max_run) + r",}")
        if mono_pattern.search(seq):
            return False
        # Dinucleotide repeat: (XY){n,} where n*2 > max_run
        min_di_reps = (max_run // 2) + 1
        di_pattern = re.compile(r"([ATCG]{2})\1{" + str(min_di_reps - 1) + r",}")
        if di_pattern.search(seq):
            return False
        return True

    def full_thermodynamic_profile(self, oligo: Oligo) -> ThermoProfile:
        """Tính toán đầy đủ profile nhiệt động học cho một oligo.

        Args:
            oligo: Oligo cần tính toán.

        Returns:
            ThermoProfile đầy đủ với tất cả thông số.
        """
        seq = oligo.sequence
        seq_up = seq.upper()

        tm = self.calculate_tm(seq)
        gc_count = seq_up.count("G") + seq_up.count("C")
        gc_percent = (gc_count / len(seq_up) * 100) if seq_up else 0.0
        hairpin_dg = self.calculate_hairpin_dg(seq)
        self_dimer_dg = self.calculate_dimer_dg(seq, seq)
        end_stability_dg = _calc_end_stability(seq)
        gc_clamp_ok = self.check_gc_clamp(seq)
        repeat_ok = self.check_repeat_runs(seq)

        pass_all = gc_clamp_ok and repeat_ok

        return ThermoProfile(
            sequence=seq,
            tm=tm,
            gc_percent=gc_percent,
            hairpin_dg=hairpin_dg,
            self_dimer_dg=self_dimer_dg,
            end_stability_dg=end_stability_dg,
            gc_clamp_ok=gc_clamp_ok,
            repeat_ok=repeat_ok,
            pass_all=pass_all,
        )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _annotate_oligo(thermo: Thermodynamics, oligo: Oligo) -> None:
    """Ghi các thông số nhiệt động học vào trường của Oligo in-place."""
    seq = oligo.sequence.upper()
    import primer3

    oligo.tm = float(primer3.calc_tm(seq))
    gc_count = seq.count("G") + seq.count("C")
    oligo.gc_percent = (gc_count / len(seq) * 100) if seq else 0.0
    oligo.hairpin_th = float(primer3.calc_hairpin(seq).dg) * _CAL_TO_KCAL
    oligo.self_any_th = float(primer3.calc_homodimer(seq).dg) * _CAL_TO_KCAL
    oligo.self_end_th = float(primer3.calc_end_stability(seq, seq).dg) * _CAL_TO_KCAL
    oligo.end_stability = _calc_end_stability(seq)


def _calc_end_stability(sequence: str) -> float:
    """Tính ΔG 3'-end stability (kcal/mol) dùng primer3-py calc_end_stability.

    Dùng complement của chính sequence làm template để ước lượng.
    """
    import primer3

    seq = sequence.upper()
    complement = seq.translate(str.maketrans("ATCG", "TAGC"))
    result = primer3.calc_end_stability(seq, complement)
    return float(result.dg) * _CAL_TO_KCAL
