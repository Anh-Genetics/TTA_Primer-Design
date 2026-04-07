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

_COMPLEMENT = str.maketrans("ATCGatcg", "TAGCtagc")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _extract_amplicon(sequence: str, pair: PrimerPair) -> str | None:
    """Extract amplicon sequence from template using primer positions."""
    left_start = pair.left_primer.start
    left_len = pair.left_primer.length
    right_start = pair.right_primer.start
    right_len = pair.right_primer.length
    if left_len == 0 or right_len == 0:
        return None
    amp_end = right_start + 1
    if amp_end <= left_start:
        return None
    return sequence[left_start:amp_end]


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
        """
        mode = self.config.pipeline.mode
        if mode not in ("taqman", "qpcr"):
            return primer_pairs
        sequence = processed_seq.sequence
        for pair in primer_pairs:
            if pair.probe is not None:
                continue
            amplicon = _extract_amplicon(sequence, pair)
            if amplicon:
                probe = self.design_taqman_probe(amplicon, pair)
                if probe is not None:
                    pair.probe = probe
        return primer_pairs

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
        """
        left_len = pair.left_primer.length
        right_len = pair.right_primer.length
        amp_len = len(amplicon)
        avail_start = left_len
        avail_end = amp_len - right_len
        if avail_end - avail_start < 18:
            return None
        best_probe = None
        best_tm_diff = float("inf")
        avg_primer_tm = (pair.left_primer.tm + pair.right_primer.tm) / 2.0
        target_probe_tm = avg_primer_tm + 9.0
        for probe_len in range(18, 31):
            for pos in range(avail_start, avail_end - probe_len + 1):
                probe_seq = self.select_best_strand(amplicon, pos, probe_len)
                try:
                    tm = self.calculate_probe_tm(probe_seq, {})
                except Exception:
                    continue
                ok, _ = self.check_probe_rules(Oligo(sequence=probe_seq, tm=tm))
                if ok:
                    tm_diff = abs(tm - target_probe_tm)
                    if tm_diff < best_tm_diff:
                        best_tm_diff = tm_diff
                        best_probe = Oligo(
                            sequence=probe_seq,
                            start=pair.left_primer.start + pos,
                            tm=tm,
                            gc_percent=(probe_seq.count("G") + probe_seq.count("C"))
                            / len(probe_seq)
                            * 100.0,
                        )
        return best_probe

    def check_probe_rules(self, probe: Oligo) -> tuple[bool, list[str]]:
        """Kiểm tra probe theo các quy tắc TaqMan.

        Args:
            probe: Oligo probe cần kiểm tra.

        Returns:
            (pass_flag, list_of_violations)
        """
        violations = []
        seq = probe.sequence.upper()
        if seq.startswith("G"):
            violations.append("5_prime_G")
        for base in ("G", "C"):
            if base * 5 in seq:
                violations.append(f"poly_{base}_run")
        gc_pct = (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0
        if not (40 <= gc_pct <= 65):
            violations.append(f"gc_out_of_range_{gc_pct:.0f}pct")
        return len(violations) == 0, violations

    def select_best_strand(self, amplicon: str, pos: int, length: int) -> str:
        """Chọn strand tốt hơn (ít G hơn) cho probe.

        Args:
            amplicon: Chuỗi amplicon.
            pos: Vị trí bắt đầu probe trên amplicon (0-based).
            length: Độ dài probe.

        Returns:
            Chuỗi probe trên strand tốt hơn.
        """
        forward = amplicon[pos : pos + length].upper()
        reverse = _reverse_complement(forward)
        # In a reverse complement: reverse.count("G") == forward.count("C")
        # So forward.count("G") > forward.count("C") ⟺ G(fwd) > G(rev)
        g_fwd = forward.count("G")
        g_rev = reverse.count("G")  # equals forward.count("C")
        if g_fwd > g_rev:
            return reverse
        if g_fwd < g_rev:
            return forward
        # Equal G count — prefer the strand that does not start with G
        if forward.startswith("G") and not reverse.startswith("G"):
            return reverse
        return forward

    def calculate_probe_tm(self, sequence: str, params: dict) -> float:
        """Tính Tm của probe.

        Args:
            sequence: Chuỗi probe.
            params: Dict tham số nhiệt động học.

        Returns:
            Tm (°C).
        """
        try:
            import primer3
        except ImportError as exc:
            raise ImportError("primer3-py is required for probe Tm calculation") from exc
        return primer3.calc_tm(
            sequence.upper(),
            mv_conc=float(params.get("mv_conc", 50.0)),
            dv_conc=float(params.get("dv_conc", 1.5)),
            dntp_conc=float(params.get("dntp_conc", 0.6)),
            dna_conc=float(params.get("dna_conc", 250.0)),
        )
