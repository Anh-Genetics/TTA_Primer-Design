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
"""

from __future__ import annotations

import logging
import re

import primer3

from tta_primer_design.config import AppConfig
from tta_primer_design.models import Oligo, PrimerPair, ProcessedSequence

logger = logging.getLogger("tta_primer_design.modules.probe_designer")

# Probe length range (bp)
_PROBE_MIN_LEN = 18
_PROBE_MAX_LEN = 30

# TaqMan probe Tm requirements
_PROBE_TM_MIN = 65.0  # °C
_PROBE_TM_MAX = 72.0  # °C
_PROBE_TM_DELTA = 8.0  # Probe Tm must exceed primer Tm by this much

# Hairpin ΔG threshold (kcal/mol) — avoid if more negative than this
_HAIRPIN_DG_THRESHOLD = -2000.0  # cal/mol (primer3 returns cal/mol, not kcal/mol)

# Minimum GC content for probe
_PROBE_GC_MIN = 40.0  # %
_PROBE_GC_MAX = 80.0  # %

# Max consecutive identical bases (poly-run)
_MAX_POLY_RUN = 4

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def _gc_percent(seq: str) -> float:
    """Calculate GC content percentage."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100.0


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
        seq = processed_seq.sequence
        for pair in primer_pairs:
            amplicon = pair.amplicon_sequence
            if not amplicon:
                left_start = pair.left_primer.start
                amplicon_size = pair.amplicon_size or (
                    pair.right_primer.start + len(pair.right_primer.sequence) - left_start
                )
                amplicon = seq[left_start : left_start + amplicon_size]
                if not amplicon:
                    logger.warning("Cannot extract amplicon for pair %s", pair.pair_id)
                    continue

            probe = self.design_taqman_probe(amplicon, pair)
            if probe is not None:
                pair.probe = probe
                logger.debug("Designed probe for pair %s: %s", pair.pair_id, probe.sequence)
            else:
                logger.warning("No valid probe found for pair %s", pair.pair_id)

        return primer_pairs

    def design_taqman_probe(
        self,
        amplicon: str,
        pair: PrimerPair,
    ) -> Oligo | None:
        """Thiết kế TaqMan probe cho một amplicon.

        Sliding window qua phần nội bộ của amplicon (không overlap với primer).
        Chọn probe có điểm tốt nhất dựa trên Tm và GC%.

        Args:
            amplicon: Chuỗi amplicon (bao gồm cả primer sequences).
            pair: PrimerPair tương ứng.

        Returns:
            Oligo probe hoặc None nếu không thiết kế được.
        """
        left_len = len(pair.left_primer.sequence)
        right_len = len(pair.right_primer.sequence)
        amplicon_len = len(amplicon)

        # Interior region: exclude primer binding sites
        interior_start = left_len
        interior_end = amplicon_len - right_len

        if interior_end - interior_start < _PROBE_MIN_LEN:
            logger.debug(
                "Amplicon interior too short for probe: %d bp", interior_end - interior_start
            )
            return None

        # Average primer Tm for threshold check
        avg_primer_tm = (pair.left_primer.tm + pair.right_primer.tm) / 2.0

        best_probe: Oligo | None = None
        best_score = -1.0

        params = {
            "mv_conc": 50.0,
            "dv_conc": 1.5,
            "dntp_conc": 0.25,
            "dna_conc": 250.0,
        }

        for probe_len in range(_PROBE_MIN_LEN, _PROBE_MAX_LEN + 1):
            for pos in range(interior_start, interior_end - probe_len + 1):
                candidate_seq = self.select_best_strand(amplicon, pos, probe_len)
                tm = self.calculate_probe_tm(candidate_seq, params)

                # Tm probe must exceed avg primer Tm by at least _PROBE_TM_DELTA
                if tm < avg_primer_tm + _PROBE_TM_DELTA:
                    continue
                if not (_PROBE_TM_MIN <= tm <= _PROBE_TM_MAX):
                    continue

                gc = _gc_percent(candidate_seq)
                if not (_PROBE_GC_MIN <= gc <= _PROBE_GC_MAX):
                    continue

                # Position of probe on template (amplicon-relative)
                probe_start_on_template = pair.left_primer.start + pos

                oligo = Oligo(
                    sequence=candidate_seq,
                    start=probe_start_on_template,
                    tm=tm,
                    gc_percent=gc,
                )

                passed, violations = self.check_probe_rules(oligo)
                if not passed:
                    logger.debug("Probe candidate failed rules: %s", violations)
                    continue

                # Score: prefer Tm in middle of range, high GC within range
                tm_score = 1.0 - abs(tm - 68.0) / 4.0  # 68°C is ideal
                gc_score = 1.0 - abs(gc - 55.0) / 25.0  # 55% GC is ideal
                score = tm_score * 0.6 + gc_score * 0.4

                if score > best_score:
                    best_score = score
                    best_probe = oligo

        return best_probe

    def check_probe_rules(self, probe: Oligo) -> tuple[bool, list[str]]:
        """Kiểm tra probe theo các quy tắc TaqMan.

        Kiểm tra:
            1. Không có G ở đầu 5'
            2. Không có polyG/polyC > 4 nucleotide liên tiếp
            3. ΔG hairpin không quá âm (> _HAIRPIN_DG_THRESHOLD)

        Args:
            probe: Oligo probe cần kiểm tra.

        Returns:
            (pass_flag, list_of_violations)
        """
        seq = probe.sequence.upper()
        violations: list[str] = []

        # Rule 1: No G at 5' end
        if seq.startswith("G"):
            violations.append("5' end is G (quenches fluorophore)")

        # Rule 2: No polyG or polyC > max_run
        pattern = rf"([GC])\1{{{_MAX_POLY_RUN},}}"
        if re.search(pattern, seq):
            violations.append(f"polyG/polyC run > {_MAX_POLY_RUN} bases")

        # Rule 3: Hairpin ΔG check (primer3 returns cal/mol)
        try:
            hairpin = primer3.calc_hairpin(seq)
            if hairpin.dg < _HAIRPIN_DG_THRESHOLD:
                violations.append(f"hairpin ΔG = {hairpin.dg / 1000:.2f} kcal/mol < threshold")
        except Exception as exc:
            logger.debug("Hairpin calculation failed: %s", exc)

        return len(violations) == 0, violations

    def select_best_strand(self, amplicon: str, pos: int, length: int) -> str:
        """Chọn strand tốt hơn (ít G hơn ở 5') cho probe.

        Nếu G > C trong probe → dùng reverse complement (quy tắc TaqMan số 4).

        Args:
            amplicon: Chuỗi amplicon.
            pos: Vị trí bắt đầu probe trên amplicon (0-based).
            length: Độ dài probe.

        Returns:
            Chuỗi probe trên strand tốt hơn.
        """
        fwd = amplicon[pos : pos + length].upper()
        rc = _reverse_complement(fwd)

        g_count_fwd = fwd.count("G")
        c_count_fwd = fwd.count("C")

        # Use RC if G > C in forward strand
        return rc if g_count_fwd > c_count_fwd else fwd

    def calculate_probe_tm(self, sequence: str, params: dict) -> float:
        """Tính Tm của probe dùng SantaLucia 1998 (primer3-py).

        Args:
            sequence: Chuỗi probe.
            params: Dict tham số nhiệt động học (mv_conc, dv_conc, dntp_conc, dna_conc).

        Returns:
            Tm (°C).
        """
        return primer3.calc_tm(
            sequence,
            mv_conc=params.get("mv_conc", 50.0),
            dv_conc=params.get("dv_conc", 1.5),
            dntp_conc=params.get("dntp_conc", 0.25),
            dna_conc=params.get("dna_conc", 250.0),
        )
