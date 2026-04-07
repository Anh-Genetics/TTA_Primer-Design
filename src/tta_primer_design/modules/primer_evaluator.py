"""Module — PrimerEvaluator: kiểm thử và đánh giá một cặp mồi sẵn có.

Chức năng:
    1. Tính nhiệt động học (Tm, GC%, GC clamp, hairpin, self-dimer, hetero-dimer)
       bằng primer3-py.
    2. Kiểm tra tính đặc hiệu (specificity) qua NCBI BLAST API (Biopython NCBIWWW).
    3. Dự đoán amplicon off-target từ kết quả BLAST.
    4. Tổng hợp ``EvaluationReport`` với khuyến nghị PASS / WARNING / FAIL.

Sử dụng:
    evaluator = PrimerEvaluator(config)
    report = evaluator.evaluate(
        pair_id="my_primers",
        left_seq="GCACTGACCTCCCACTTCAA",
        right_seq="TTGCTGATCCACATCTGCTG",
        organism="Homo sapiens",
        run_blast=True,
    )
    print(report.overall_recommendation)

Ngưỡng mặc định:
    Primer  Tm:     55–65 °C
    Probe   Tm:     65–72 °C
    GC%:            40–60 %
    Hairpin ΔG:     > -9 kcal/mol
    Self-dimer ΔG:  > -9 kcal/mol
    GC clamp 3':    1–3 G/C trong 5 bp cuối
    Repeat run:     ≤ 4 base giống nhau liên tiếp
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from io import StringIO
from typing import Any

from tta_primer_design.config import AppConfig
from tta_primer_design.modules.blast_specificity import (
    BlastHit,
    OffTargetAmplicon,
    SpecificityResult,
)

logger = logging.getLogger("tta_primer_design.modules.primer_evaluator")

# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------

_TM_PRIMER_MIN: float = 55.0
_TM_PRIMER_MAX: float = 65.0
_TM_PROBE_MIN: float = 65.0
_TM_PROBE_MAX: float = 72.0
_GC_MIN: float = 40.0
_GC_MAX: float = 60.0
_HAIRPIN_DG_THRESHOLD: float = -9.0  # kcal/mol — phải > ngưỡng này
_DIMER_DG_THRESHOLD: float = -9.0  # kcal/mol — phải > ngưỡng này
_GC_CLAMP_WINDOW: int = 5
_GC_CLAMP_MIN: int = 1
_GC_CLAMP_MAX: int = 3
_MAX_REPEAT_RUN: int = 4

# BLAST defaults
_BLAST_WORD_SIZE: int = 7
_BLAST_EVALUE: float = 1000.0
_BLAST_HITLIST_SIZE: int = 50
_BLAST_RETRY_DELAY: float = 5.0  # giây


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass
class ThermoProfile:
    """Profile nhiệt động học cho một oligo đơn.

    Attributes:
        sequence: Chuỗi oligo (đã uppercase).
        tm: Nhiệt độ nóng chảy (°C).
        gc_percent: Tỷ lệ GC (%).
        hairpin_dg: ΔG cấu trúc hairpin (kcal/mol); âm = có hairpin.
        homodimer_dg: ΔG self-dimer (kcal/mol).
        gc_clamp_ok: True nếu GC clamp ở 3' end hợp lệ.
        repeat_ok: True nếu không có repeat run > _MAX_REPEAT_RUN.
        pass_all: True nếu tất cả thông số trong ngưỡng cho phép.
        warnings: Danh sách cảnh báo cụ thể.
    """

    sequence: str
    tm: float = 0.0
    gc_percent: float = 0.0
    hairpin_dg: float = 0.0
    homodimer_dg: float = 0.0
    gc_clamp_ok: bool = True
    repeat_ok: bool = True
    pass_all: bool = True
    warnings: list[str] = field(default_factory=list)


@dataclass
class EvaluationReport:
    """Kết quả đánh giá toàn diện một cặp primer sẵn có.

    Attributes:
        pair_id: Tên/ID cặp primer.
        left_sequence: Trình tự primer xuôi (5'→3').
        right_sequence: Trình tự primer ngược (5'→3').
        probe_sequence: Trình tự probe (None nếu không có).
        left_thermo: Profile nhiệt động học của primer trái.
        right_thermo: Profile nhiệt động học của primer phải.
        probe_thermo: Profile nhiệt động học của probe (None nếu không có probe).
        heterodimer_dg: ΔG hetero-dimer giữa left và right (kcal/mol).
        specificity: Kết quả BLAST specificity (None nếu chưa chạy BLAST).
        blast_performed: True nếu đã thực hiện BLAST.
        pass_thermodynamics: True nếu cả 2 primer (và probe nếu có) đạt thermo.
        overall_recommendation: "PASS" | "WARNING" | "FAIL".
        summary_warnings: Danh sách tóm tắt các vấn đề phát hiện.
    """

    pair_id: str
    left_sequence: str
    right_sequence: str
    probe_sequence: str | None
    left_thermo: ThermoProfile
    right_thermo: ThermoProfile
    probe_thermo: ThermoProfile | None
    heterodimer_dg: float
    specificity: SpecificityResult | None = None
    blast_performed: bool = False
    pass_thermodynamics: bool = True
    overall_recommendation: str = "PASS"
    summary_warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------


class PrimerEvaluator:
    """Kiểm thử và đánh giá toàn diện một cặp mồi/probe sẵn có.

    Sử dụng primer3-py để tính nhiệt động học và Biopython NCBIWWW để chạy
    NCBI BLAST API kiểm tra specificity.

    Args:
        config: AppConfig đã load.
        max_amplicon_size: Kích thước off-target amplicon tối đa cần flag (bp).

    Example::

        evaluator = PrimerEvaluator(config)
        report = evaluator.evaluate(
            pair_id="ACTB_qPCR",
            left_seq="GCACTGACCTCCCACTTCAA",
            right_seq="TTGCTGATCCACATCTGCTG",
            organism="Homo sapiens",
            run_blast=True,
        )
        print(report.overall_recommendation)
    """

    def __init__(self, config: AppConfig, max_amplicon_size: int = 4000) -> None:
        self.config = config
        self.max_amplicon_size = max_amplicon_size

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def evaluate(
        self,
        pair_id: str,
        left_seq: str,
        right_seq: str,
        probe_seq: str | None = None,
        organism: str = "Homo sapiens",
        run_blast: bool = True,
    ) -> EvaluationReport:
        """Đánh giá toàn diện một cặp primer.

        Luồng xử lý:
            1. Chuẩn hoá chuỗi (uppercase, strip).
            2. Tính nhiệt động học cho từng oligo bằng primer3-py.
            3. Tính hetero-dimer ΔG giữa left và right primer.
            4. (Tuỳ chọn) Chạy NCBI BLAST API để kiểm tra specificity.
            5. Tổng hợp EvaluationReport với khuyến nghị cuối cùng.

        Args:
            pair_id: Tên/ID cặp primer.
            left_seq: Trình tự primer xuôi (5'→3').
            right_seq: Trình tự primer ngược (5'→3').
            probe_seq: Trình tự probe tuỳ chọn.
            organism: Tên loài NCBI (dùng cho BLAST entrez_query).
            run_blast: True để chạy NCBI BLAST (cần internet, có thể chậm).

        Returns:
            EvaluationReport đầy đủ.

        Raises:
            ValueError: Nếu chuỗi primer rỗng hoặc chứa ký tự không hợp lệ.
            ImportError: Nếu primer3-py chưa được cài.
        """
        left_seq = self._validate_sequence(left_seq, "left_primer")
        right_seq = self._validate_sequence(right_seq, "right_primer")
        if probe_seq is not None:
            probe_seq = self._validate_sequence(probe_seq, "probe")

        logger.info(
            "Evaluating primer pair '%s': L=%s | R=%s",
            pair_id,
            left_seq,
            right_seq,
        )

        # 1. Thermodynamics
        left_thermo = self._calc_thermodynamics(left_seq, is_probe=False)
        right_thermo = self._calc_thermodynamics(right_seq, is_probe=False)
        probe_thermo = (
            self._calc_thermodynamics(probe_seq, is_probe=True) if probe_seq else None
        )
        heterodimer_dg = self._calc_heterodimer_dg(left_seq, right_seq)

        # 2. BLAST specificity
        specificity: SpecificityResult | None = None
        blast_performed = False
        if run_blast:
            try:
                specificity = self._blast_pair(
                    pair_id=pair_id,
                    left_seq=left_seq,
                    right_seq=right_seq,
                    probe_seq=probe_seq,
                    organism=organism,
                )
                blast_performed = True
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "BLAST failed for pair '%s': %s — continuing without specificity data.",
                    pair_id,
                    exc,
                )

        # 3. Build report
        return self._build_report(
            pair_id=pair_id,
            left_seq=left_seq,
            right_seq=right_seq,
            probe_seq=probe_seq,
            left_thermo=left_thermo,
            right_thermo=right_thermo,
            probe_thermo=probe_thermo,
            heterodimer_dg=heterodimer_dg,
            specificity=specificity,
            blast_performed=blast_performed,
        )

    # ------------------------------------------------------------------
    # Thermodynamics helpers
    # ------------------------------------------------------------------

    def _calc_thermodynamics(self, sequence: str, is_probe: bool = False) -> ThermoProfile:
        """Tính toán và đánh giá profile nhiệt động học cho một oligo.

        Args:
            sequence: Chuỗi oligo (uppercase, IUPAC DNA).
            is_probe: True nếu là probe (dùng ngưỡng Tm khác).

        Returns:
            ThermoProfile đầy đủ.

        Raises:
            ImportError: Nếu primer3-py chưa được cài.
        """
        try:
            import primer3
        except ImportError as exc:
            raise ImportError(
                "primer3-py chưa được cài. Chạy: pip install primer3-py"
            ) from exc

        tm = float(primer3.calc_tm(sequence))
        gc_percent = self._calc_gc_percent(sequence)

        hairpin_result = primer3.calc_hairpin(sequence)
        hairpin_dg = hairpin_result.dg / 1000.0  # cal/mol → kcal/mol

        homodimer_result = primer3.calc_homodimer(sequence)
        homodimer_dg = homodimer_result.dg / 1000.0  # cal/mol → kcal/mol

        gc_clamp_ok = self._check_gc_clamp(sequence)
        repeat_ok = self._check_repeat_runs(sequence)

        # Evaluate against thresholds
        warnings: list[str] = []
        pass_all = True

        tm_min = _TM_PROBE_MIN if is_probe else _TM_PRIMER_MIN
        tm_max = _TM_PROBE_MAX if is_probe else _TM_PRIMER_MAX
        label = "probe" if is_probe else "primer"

        if not (tm_min <= tm <= tm_max):
            warnings.append(
                f"Tm={tm:.1f}°C ngoài khoảng [{tm_min}–{tm_max}°C] cho {label}"
            )
            pass_all = False

        if not (_GC_MIN <= gc_percent <= _GC_MAX):
            warnings.append(
                f"GC%={gc_percent:.1f}% ngoài khoảng [{_GC_MIN}–{_GC_MAX}%]"
            )
            pass_all = False

        if hairpin_dg < _HAIRPIN_DG_THRESHOLD:
            warnings.append(
                f"Hairpin ΔG={hairpin_dg:.2f} kcal/mol < {_HAIRPIN_DG_THRESHOLD} (có cấu trúc hairpin)"
            )
            pass_all = False

        if homodimer_dg < _DIMER_DG_THRESHOLD:
            warnings.append(
                f"Self-dimer ΔG={homodimer_dg:.2f} kcal/mol < {_DIMER_DG_THRESHOLD}"
            )
            pass_all = False

        if not gc_clamp_ok:
            warnings.append(
                f"GC clamp không hợp lệ: cần {_GC_CLAMP_MIN}–{_GC_CLAMP_MAX} G/C "
                f"trong {_GC_CLAMP_WINDOW} bp cuối 3'"
            )
            pass_all = False

        if not repeat_ok:
            warnings.append(f"Có ≥{_MAX_REPEAT_RUN + 1} base giống nhau liên tiếp")
            pass_all = False

        return ThermoProfile(
            sequence=sequence,
            tm=tm,
            gc_percent=gc_percent,
            hairpin_dg=hairpin_dg,
            homodimer_dg=homodimer_dg,
            gc_clamp_ok=gc_clamp_ok,
            repeat_ok=repeat_ok,
            pass_all=pass_all,
            warnings=warnings,
        )

    def _calc_heterodimer_dg(self, seq1: str, seq2: str) -> float:
        """Tính ΔG hetero-dimer giữa 2 oligo (kcal/mol).

        Args:
            seq1: Chuỗi oligo thứ nhất.
            seq2: Chuỗi oligo thứ hai.

        Returns:
            ΔG hetero-dimer (kcal/mol).

        Raises:
            ImportError: Nếu primer3-py chưa được cài.
        """
        try:
            import primer3
        except ImportError as exc:
            raise ImportError(
                "primer3-py chưa được cài. Chạy: pip install primer3-py"
            ) from exc

        result = primer3.calc_heterodimer(seq1, seq2)
        return result.dg / 1000.0  # cal/mol → kcal/mol

    def _calc_gc_percent(self, sequence: str) -> float:
        """Tính tỷ lệ GC (%).

        Args:
            sequence: Chuỗi nucleotide (uppercase).

        Returns:
            GC% (0.0–100.0).
        """
        if not sequence:
            return 0.0
        gc = sum(1 for base in sequence if base in "GC")
        return gc / len(sequence) * 100.0

    def _check_gc_clamp(
        self,
        sequence: str,
        window: int = _GC_CLAMP_WINDOW,
        min_gc: int = _GC_CLAMP_MIN,
        max_gc: int = _GC_CLAMP_MAX,
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
        tail = sequence[-window:]
        gc_count = sum(1 for b in tail if b in "GC")
        return min_gc <= gc_count <= max_gc

    def _check_repeat_runs(self, sequence: str, max_run: int = _MAX_REPEAT_RUN) -> bool:
        """Kiểm tra không có repeat run quá dài.

        Args:
            sequence: Chuỗi oligo.
            max_run: Số base giống nhau tối đa cho phép liên tiếp.

        Returns:
            True nếu không có run nào vượt max_run.
        """
        if not sequence:
            return True
        # Tìm run của bất kỳ base nào: A+, T+, C+, G+
        pattern = re.compile(r"([ATCG])\1{" + str(max_run) + r",}", re.IGNORECASE)
        return pattern.search(sequence) is None

    # ------------------------------------------------------------------
    # BLAST helpers
    # ------------------------------------------------------------------

    def _blast_pair(
        self,
        pair_id: str,
        left_seq: str,
        right_seq: str,
        probe_seq: str | None,
        organism: str,
    ) -> SpecificityResult:
        """BLAST cả hai primer (và probe nếu có) qua NCBI BLAST API.

        Args:
            pair_id: ID cặp primer.
            left_seq: Primer trái.
            right_seq: Primer phải.
            probe_seq: Probe (tuỳ chọn).
            organism: Tên loài cho entrez_query.

        Returns:
            SpecificityResult với danh sách hits và off-target amplicons.
        """
        logger.info("Running NCBI BLAST for pair '%s' (organism: %s)…", pair_id, organism)

        left_hits = self._blast_oligo(left_seq, organism)
        logger.debug("Left primer: %d BLAST hits", len(left_hits))

        right_hits = self._blast_oligo(right_seq, organism)
        logger.debug("Right primer: %d BLAST hits", len(right_hits))

        probe_hits: list[BlastHit] | None = None
        if probe_seq:
            probe_hits = self._blast_oligo(probe_seq, organism)
            logger.debug("Probe: %d BLAST hits", len(probe_hits))

        off_targets = self._predict_offtarget_amplicons(left_hits, right_hits)

        # Specificity score: 100 nếu không có off-target, giảm dần
        if off_targets:
            score = max(0.0, 100.0 - len(off_targets) * 20.0)
            is_specific = False
        else:
            score = 100.0
            is_specific = True

        return SpecificityResult(
            primer_pair_id=pair_id,
            is_specific=is_specific,
            off_target_amplicons=off_targets,
            specificity_score=score,
            blast_hits_left=left_hits,
            blast_hits_right=right_hits,
            blast_hits_probe=probe_hits,
        )

    def _blast_oligo(self, sequence: str, organism: str) -> list[BlastHit]:
        """BLAST một oligo qua NCBI BLAST Web API (blastn/nt).

        Dùng Biopython ``NCBIWWW.qblast`` với entrez_query lọc theo loài.
        Tự động retry 1 lần nếu gặp lỗi mạng thoáng qua.

        Args:
            sequence: Chuỗi oligo.
            organism: Tên loài NCBI.

        Returns:
            Danh sách BlastHit.

        Raises:
            ImportError: Nếu Biopython chưa được cài.
            RuntimeError: Nếu BLAST request thất bại sau retry.
        """
        try:
            from Bio.Blast import NCBIWWW
        except ImportError as exc:
            raise ImportError(
                "Biopython chưa được cài. Chạy: pip install biopython"
            ) from exc

        entrez_query = f"({organism}[Organism])"
        database = self.config.blast.database or "nt"

        for attempt in (1, 2):
            try:
                handle = NCBIWWW.qblast(
                    "blastn",
                    database,
                    sequence,
                    entrez_query=entrez_query,
                    word_size=_BLAST_WORD_SIZE,
                    expect=_BLAST_EVALUE,
                    hitlist_size=_BLAST_HITLIST_SIZE,
                    format_type="XML",
                )
                xml_data = handle.read()
                handle.close()
                break
            except Exception as exc:
                if attempt == 2:
                    raise RuntimeError(
                        f"NCBI BLAST request thất bại sau 2 lần thử: {exc}"
                    ) from exc
                logger.warning(
                    "BLAST attempt %d failed (%s) — retrying in %ss…",
                    attempt,
                    exc,
                    _BLAST_RETRY_DELAY,
                )
                time.sleep(_BLAST_RETRY_DELAY)

        return self._parse_blast_xml(xml_data)

    def _parse_blast_xml(self, xml_data: str | bytes) -> list[BlastHit]:
        """Parse kết quả BLAST XML thành danh sách BlastHit.

        Args:
            xml_data: Nội dung XML trả về từ NCBI BLAST.

        Returns:
            Danh sách BlastHit (mỗi HSP tốt nhất của mỗi hit).
        """
        try:
            from Bio.Blast import NCBIXML
        except ImportError as exc:
            raise ImportError("Biopython chưa được cài.") from exc

        if isinstance(xml_data, bytes):
            xml_str = xml_data.decode("utf-8", errors="replace")
        else:
            xml_str = xml_data

        hits: list[BlastHit] = []
        try:
            record = NCBIXML.read(StringIO(xml_str))
        except Exception as exc:
            logger.error("Không thể parse BLAST XML: %s", exc)
            return hits

        for alignment in record.alignments:
            if not alignment.hsps:
                continue
            # Lấy HSP tốt nhất (bit_score cao nhất)
            best_hsp = max(alignment.hsps, key=lambda h: h.bits)

            mismatches_3prime = self._count_3prime_mismatches(best_hsp)

            hits.append(
                BlastHit(
                    subject_id=alignment.accession,
                    subject_title=alignment.title[:200],
                    identity=best_hsp.identities / best_hsp.align_length * 100.0,
                    alignment_length=best_hsp.align_length,
                    mismatches=best_hsp.align_length - best_hsp.identities,
                    gaps=best_hsp.gaps,
                    query_start=best_hsp.query_start,
                    query_end=best_hsp.query_end,
                    subject_start=best_hsp.sbjct_start,
                    subject_end=best_hsp.sbjct_end,
                    evalue=best_hsp.expect,
                    bit_score=best_hsp.bits,
                    mismatches_3prime=mismatches_3prime,
                )
            )

        logger.debug("Parsed %d BLAST hits", len(hits))
        return hits

    def _count_3prime_mismatches(self, hsp: Any, window: int = 3) -> int:
        """Đếm số mismatch trong 3 bp cuối của query (3' end primer).

        Args:
            hsp: Biopython HSP object (có hsp.query và hsp.sbjct).
            window: Số bp ở 3' cần kiểm tra.

        Returns:
            Số mismatch ở 3' end.
        """
        try:
            query_tail = hsp.query[-window:].upper()
            sbjct_tail = hsp.sbjct[-window:].upper()
            mismatches = sum(
                1
                for q, s in zip(query_tail, sbjct_tail, strict=False)
                if q != s and q != "-" and s != "-"
            )
            return mismatches
        except Exception:
            return 0

    def _predict_offtarget_amplicons(
        self,
        left_hits: list[BlastHit],
        right_hits: list[BlastHit],
    ) -> list[OffTargetAmplicon]:
        """Dự đoán các off-target amplicon từ BLAST hits của left và right primer.

        Logic:
            - Left primer: tìm hit trên chuỗi xuôi (sbjct_start < sbjct_end).
            - Right primer: tìm hit trên chuỗi ngược (sbjct_start > sbjct_end).
            - Nếu cùng subject_id và khoảng cách < max_amplicon_size → off-target.

        Args:
            left_hits: BLAST hits của primer trái.
            right_hits: BLAST hits của primer phải.

        Returns:
            Danh sách OffTargetAmplicon tiềm năng.
        """
        # Index right hits theo subject_id: chỉ giữ minus-strand hits
        right_minus: dict[str, list[BlastHit]] = {}
        for h in right_hits:
            if h.subject_start > h.subject_end:  # minus strand
                right_minus.setdefault(h.subject_id, []).append(h)

        off_targets: list[OffTargetAmplicon] = []
        for lh in left_hits:
            if lh.subject_start >= lh.subject_end:
                # Bỏ qua minus strand hits của left primer
                continue
            rh_list = right_minus.get(lh.subject_id, [])
            for rh in rh_list:
                # rh.subject_start > rh.subject_end (minus strand)
                # vị trí 3' của right primer trên + strand ≈ rh.subject_end
                # vị trí 5' của right primer trên + strand ≈ rh.subject_start
                right_pos_start = rh.subject_end  # nhỏ hơn (5' end trên + strand)
                left_pos_end = lh.subject_end

                if right_pos_start <= left_pos_end:
                    continue  # overlap — không tạo amplicon

                amplicon_size = right_pos_start - lh.subject_start + 1
                if amplicon_size <= self.max_amplicon_size:
                    logger.debug(
                        "Off-target amplicon on %s: size=%d (left=%d–%d, right=%d–%d)",
                        lh.subject_id,
                        amplicon_size,
                        lh.subject_start,
                        lh.subject_end,
                        rh.subject_end,
                        rh.subject_start,
                    )
                    off_targets.append(
                        OffTargetAmplicon(
                            subject_id=lh.subject_id,
                            amplicon_size=amplicon_size,
                            left_hit=lh,
                            right_hit=rh,
                        )
                    )

        return off_targets

    # ------------------------------------------------------------------
    # Report builder
    # ------------------------------------------------------------------

    def _build_report(
        self,
        pair_id: str,
        left_seq: str,
        right_seq: str,
        probe_seq: str | None,
        left_thermo: ThermoProfile,
        right_thermo: ThermoProfile,
        probe_thermo: ThermoProfile | None,
        heterodimer_dg: float,
        specificity: SpecificityResult | None,
        blast_performed: bool,
    ) -> EvaluationReport:
        """Tổng hợp EvaluationReport và xác định khuyến nghị cuối cùng.

        Args:
            Các tham số đầu vào tương ứng với từng thành phần kết quả.

        Returns:
            EvaluationReport hoàn chỉnh với overall_recommendation.
        """
        warnings: list[str] = []
        is_fail = False

        # --- Thermodynamics ---
        for prefix, profile in [("Left primer", left_thermo), ("Right primer", right_thermo)]:
            for w in profile.warnings:
                warnings.append(f"{prefix}: {w}")
            if not profile.pass_all:
                is_fail = True

        if probe_thermo:
            for w in probe_thermo.warnings:
                warnings.append(f"Probe: {w}")
            if not probe_thermo.pass_all:
                is_fail = True

        # Hetero-dimer check
        if heterodimer_dg < _DIMER_DG_THRESHOLD:
            warnings.append(
                f"Hetero-dimer (left×right) ΔG={heterodimer_dg:.2f} kcal/mol "
                f"< {_DIMER_DG_THRESHOLD} — nguy cơ primer dimer"
            )
            is_fail = True

        pass_thermodynamics = not is_fail

        # --- Specificity ---
        if blast_performed and specificity is not None:
            if specificity.off_target_amplicons:
                n = len(specificity.off_target_amplicons)
                warnings.append(
                    f"BLAST: phát hiện {n} off-target amplicon(s) tiềm năng — "
                    f"specificity score={specificity.specificity_score:.0f}/100"
                )
                is_fail = True
            else:
                logger.info(
                    "BLAST: không có off-target amplicon — specificity score=%.0f/100",
                    specificity.specificity_score,
                )
        elif not blast_performed:
            # Thông tin: BLAST bị bỏ qua theo yêu cầu — không ảnh hưởng đến khuyến nghị
            logger.info("BLAST không được thực hiện — không có dữ liệu specificity.")

        # --- Overall recommendation ---
        if is_fail:
            recommendation = "FAIL"
        elif warnings:
            recommendation = "WARNING"
        else:
            recommendation = "PASS"

        return EvaluationReport(
            pair_id=pair_id,
            left_sequence=left_seq,
            right_sequence=right_seq,
            probe_sequence=probe_seq,
            left_thermo=left_thermo,
            right_thermo=right_thermo,
            probe_thermo=probe_thermo,
            heterodimer_dg=heterodimer_dg,
            specificity=specificity,
            blast_performed=blast_performed,
            pass_thermodynamics=pass_thermodynamics,
            overall_recommendation=recommendation,
            summary_warnings=warnings,
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    @staticmethod
    def _validate_sequence(seq: str, name: str) -> str:
        """Chuẩn hoá và validate chuỗi nucleotide.

        Args:
            seq: Chuỗi đầu vào.
            name: Tên trường để hiển thị trong thông báo lỗi.

        Returns:
            Chuỗi đã uppercase và strip whitespace.

        Raises:
            ValueError: Nếu chuỗi rỗng hoặc chứa ký tự không hợp lệ.
        """
        seq = seq.strip().upper()
        if not seq:
            raise ValueError(f"{name}: chuỗi không được rỗng")
        invalid = set(seq) - set("ATCGRYSWKMBDHVN")
        if invalid:
            raise ValueError(
                f"{name}: chứa ký tự không hợp lệ: {', '.join(sorted(invalid))}"
            )
        return seq
