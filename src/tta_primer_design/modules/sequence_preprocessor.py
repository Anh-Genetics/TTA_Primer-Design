"""Module 3 — SequencePreprocessor: tiền xử lý sequence trước Primer3.

Chức năng:
    - Mask low-complexity regions (DUST algorithm đơn giản)
    - Validate sequence (chỉ chấp nhận ATCGN)
    - Tính GC content và complexity score
    - Xác định SEQUENCE_EXCLUDED_REGION (introns, repeats)
    - Xác định SEQUENCE_TARGET (vùng amplicon bắt buộc phủ)
    - Tính toán vị trí exon-exon junction
"""

from __future__ import annotations

import logging
import re

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget, ProcessedSequence

logger = logging.getLogger("tta_primer_design.modules.sequence_preprocessor")

_VALID_BASES = re.compile(r"^[ATCGNatcgn]+$")


class SequencePreprocessor:
    """Tiền xử lý sequence trước khi đưa vào Primer3.

    Args:
        config: AppConfig đã load.

    Example::

        preprocessor = SequencePreprocessor(config)
        proc_seq = preprocessor.process(seq_record, target)
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def process(self, seq_record: object, target: DesignTarget) -> ProcessedSequence:
        """Tiền xử lý sequence.

        Luồng xử lý:
            1. Lấy chuỗi str từ seq_record (hỗ trợ str, BioPython SeqRecord/Seq).
            2. Validate ký tự hợp lệ.
            3. Mask low-complexity regions.
            4. Tính GC content và complexity score.
            5. Xác định excluded_regions, target_regions, included_region từ target.
            6. Xác định exon-exon junction nếu target.exon_junction và seq_record có exon_coords.

        Args:
            seq_record: BioPython SeqRecord (hoặc chuỗi str thuần).
            target: DesignTarget tương ứng.

        Returns:
            ProcessedSequence sẵn sàng cho Primer3.

        Raises:
            ValueError: Nếu sequence chứa ký tự không hợp lệ.
        """
        # --- 1. Extract string sequence ---
        raw_seq = _extract_sequence_str(seq_record)

        # --- 2. Validate ---
        self.validate_sequence(raw_seq)

        # --- 3. Mask low-complexity ---
        masked_seq = self.mask_low_complexity(raw_seq)

        # --- 4. GC content & complexity score ---
        gc = self.calculate_gc_content(masked_seq)
        complexity = _calculate_complexity_score(masked_seq)

        # --- 5. Regions from target ---
        excluded: list[tuple[int, int]] = list(target.region_exclude or [])
        target_regions: list[tuple[int, int]] = []
        included_region: tuple[int, int] | None = None

        if target.region_include:
            start, end = target.region_include
            included_region = (start, end - start)  # convert to (start, length) for Primer3
            target_regions.append(target.region_include)

        # --- 6. Exon-exon junctions ---
        exon_junctions: list[int] = []
        if target.exon_junction:
            exon_junctions = _extract_exon_junctions(seq_record)

        logger.debug(
            "process: id=%s len=%d gc=%.2f complexity=%.2f excluded=%d junctions=%d",
            target.target_id,
            len(masked_seq),
            gc,
            complexity,
            len(excluded),
            len(exon_junctions),
        )

        return ProcessedSequence(
            sequence=masked_seq,
            excluded_regions=excluded,
            target_regions=target_regions,
            included_region=included_region,
            exon_junctions=exon_junctions,
            gc_content=gc,
            complexity_score=complexity,
        )

    def validate_sequence(self, sequence: str) -> bool:
        """Kiểm tra sequence chỉ chứa các base hợp lệ (ATCGN).

        Args:
            sequence: Chuỗi nucleotide cần kiểm tra.

        Returns:
            True nếu hợp lệ.

        Raises:
            ValueError: Nếu sequence chứa ký tự không hợp lệ.
        """
        if not sequence:
            raise ValueError("Sequence is empty")
        if not _VALID_BASES.match(sequence):
            invalid = set(sequence.upper()) - set("ATCGN")
            raise ValueError(f"Sequence contains invalid characters: {sorted(invalid)}")
        return True

    def calculate_gc_content(self, sequence: str) -> float:
        """Tính tỷ lệ GC content.

        Args:
            sequence: Chuỗi nucleotide.

        Returns:
            Tỷ lệ GC (0.0–1.0).
        """
        seq = sequence.upper()
        if not seq:
            return 0.0
        gc = seq.count("G") + seq.count("C")
        return gc / len(seq)

    def mask_low_complexity(self, sequence: str, window: int = 12, threshold: float = 0.7) -> str:
        """Mask các vùng low-complexity bằng 'N'.

        Sử dụng thuật toán đơn giản: sliding window, nếu 1 base chiếm
        > threshold trong cửa sổ → toàn bộ cửa sổ bị đánh dấu để mask.
        Các vị trí bị đánh dấu được thay bằng 'N'.

        Args:
            sequence: Chuỗi nucleotide gốc.
            window: Kích thước cửa sổ trượt.
            threshold: Ngưỡng tỷ lệ base thống trị (0.0–1.0).

        Returns:
            Chuỗi đã mask (cùng độ dài với input).
        """
        seq = sequence.upper()
        n = len(seq)
        if n == 0:
            return sequence
        if window > n:
            window = n

        to_mask = bytearray(n)  # 0 = keep, 1 = mask

        for i in range(n - window + 1):
            win = seq[i : i + window]
            for base in "ATCG":
                if win.count(base) / window > threshold:
                    for j in range(i, i + window):
                        to_mask[j] = 1
                    break

        result = list(seq)
        for i, flag in enumerate(to_mask):
            if flag:
                result[i] = "N"
        return "".join(result)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _extract_sequence_str(seq_record: object) -> str:
    """Lấy chuỗi str từ SeqRecord, Seq, hoặc str thuần."""
    if isinstance(seq_record, str):
        return seq_record.upper()
    # BioPython SeqRecord
    if hasattr(seq_record, "seq"):
        return str(seq_record.seq).upper()
    # BioPython Seq
    if hasattr(seq_record, "__str__") and not isinstance(seq_record, (int, float, bytes)):
        return str(seq_record).upper()
    raise TypeError(f"Unsupported seq_record type: {type(seq_record)}")


def _calculate_complexity_score(sequence: str) -> float:
    """Tính complexity score dựa trên entropy Shannon đơn giản (0.0–1.0).

    Score = 1.0 nghĩa là tối đa entropy (bộ 4 base đều nhau).
    Score = 0.0 nghĩa là chỉ có 1 base (low complexity).
    """
    seq = sequence.upper().replace("N", "")
    if not seq:
        return 0.0
    import math

    n = len(seq)
    counts = {b: seq.count(b) for b in "ATCG"}
    entropy = 0.0
    for c in counts.values():
        if c > 0:
            p = c / n
            entropy -= p * math.log2(p)
    # Normalize by max entropy (log2(4) = 2.0)
    return entropy / 2.0


def _extract_exon_junctions(seq_record: object) -> list[int]:
    """Trích xuất vị trí exon-exon junction từ SeqRecord (nếu có).

    Tìm trong SeqRecord.features (type="exon") hoặc attribute exon_coords.
    Trả về danh sách vị trí junction (vị trí kết thúc exon, 0-based).
    """
    junctions: list[int] = []

    # Nếu seq_record có exon_coords attribute (từ SequenceFetcher)
    if hasattr(seq_record, "exon_coords"):
        coords: list[tuple[int, int]] = seq_record.exon_coords
        # Junction = điểm cuối mỗi exon (trừ exon cuối)
        for _start, end in coords[:-1]:
            junctions.append(end)
        return junctions

    # Thử parse từ BioPython SeqRecord.features
    if hasattr(seq_record, "features"):
        exon_ends = []
        for feat in seq_record.features:
            if feat.type == "exon":
                exon_ends.append(int(feat.location.end))
        exon_ends.sort()
        junctions = exon_ends[:-1] if len(exon_ends) > 1 else []

    return junctions
