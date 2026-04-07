"""Module 3 — SequencePreprocessor: tiền xử lý sequence trước Primer3.

Chức năng:
    - Mask low-complexity regions (DUST algorithm đơn giản)
    - Validate sequence (chỉ chấp nhận ATCGN)
    - Tính GC content và complexity score
    - Xác định SEQUENCE_EXCLUDED_REGION (introns, repeats)
    - Xác định SEQUENCE_TARGET (vùng amplicon bắt buộc phủ)
    - Tính toán vị trí exon-exon junction

TODO (Sprint 2):
    - Implement DUST mask
    - Implement repeat masking (simple k-mer based)
    - Integrate với exon coords từ SequenceFetcher
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

        Args:
            seq_record: BioPython SeqRecord (hoặc chuỗi str thuần).
            target: DesignTarget tương ứng.

        Returns:
            ProcessedSequence sẵn sàng cho Primer3.

        Raises:
            NotImplementedError: Chưa implement đầy đủ (Sprint 2).
        """
        raise NotImplementedError("process chưa được implement (Sprint 2)")

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
        > threshold → mask thành 'N'.

        Args:
            sequence: Chuỗi nucleotide gốc.
            window: Kích thước cửa sổ trượt.
            threshold: Ngưỡng tỷ lệ base thống trị.

        Returns:
            Chuỗi đã mask.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("mask_low_complexity chưa được implement (Sprint 2)")
