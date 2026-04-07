"""Module 6 — BlastSpecificity: kiểm tra tính đặc hiệu của primer/probe.

Ba phương án:
    A) Dùng kết quả sẵn từ NCBI Primer-BLAST (ưu tiên)
    B) Chạy BLAST+ local cho batch lớn
    C) Gọi NCBI BLAST API cho từng oligo

Tiêu chí đánh giá:
    - Số off-target hits có thể tạo amplicon
    - Max amplicon size từ off-target (nếu > cutoff → FAIL)
    - Perfect match (mismatch = 0) trên non-target sequences
    - 3' end mismatch (≥1 mismatch ở 3 bp cuối = acceptable)

TODO (Sprint 3):
    - Implement NCBI BLAST API wrapper
    - Implement local BLAST+ subprocess runner
    - Implement off-target amplicon prediction
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from tta_primer_design.config import AppConfig
from tta_primer_design.models import PrimerPair

logger = logging.getLogger("tta_primer_design.modules.blast_specificity")


@dataclass
class BlastHit:
    """Một BLAST hit cho một oligo.

    Attributes:
        subject_id: ID sequence đích trong database.
        subject_title: Tên/mô tả sequence đích.
        identity: % identity.
        alignment_length: Độ dài alignment.
        mismatches: Số mismatch.
        gaps: Số gap.
        query_start: Vị trí bắt đầu trên query.
        query_end: Vị trí kết thúc trên query.
        subject_start: Vị trí bắt đầu trên subject.
        subject_end: Vị trí kết thúc trên subject.
        evalue: E-value.
        bit_score: Bit score.
        mismatches_3prime: Số mismatch ở 3 bp cuối primer.
    """

    subject_id: str
    subject_title: str = ""
    identity: float = 0.0
    alignment_length: int = 0
    mismatches: int = 0
    gaps: int = 0
    query_start: int = 0
    query_end: int = 0
    subject_start: int = 0
    subject_end: int = 0
    evalue: float = 0.0
    bit_score: float = 0.0
    mismatches_3prime: int = 0


@dataclass
class OffTargetAmplicon:
    """Một amplicon off-target tiềm năng (2 primer hit gần nhau).

    Attributes:
        subject_id: ID sequence đích.
        amplicon_size: Kích thước amplicon ước tính.
        left_hit: BLAST hit của primer trái.
        right_hit: BLAST hit của primer phải.
    """

    subject_id: str
    amplicon_size: int
    left_hit: BlastHit
    right_hit: BlastHit


@dataclass
class SpecificityResult:
    """Kết quả kiểm tra specificity cho một PrimerPair.

    Attributes:
        primer_pair_id: ID cặp primer.
        is_specific: True nếu pass specificity check.
        off_target_amplicons: Danh sách amplicon off-target.
        specificity_score: Điểm specificity (0–100).
        blast_hits_left: BLAST hits của primer trái.
        blast_hits_right: BLAST hits của primer phải.
        blast_hits_probe: BLAST hits của probe (nếu có).
    """

    primer_pair_id: str
    is_specific: bool = True
    off_target_amplicons: list[OffTargetAmplicon] = field(default_factory=list)
    specificity_score: float = 100.0
    blast_hits_left: list[BlastHit] = field(default_factory=list)
    blast_hits_right: list[BlastHit] = field(default_factory=list)
    blast_hits_probe: list[BlastHit] | None = None


class BlastSpecificity:
    """Kiểm tra specificity của primer/probe bằng BLAST.

    Args:
        config: AppConfig đã load.

    Example::

        checker = BlastSpecificity(config)
        results = checker.check_all(primer_pairs, organism="Homo sapiens")
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def check_all(
        self,
        primer_pairs: list[PrimerPair],
        organism: str = "Homo sapiens",
    ) -> list[PrimerPair]:
        """Kiểm tra specificity cho tất cả primer pairs.

        Gán ``specificity_result`` vào từng PrimerPair.

        Args:
            primer_pairs: Danh sách PrimerPair cần kiểm tra.
            organism: Tên loài.

        Returns:
            Danh sách PrimerPair đã có specificity_result.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("check_all chưa được implement (Sprint 3)")

    def check_pair(
        self,
        pair: PrimerPair,
        organism: str = "Homo sapiens",
    ) -> SpecificityResult:
        """Kiểm tra specificity cho một PrimerPair.

        Args:
            pair: PrimerPair cần kiểm tra.
            organism: Tên loài.

        Returns:
            SpecificityResult.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("check_pair chưa được implement (Sprint 3)")
