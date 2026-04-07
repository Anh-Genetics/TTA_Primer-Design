"""Module 8 — SNPChecker: kiểm tra SNP trong vùng primer/probe.

Nguồn dữ liệu:
    A) NCBI Variation Services API (online)
       Endpoint: https://api.ncbi.nlm.nih.gov/variation/v0/
    B) dbSNP VCF file local (khuyến nghị cho batch)
       Dùng pysam hoặc cyvcf2 để parse

Logic:
    1. Lấy genomic coordinates của primer/probe
    2. Query dbSNP cho region đó
    3. Flag SNPs trong primer với MAF > threshold (default: 0.01)
    4. Cảnh báo đặc biệt SNP ở 3' end primer
    5. Cảnh báo tất cả SNP trong probe

Tiêu chí đánh giá:
    - "PASS"    : không có SNP
    - "WARNING" : có SNP nhưng MAF < threshold hoặc ở giữa primer
    - "FAIL"    : SNP ở 3' end hoặc MAF > threshold trong probe

TODO (Sprint 3):
    - Implement NCBI Variation API wrapper
    - Implement VCF-based lookup
    - Implement coordinate mapping (transcript → genomic)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from tta_primer_design.config import AppConfig
from tta_primer_design.models import PrimerPair

logger = logging.getLogger("tta_primer_design.modules.snp_checker")

_DEFAULT_MAF_THRESHOLD = 0.01


@dataclass
class SNPInfo:
    """Thông tin về một SNP trong vùng oligo.

    Attributes:
        rsid: ID dbSNP (ví dụ: rs12345).
        position_in_oligo: Vị trí SNP trong oligo (0-based từ 5').
        maf: Minor Allele Frequency (0.0–1.0).
        alleles: Allele string (ví dụ: "A/G").
        is_3prime: True nếu SNP nằm ở 3 bp cuối primer.
    """

    rsid: str
    position_in_oligo: int = 0
    maf: float = 0.0
    alleles: str = ""
    is_3prime: bool = False


@dataclass
class SNPResult:
    """Kết quả kiểm tra SNP cho một oligo.

    Attributes:
        oligo_id: ID oligo (primer/probe).
        has_snp: True nếu có bất kỳ SNP nào.
        snp_list: Danh sách tất cả SNP.
        critical_snps: SNP nghiêm trọng (3' end hoặc trong probe).
        recommendation: "PASS" | "WARNING" | "FAIL".
    """

    oligo_id: str
    has_snp: bool = False
    snp_list: list[SNPInfo] = field(default_factory=list)
    critical_snps: list[SNPInfo] = field(default_factory=list)
    recommendation: str = "PASS"


class SNPChecker:
    """Kiểm tra SNP trong vùng primer/probe từ dbSNP.

    Args:
        config: AppConfig đã load.
        maf_threshold: Ngưỡng MAF để flag SNP (default: 0.01).

    Example::

        checker = SNPChecker(config)
        pairs_flagged = checker.check_all(primer_pairs, seq_record)
    """

    def __init__(
        self,
        config: AppConfig,
        maf_threshold: float = _DEFAULT_MAF_THRESHOLD,
    ) -> None:
        self.config = config
        self.maf_threshold = maf_threshold

    def check_all(
        self,
        primer_pairs: list[PrimerPair],
        seq_record: object,
    ) -> list[PrimerPair]:
        """Kiểm tra SNP cho tất cả primer pairs.

        Gán flags vào ``pair.snp_flags``.

        Args:
            primer_pairs: Danh sách PrimerPair.
            seq_record: BioPython SeqRecord (có chứa accession/coords).

        Returns:
            Danh sách PrimerPair đã được gán snp_flags.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("check_all chưa được implement (Sprint 3)")

    def check_oligo(
        self,
        oligo_sequence: str,
        genomic_start: int,
        genomic_end: int,
        is_probe: bool = False,
    ) -> SNPResult:
        """Kiểm tra SNP cho một oligo.

        Args:
            oligo_sequence: Chuỗi oligo.
            genomic_start: Vị trí genomic bắt đầu (0-based).
            genomic_end: Vị trí genomic kết thúc.
            is_probe: True nếu là probe (strictter rules).

        Returns:
            SNPResult.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("check_oligo chưa được implement (Sprint 3)")
