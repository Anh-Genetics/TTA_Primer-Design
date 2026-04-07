"""Module 2 — SequenceFetcher: lấy sequence từ NCBI Entrez.

Chức năng chính:
    - Fetch FASTA sequence (efetch db=nucleotide rettype=fasta)
    - Fetch GenBank record để lấy thông tin exon/CDS
    - Cache local (tránh re-fetch — sử dụng file-based cache đơn giản)
    - Hỗ trợ retry với exponential backoff khi gặp lỗi mạng/429

TODO (Sprint 1):
    - Implement ``fetch_fasta`` với Biopython Entrez
    - Implement ``fetch_genbank`` và ``extract_exon_coords``
    - Implement SQLite cache
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from tta_primer_design.config import AppConfig

if TYPE_CHECKING:
    pass  # Biopython SeqRecord chỉ dùng khi có biopython

logger = logging.getLogger("tta_primer_design.modules.sequence_fetcher")


class SequenceFetcher:
    """Lấy sequence từ NCBI Entrez và quản lý cache.

    Args:
        config: AppConfig đã load.

    Example::

        fetcher = SequenceFetcher(config)
        record = fetcher.fetch_fasta("NM_001101")
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config
        self._cache: dict[str, object] = {}

    def fetch_fasta(self, accession: str) -> object:
        """Fetch FASTA record từ NCBI.

        Args:
            accession: NCBI accession (ví dụ: NM_001101).

        Returns:
            BioPython SeqRecord.

        Raises:
            NotImplementedError: Chưa implement.
        """
        raise NotImplementedError("fetch_fasta chưa được implement (Sprint 1)")

    def fetch_genbank(self, accession: str) -> object:
        """Fetch GenBank record (có chứa feature annotations).

        Args:
            accession: NCBI accession.

        Returns:
            BioPython SeqRecord với features đầy đủ.

        Raises:
            NotImplementedError: Chưa implement.
        """
        raise NotImplementedError("fetch_genbank chưa được implement (Sprint 1)")

    def extract_exon_coords(self, genbank_record: object) -> list[tuple[int, int]]:
        """Trích xuất toạ độ exon từ GenBank record.

        Args:
            genbank_record: BioPython SeqRecord (GenBank format).

        Returns:
            Danh sách (start, end) 0-based cho mỗi exon.

        Raises:
            NotImplementedError: Chưa implement.
        """
        raise NotImplementedError("extract_exon_coords chưa được implement (Sprint 1)")

    def get_mrna_sequence(self, accession: str) -> object:
        """Lấy chuỗi mRNA (chỉ exon).

        Args:
            accession: NCBI accession.

        Returns:
            BioPython SeqRecord.

        Raises:
            NotImplementedError: Chưa implement.
        """
        raise NotImplementedError("get_mrna_sequence chưa được implement (Sprint 1)")

    def cache_sequence(self, accession: str, record: object) -> None:
        """Lưu record vào cache.

        Args:
            accession: Key cache.
            record: BioPython SeqRecord cần lưu.
        """
        self._cache[accession] = record
        logger.debug("Cached sequence for accession: %s", accession)

    def load_from_cache(self, accession: str) -> object | None:
        """Lấy record từ cache.

        Args:
            accession: Key cache.

        Returns:
            SeqRecord nếu tìm thấy, None nếu không có.
        """
        record = self._cache.get(accession)
        if record is not None:
            logger.debug("Cache hit for accession: %s", accession)
        return record
