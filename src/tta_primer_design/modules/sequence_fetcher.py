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
import time
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
        """
        cached = self.load_from_cache(accession)
        if cached is not None:
            return cached
        try:
            from Bio import Entrez, SeqIO
        except ImportError as exc:
            raise ImportError("biopython is required for NCBI fetch") from exc
        Entrez.email = self.config.ncbi.email
        if self.config.ncbi.api_key:
            Entrez.api_key = self.config.ncbi.api_key
        logger.info("Fetching FASTA for accession: %s", accession)
        last_exc: Exception | None = None
        for attempt in range(self.config.ncbi.retries):
            try:
                if attempt > 0:
                    time.sleep(min(2**attempt, 30))
                handle = Entrez.efetch(
                    db="nucleotide", id=accession, rettype="fasta", retmode="text"
                )
                record = SeqIO.read(handle, "fasta")
                handle.close()
                self.cache_sequence(accession, record)
                return record
            except Exception as exc:
                last_exc = exc
                logger.warning(
                    "Attempt %d/%d failed for %s: %s",
                    attempt + 1,
                    self.config.ncbi.retries,
                    accession,
                    exc,
                )
        raise RuntimeError(
            f"Failed to fetch {accession} after {self.config.ncbi.retries} attempts"
        ) from last_exc

    def fetch_genbank(self, accession: str) -> object:
        """Fetch GenBank record (có chứa feature annotations).

        Args:
            accession: NCBI accession.

        Returns:
            BioPython SeqRecord với features đầy đủ.
        """
        cache_key = f"gb_{accession}"
        cached = self.load_from_cache(cache_key)
        if cached is not None:
            return cached
        try:
            from Bio import Entrez, SeqIO
        except ImportError as exc:
            raise ImportError("biopython is required for NCBI fetch") from exc
        Entrez.email = self.config.ncbi.email
        if self.config.ncbi.api_key:
            Entrez.api_key = self.config.ncbi.api_key
        logger.info("Fetching GenBank for accession: %s", accession)
        last_exc: Exception | None = None
        for attempt in range(self.config.ncbi.retries):
            try:
                if attempt > 0:
                    time.sleep(min(2**attempt, 30))
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                handle.close()
                self.cache_sequence(cache_key, record)
                return record
            except Exception as exc:
                last_exc = exc
                logger.warning(
                    "Attempt %d/%d failed for %s: %s",
                    attempt + 1,
                    self.config.ncbi.retries,
                    accession,
                    exc,
                )
        raise RuntimeError(
            f"Failed to fetch {accession} after {self.config.ncbi.retries} attempts"
        ) from last_exc

    def extract_exon_coords(self, genbank_record: object) -> list[tuple[int, int]]:
        """Trích xuất toạ độ exon từ GenBank record.

        Args:
            genbank_record: BioPython SeqRecord (GenBank format).

        Returns:
            Danh sách (start, end) 0-based cho mỗi exon.
        """
        exon_coords: list[tuple[int, int]] = []
        try:
            for feature in genbank_record.features:
                if feature.type == "exon":
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    exon_coords.append((start, end))
        except AttributeError:
            pass
        return exon_coords

    def get_mrna_sequence(self, accession: str) -> object:
        """Lấy chuỗi mRNA (chỉ exon).

        Args:
            accession: NCBI accession.

        Returns:
            BioPython SeqRecord.
        """
        cached = self.load_from_cache(accession)
        if cached is not None:
            return cached
        return self.fetch_fasta(accession)

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
