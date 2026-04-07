"""Module 2 — SequenceFetcher: lấy sequence từ NCBI Entrez.

Chức năng chính:
    - Fetch FASTA sequence (efetch db=nucleotide rettype=fasta)
    - Fetch GenBank record để lấy thông tin exon/CDS
    - Cache hai tầng: in-memory dict + SQLite (tuỳ chọn)
    - Retry với exponential backoff khi gặp lỗi mạng/429 (qua ncbi_api)
"""

from __future__ import annotations

import io
import logging
import sqlite3
import time
from pathlib import Path
from typing import TYPE_CHECKING

from tta_primer_design.config import AppConfig
from tta_primer_design.utils.ncbi_api import entrez_efetch, setup_entrez

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("tta_primer_design.modules.sequence_fetcher")

_DB_SCHEMA = """
CREATE TABLE IF NOT EXISTS sequence_cache (
    accession  TEXT PRIMARY KEY,
    format     TEXT NOT NULL,
    content    TEXT NOT NULL,
    fetched_at REAL NOT NULL
);
"""


class SequenceFetcher:
    """Lấy sequence từ NCBI Entrez và quản lý cache (in-memory + SQLite).

    Args:
        config:    AppConfig đã load.
        cache_dir: Thư mục lưu SQLite database.
                   ``None`` → chỉ sử dụng in-memory cache.

    Example::

        fetcher = SequenceFetcher(config)
        record = fetcher.fetch_fasta("NM_001101")

        # Với SQLite cache bền vững:
        fetcher = SequenceFetcher(config, cache_dir="data/cache")
    """

    def __init__(
        self,
        config: AppConfig,
        cache_dir: Path | str | None = None,
    ) -> None:
        self.config = config
        self._cache: dict[str, object] = {}
        self._db_conn: sqlite3.Connection | None = None

        setup_entrez(config.ncbi.email, config.ncbi.api_key)

        if cache_dir is not None:
            cache_path = Path(cache_dir)
            cache_path.mkdir(parents=True, exist_ok=True)
            db_file = cache_path / "sequences.db"
            self._db_conn = sqlite3.connect(str(db_file), check_same_thread=False)
            self._db_conn.execute(_DB_SCHEMA)
            self._db_conn.commit()
            logger.debug("SQLite cache opened: %s", db_file)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fetch_fasta(self, accession: str) -> SeqRecord:
        """Fetch FASTA record từ NCBI.

        Kết quả được lưu vào cache (in-memory và SQLite nếu có).
        Lần gọi tiếp theo với cùng accession sẽ trả về từ cache mà không
        gọi mạng.

        Args:
            accession: NCBI accession (ví dụ: ``NM_001101``).

        Returns:
            BioPython ``SeqRecord`` (FASTA).

        Raises:
            RuntimeError: Nếu fetch thất bại sau tất cả các lần thử.
        """
        from Bio import SeqIO

        cached = self._load_record(accession)
        if cached is not None:
            return cached

        handle = entrez_efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta",
            retmode="text",
            retries=self.config.ncbi.retries,
            timeout=self.config.ncbi.timeout,
        )
        try:
            record = SeqIO.read(handle, "fasta")
        finally:
            handle.close()

        self._store_record(accession, record, fmt="fasta")
        return record

    def fetch_genbank(self, accession: str) -> SeqRecord:
        """Fetch GenBank record (có chứa feature annotations).

        Args:
            accession: NCBI accession.

        Returns:
            BioPython ``SeqRecord`` với đầy đủ features (exon, CDS, …).

        Raises:
            RuntimeError: Nếu fetch thất bại sau tất cả các lần thử.
        """
        from Bio import SeqIO

        cached = self._load_record(accession)
        if cached is not None:
            return cached

        handle = entrez_efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text",
            retries=self.config.ncbi.retries,
            timeout=self.config.ncbi.timeout,
        )
        try:
            record = SeqIO.read(handle, "genbank")
        finally:
            handle.close()

        self._store_record(accession, record, fmt="genbank")
        return record

    def extract_exon_coords(self, genbank_record: SeqRecord) -> list[tuple[int, int]]:
        """Trích xuất toạ độ exon từ GenBank record.

        Thứ tự ưu tiên:
        1. Feature có type ``exon``.
        2. Sub-parts của feature ``mRNA`` (joined location).
        3. Sub-parts của feature ``CDS``.
        4. Fallback: toàn bộ sequence coi như một exon.

        Args:
            genbank_record: BioPython ``SeqRecord`` (GenBank format).

        Returns:
            Danh sách ``(start, end)`` 0-based (Python slice convention),
            đã sắp xếp theo start.
        """
        coords: list[tuple[int, int]] = []

        # 1. explicit exon features
        for feat in genbank_record.features:
            if feat.type == "exon":
                for part in feat.location.parts:
                    coords.append((int(part.start), int(part.end)))

        if coords:
            coords.sort()
            return coords

        # 2. mRNA feature sub-parts (joined locations)
        for feat in genbank_record.features:
            if feat.type == "mRNA":
                for part in feat.location.parts:
                    coords.append((int(part.start), int(part.end)))

        if coords:
            coords.sort()
            return coords

        # 3. CDS sub-parts
        for feat in genbank_record.features:
            if feat.type == "CDS":
                for part in feat.location.parts:
                    coords.append((int(part.start), int(part.end)))

        if coords:
            coords.sort()
            return coords

        # 4. fallback: whole sequence as single exon
        seq_len = len(genbank_record.seq)
        logger.debug(
            "No exon/mRNA/CDS features in %s; treating entire sequence as one exon",
            genbank_record.id,
        )
        return [(0, seq_len)]

    def get_mrna_sequence(self, accession: str) -> SeqRecord:
        """Lấy chuỗi mRNA bằng cách ghép các exon từ GenBank record.

        Args:
            accession: NCBI accession.

        Returns:
            BioPython ``SeqRecord`` chứa chuỗi mRNA (chỉ exon, đã ghép).

        Raises:
            RuntimeError: Nếu fetch thất bại.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        gb_record = self.fetch_genbank(accession)
        exon_coords = self.extract_exon_coords(gb_record)

        mrna_seq = Seq("".join(str(gb_record.seq[s:e]) for s, e in exon_coords))
        mrna_record = SeqRecord(
            mrna_seq,
            id=accession,
            name=gb_record.name,
            description=f"mRNA (exon-joined) for {accession}",
        )
        return mrna_record

    def cache_sequence(self, accession: str, record: object) -> None:
        """Lưu record vào cache (in-memory và SQLite nếu có).

        Args:
            accession: Key cache.
            record: BioPython ``SeqRecord`` cần lưu.
        """
        self._cache[accession] = record
        logger.debug("Cached sequence for accession: %s", accession)

        if self._db_conn is not None:
            # Try GenBank first (preserves features); fall back to FASTA.
            for fmt in ("genbank", "fasta"):
                try:
                    self._store_record(accession, record, fmt=fmt)
                    break
                except Exception as exc:
                    logger.debug("Could not store %s as %s: %s", accession, fmt, exc)

    def load_from_cache(self, accession: str) -> object | None:
        """Lấy record từ cache (in-memory trước, rồi SQLite).

        Args:
            accession: Key cache.

        Returns:
            ``SeqRecord`` nếu tìm thấy, ``None`` nếu không có.
        """
        record = self._cache.get(accession)
        if record is not None:
            logger.debug("Cache hit (memory) for accession: %s", accession)
            return record

        if self._db_conn is not None:
            loaded = self._load_record(accession)
            if loaded is not None:
                self._cache[accession] = loaded
                return loaded

        return None

    def close(self) -> None:
        """Đóng kết nối SQLite cache (nếu có)."""
        if self._db_conn is not None:
            self._db_conn.close()
            self._db_conn = None
            logger.debug("SQLite cache connection closed")

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _store_record(self, accession: str, record: object, fmt: str) -> None:
        """Serialize *record* và lưu vào in-memory cache + SQLite."""
        from Bio import SeqIO

        # Always update in-memory cache.
        self._cache[accession] = record

        if self._db_conn is None:
            return

        buf = io.StringIO()
        SeqIO.write(record, buf, fmt)  # type: ignore[arg-type]
        content = buf.getvalue()

        self._db_conn.execute(
            """
            INSERT OR REPLACE INTO sequence_cache (accession, format, content, fetched_at)
            VALUES (?, ?, ?, ?)
            """,
            (accession, fmt, content, time.time()),
        )
        self._db_conn.commit()
        logger.debug("Stored %s in SQLite cache (format=%s)", accession, fmt)

    def _load_record(self, accession: str) -> SeqRecord | None:
        """Tải SeqRecord từ in-memory cache hoặc SQLite.

        Returns ``None`` if the accession is not found in either layer.
        """
        from Bio import SeqIO

        # In-memory first.
        record = self._cache.get(accession)
        if record is not None:
            logger.debug("Cache hit (memory) for %s", accession)
            return record  # type: ignore[return-value]

        # SQLite.
        if self._db_conn is not None:
            row = self._db_conn.execute(
                "SELECT content, format FROM sequence_cache WHERE accession = ?",
                (accession,),
            ).fetchone()
            if row is not None:
                content, stored_fmt = row
                record = SeqIO.read(io.StringIO(content), stored_fmt)
                self._cache[accession] = record  # promote to memory
                logger.debug("Cache hit (SQLite) for %s (format=%s)", accession, stored_fmt)
                return record  # type: ignore[return-value]

        return None
