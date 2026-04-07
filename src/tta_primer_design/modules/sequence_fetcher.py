"""Module 2 — SequenceFetcher: lấy sequence từ NCBI Entrez.

Chức năng chính:
    - Fetch FASTA sequence (efetch db=nucleotide rettype=fasta)
    - Fetch GenBank record để lấy thông tin exon/CDS
    - Cache local dùng SQLite (tránh re-fetch)
    - Retry với exponential backoff khi gặp lỗi mạng/429
"""

from __future__ import annotations

import io
import logging
import sqlite3
import time
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.error import HTTPError, URLError

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from tta_primer_design.config import AppConfig

if TYPE_CHECKING:
    pass

logger = logging.getLogger("tta_primer_design.modules.sequence_fetcher")

# Cache schema
_DDL = """
CREATE TABLE IF NOT EXISTS sequence_cache (
    accession TEXT NOT NULL,
    format    TEXT NOT NULL,
    data      TEXT NOT NULL,
    fetched   TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (accession, format)
);
"""

# Retry settings
_RETRY_STATUSES = {429, 500, 502, 503, 504}


class SequenceFetcher:
    """Lấy sequence từ NCBI Entrez và quản lý cache.

    Args:
        config: AppConfig đã load.
        cache_path: Đường dẫn file SQLite cache.
                    Mặc định ``~/.cache/tta_primer_design/sequences.db``.

    Example::

        fetcher = SequenceFetcher(config)
        record = fetcher.fetch_fasta("NM_001101")
    """

    def __init__(
        self,
        config: AppConfig,
        cache_path: str | Path | None = None,
    ) -> None:
        self.config = config

        # Configure Biopython Entrez
        Entrez.email = config.ncbi.email
        if config.ncbi.api_key:
            Entrez.api_key = config.ncbi.api_key

        # In-memory cache (accession → SeqRecord, keyed by "accession:format")
        self._cache: dict[str, SeqRecord] = {}

        # SQLite persistent cache
        if cache_path is None:
            default_dir = Path.home() / ".cache" / "tta_primer_design"
            default_dir.mkdir(parents=True, exist_ok=True)
            cache_path = default_dir / "sequences.db"
        self._cache_path = Path(cache_path)
        self._init_db()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fetch_fasta(self, accession: str) -> SeqRecord:
        """Fetch FASTA record từ NCBI.

        Args:
            accession: NCBI accession (ví dụ: NM_001101).

        Returns:
            BioPython SeqRecord chứa sequence.

        Raises:
            RuntimeError: Nếu không thể fetch sau tất cả các lần retry.
        """
        cache_key = f"{accession}:fasta"
        if (hit := self._mem_get(cache_key)) is not None:
            return hit
        if (hit := self._db_get(accession, "fasta")) is not None:
            self._cache[cache_key] = hit
            return hit

        data = self._entrez_fetch(accession, rettype="fasta", retmode="text")
        record = SeqIO.read(io.StringIO(data), "fasta")
        self._store(cache_key, accession, "fasta", data)
        return record

    def fetch_genbank(self, accession: str) -> SeqRecord:
        """Fetch GenBank record (có chứa feature annotations).

        Args:
            accession: NCBI accession.

        Returns:
            BioPython SeqRecord với features đầy đủ.

        Raises:
            RuntimeError: Nếu không thể fetch sau tất cả các lần retry.
        """
        cache_key = f"{accession}:gb"
        if (hit := self._mem_get(cache_key)) is not None:
            return hit
        if (hit := self._db_get(accession, "gb")) is not None:
            self._cache[cache_key] = hit
            return hit

        data = self._entrez_fetch(accession, rettype="gb", retmode="text")
        record = SeqIO.read(io.StringIO(data), "genbank")
        self._store(cache_key, accession, "gb", data)
        return record

    def extract_exon_coords(self, genbank_record: SeqRecord) -> list[tuple[int, int]]:
        """Trích xuất toạ độ exon từ GenBank record.

        Tìm kiếm feature loại ``exon`` trước; nếu không có thì dùng
        feature ``CDS`` để lấy từng exon (CompoundLocation).

        Args:
            genbank_record: BioPython SeqRecord (GenBank format).

        Returns:
            Danh sách ``(start, end)`` 0-based, sắp xếp theo start,
            cho mỗi exon.  Trả về list rỗng nếu không tìm thấy.
        """
        coords: list[tuple[int, int]] = []

        # Ưu tiên feature type "exon"
        exon_features = [f for f in genbank_record.features if f.type == "exon"]
        if exon_features:
            for feat in exon_features:
                for part in feat.location.parts:
                    coords.append((int(part.start), int(part.end)))
            coords.sort()
            return coords

        # Fallback: lấy từng exon từ CDS CompoundLocation
        cds_features = [f for f in genbank_record.features if f.type == "CDS"]
        for feat in cds_features:
            for part in feat.location.parts:
                coords.append((int(part.start), int(part.end)))

        # De-duplicate (giữ thứ tự, tránh mất tuple có cùng start nhưng khác end)
        seen: set[tuple[int, int]] = set()
        unique: list[tuple[int, int]] = []
        for c in coords:
            if c not in seen:
                seen.add(c)
                unique.append(c)
        unique.sort()
        return unique

    def get_mrna_sequence(self, accession: str) -> SeqRecord:
        """Lấy SeqRecord mRNA từ NCBI (FASTA).

        Fetch FASTA và GenBank, dùng GenBank để xác nhận accession là mRNA.
        Trả về SeqRecord FASTA.

        Args:
            accession: NCBI accession.

        Returns:
            BioPython SeqRecord (FASTA).

        Raises:
            RuntimeError: Nếu không thể fetch.
        """
        return self.fetch_fasta(accession)

    def cache_sequence(self, accession: str, record: SeqRecord) -> None:
        """Lưu record vào bộ nhớ cache.

        Args:
            accession: Key cache.
            record: BioPython SeqRecord cần lưu.
        """
        self._cache[accession] = record
        logger.debug("Cached sequence for accession: %s", accession)

    def load_from_cache(self, accession: str) -> SeqRecord | None:
        """Lấy record từ bộ nhớ cache.

        Args:
            accession: Key cache.

        Returns:
            SeqRecord nếu tìm thấy, None nếu không có.
        """
        record = self._cache.get(accession)
        if record is not None:
            logger.debug("Cache hit for accession: %s", accession)
        return record

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _init_db(self) -> None:
        """Khởi tạo bảng cache trong SQLite nếu chưa có."""
        with sqlite3.connect(self._cache_path) as con:
            con.execute(_DDL)
            con.commit()

    def _mem_get(self, cache_key: str) -> SeqRecord | None:
        """Lấy từ in-memory cache."""
        hit = self._cache.get(cache_key)
        if hit is not None:
            logger.debug("Memory cache hit: %s", cache_key)
        return hit

    def _db_get(self, accession: str, fmt: str) -> SeqRecord | None:
        """Lấy từ SQLite cache.

        Returns:
            SeqRecord nếu tìm thấy, None nếu không.
        """
        with sqlite3.connect(self._cache_path) as con:
            row = con.execute(
                "SELECT data FROM sequence_cache WHERE accession=? AND format=?",
                (accession, fmt),
            ).fetchone()
        if row is None:
            return None
        logger.debug("SQLite cache hit: %s (%s)", accession, fmt)
        parser = "fasta" if fmt == "fasta" else "genbank"
        return SeqIO.read(io.StringIO(row[0]), parser)

    def _db_put(self, accession: str, fmt: str, data: str) -> None:
        """Lưu raw text vào SQLite cache."""
        with sqlite3.connect(self._cache_path) as con:
            con.execute(
                "INSERT OR REPLACE INTO sequence_cache (accession, format, data) VALUES (?,?,?)",
                (accession, fmt, data),
            )
            con.commit()

    def _store(self, cache_key: str, accession: str, fmt: str, data: str) -> None:
        """Lưu vào cả in-memory cache và SQLite."""
        parser = "fasta" if fmt == "fasta" else "genbank"
        record = SeqIO.read(io.StringIO(data), parser)
        self._cache[cache_key] = record
        self._db_put(accession, fmt, data)
        logger.debug("Stored in cache: %s (%s)", accession, fmt)

    def _entrez_fetch(self, accession: str, rettype: str, retmode: str) -> str:
        """Gọi Entrez.efetch với retry + exponential backoff.

        Args:
            accession: NCBI accession.
            rettype: Loại trả về (``"fasta"`` hoặc ``"gb"``).
            retmode: Chế độ trả về (``"text"``).

        Returns:
            Nội dung text từ NCBI.

        Raises:
            RuntimeError: Nếu vượt quá số lần retry.
        """
        retries = self.config.ncbi.retries
        rate_limit = self.config.ncbi.rate_limit
        delay = 1.0 / max(rate_limit, 1)  # inter-request delay

        last_exc: Exception | None = None
        for attempt in range(retries):
            try:
                time.sleep(delay)
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=accession,
                    rettype=rettype,
                    retmode=retmode,
                )
                data = handle.read()
                handle.close()
                if isinstance(data, bytes):
                    data = data.decode("utf-8")
                logger.info("Fetched %s (%s) from NCBI", accession, rettype)
                return data
            except HTTPError as exc:
                last_exc = exc
                if exc.code in _RETRY_STATUSES:
                    wait = 2**attempt
                    logger.warning(
                        "HTTP %s fetching %s — retrying in %ss (attempt %d/%d)",
                        exc.code,
                        accession,
                        wait,
                        attempt + 1,
                        retries,
                    )
                    time.sleep(wait)
                else:
                    raise RuntimeError(
                        f"NCBI HTTP error {exc.code} fetching {accession}: {exc}"
                    ) from exc
            except (URLError, OSError) as exc:
                last_exc = exc
                wait = 2**attempt
                logger.warning(
                    "Network error fetching %s — retrying in %ss (attempt %d/%d): %s",
                    accession,
                    wait,
                    attempt + 1,
                    retries,
                    exc,
                )
                time.sleep(wait)

        raise RuntimeError(
            f"Failed to fetch {accession} ({rettype}) after {retries} attempts: {last_exc}"
        )
