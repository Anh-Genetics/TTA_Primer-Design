"""Tests cho SequenceFetcher — fetch FASTA/GenBank từ NCBI Entrez, cache SQLite."""

from __future__ import annotations

import io
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch
from urllib.error import HTTPError, URLError

import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from tta_primer_design.config import AppConfig
from tta_primer_design.modules.sequence_fetcher import SequenceFetcher

# ---------------------------------------------------------------------------
# Fixtures — sample NCBI text responses
# ---------------------------------------------------------------------------

_FASTA_TEXT = textwrap.dedent(
    """\
    >NM_001101.5 Homo sapiens actin beta (ACTB), mRNA
    ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAAGCCGGCTTCGCGGGCG
    ACGATGCCCCGAGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCACCAGGGCGTGATGGTGGGCAT
    """
)

_GENBANK_TEXT = textwrap.dedent(
    """\
    LOCUS       NM_001101               5030 bp    mRNA    linear   PRI 01-JAN-2023
    DEFINITION  Homo sapiens actin beta (ACTB), mRNA.
    ACCESSION   NM_001101
    VERSION     NM_001101.5
    FEATURES             Location/Qualifiers
         source          1..5030
                         /organism="Homo sapiens"
         exon            1..150
                         /number=1
         exon            300..500
                         /number=2
         CDS             join(1..150,300..500)
                         /gene="ACTB"
                         /product="actin, cytoplasmic 1"
    ORIGIN
            1 atggatgatg atatcgccgc gctcgtcgtc gacaacggct ccggcatgtg caaagccggc
           61 ttcgcgggcg acgatgcccc gagggccgtc ttcccctcca tcgtggggcg ccccaggcac
          121 cagggcgtga tggtgggcat gggtcagaag gattcctatg tgggcgacga ggccgagggc
          181 aagctggagt acgagctgcc tgacggccag gtcatcacca ttggcaatga gcggttccgc
          241 tgccctgagg cactctttgg tgatgatgcc ctggagcatc tgatggccac tgccgccacc
          301 atgaagatca agatcattgc tcctcctgag cgcaagtact ccgtgtggat cggcggctcc
          361 atcctggcct cgctgtccac cttccagcag atgtggatca gcaagcagga gtatgacgag
          421 tccggcccct ccatcgtcca ccgcaaatgc ttctaagccc agggcaacgg catcaaggcc
          481 ttctaccggc tcttccagcc atcatcatgt
    //
    """
)


def _make_config(tmp_path: Path) -> AppConfig:
    cfg = AppConfig()
    cfg.ncbi.email = "test@example.com"
    cfg.ncbi.retries = 3
    cfg.ncbi.timeout = 10
    cfg.ncbi.rate_limit = 10  # fast in tests
    return cfg


def _make_fetcher(tmp_path: Path) -> SequenceFetcher:
    cfg = _make_config(tmp_path)
    return SequenceFetcher(cfg, cache_path=tmp_path / "test_cache.db")


# ---------------------------------------------------------------------------
# Tests: fetch_fasta
# ---------------------------------------------------------------------------


class TestFetchFasta:
    def test_returns_seqrecord(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch("tta_primer_design.modules.sequence_fetcher.Entrez.efetch", return_value=mock_handle):
            record = fetcher.fetch_fasta("NM_001101")

        assert isinstance(record, SeqRecord)
        assert "NM_001101" in record.id

    def test_sequence_content(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch("tta_primer_design.modules.sequence_fetcher.Entrez.efetch", return_value=mock_handle):
            record = fetcher.fetch_fasta("NM_001101")

        assert str(record.seq).startswith("ATGGATGATG")

    def test_memory_cache_hit_skips_network(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch",
            return_value=mock_handle,
        ) as mock_efetch:
            fetcher.fetch_fasta("NM_001101")
            fetcher.fetch_fasta("NM_001101")  # second call — should use cache

        mock_efetch.assert_called_once()

    def test_sqlite_cache_persists(self, tmp_path: Path) -> None:
        cache_db = tmp_path / "test_cache.db"
        cfg = _make_config(tmp_path)

        # First fetcher populates cache
        fetcher1 = SequenceFetcher(cfg, cache_path=cache_db)
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]
        with patch("tta_primer_design.modules.sequence_fetcher.Entrez.efetch", return_value=mock_handle):
            fetcher1.fetch_fasta("NM_001101")

        # Second fetcher reads from SQLite without calling NCBI
        fetcher2 = SequenceFetcher(cfg, cache_path=cache_db)
        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch"
        ) as mock_efetch:
            record = fetcher2.fetch_fasta("NM_001101")

        mock_efetch.assert_not_called()
        assert isinstance(record, SeqRecord)

    def test_bytes_response_decoded(self, tmp_path: Path) -> None:
        """Entrez đôi khi trả về bytes — phải decode."""
        fetcher = _make_fetcher(tmp_path)
        mock_handle = MagicMock()
        mock_handle.read.return_value = _FASTA_TEXT.encode("utf-8")

        with patch("tta_primer_design.modules.sequence_fetcher.Entrez.efetch", return_value=mock_handle):
            record = fetcher.fetch_fasta("NM_001101")

        assert isinstance(record, SeqRecord)


# ---------------------------------------------------------------------------
# Tests: fetch_genbank
# ---------------------------------------------------------------------------


class TestFetchGenbank:
    def test_returns_seqrecord_with_features(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        mock_handle = io.StringIO(_GENBANK_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch("tta_primer_design.modules.sequence_fetcher.Entrez.efetch", return_value=mock_handle):
            record = fetcher.fetch_genbank("NM_001101")

        assert isinstance(record, SeqRecord)
        assert len(record.features) > 0

    def test_genbank_memory_cache(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        mock_handle = io.StringIO(_GENBANK_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch",
            return_value=mock_handle,
        ) as mock_efetch:
            fetcher.fetch_genbank("NM_001101")
            fetcher.fetch_genbank("NM_001101")

        mock_efetch.assert_called_once()


# ---------------------------------------------------------------------------
# Tests: extract_exon_coords
# ---------------------------------------------------------------------------


class TestExtractExonCoords:
    def _get_record(self) -> SeqRecord:
        return SeqIO.read(io.StringIO(_GENBANK_TEXT), "genbank")

    def test_returns_list_of_tuples(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        record = self._get_record()
        coords = fetcher.extract_exon_coords(record)
        assert isinstance(coords, list)
        assert all(isinstance(c, tuple) and len(c) == 2 for c in coords)

    def test_exon_coords_values(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        record = self._get_record()
        coords = fetcher.extract_exon_coords(record)
        # Record has exon 1..150 and exon 300..500 (0-based: 0..150, 299..500)
        assert len(coords) >= 2
        starts = [c[0] for c in coords]
        assert starts == sorted(starts), "Coords should be sorted by start"

    def test_exon_coords_sorted(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        record = self._get_record()
        coords = fetcher.extract_exon_coords(record)
        assert coords == sorted(coords)

    def test_no_features_returns_empty(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        record = SeqIO.read(io.StringIO(_FASTA_TEXT), "fasta")
        coords = fetcher.extract_exon_coords(record)
        assert coords == []


# ---------------------------------------------------------------------------
# Tests: get_mrna_sequence
# ---------------------------------------------------------------------------


class TestGetMrnaSequence:
    def test_returns_seqrecord(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch("tta_primer_design.modules.sequence_fetcher.Entrez.efetch", return_value=mock_handle):
            record = fetcher.get_mrna_sequence("NM_001101")

        assert isinstance(record, SeqRecord)


# ---------------------------------------------------------------------------
# Tests: cache helpers
# ---------------------------------------------------------------------------


class TestCacheHelpers:
    def test_cache_sequence_and_load(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        fasta_record = SeqIO.read(io.StringIO(_FASTA_TEXT), "fasta")
        fetcher.cache_sequence("NM_001101", fasta_record)
        result = fetcher.load_from_cache("NM_001101")
        assert result is fasta_record

    def test_load_from_cache_miss_returns_none(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        result = fetcher.load_from_cache("NOTEXIST")
        assert result is None


# ---------------------------------------------------------------------------
# Tests: retry behaviour
# ---------------------------------------------------------------------------


def _http_error(code: int, msg: str) -> HTTPError:
    """Helper to build an HTTPError without type-ignore comments."""
    fp = io.BytesIO(b"")
    return HTTPError(url="https://eutils.ncbi.nlm.nih.gov/", code=code, msg=msg, hdrs={}, fp=fp)  # type: ignore[arg-type]


class TestRetry:
    def test_retries_on_http_429(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        fetcher.config.ncbi.retries = 3

        error_429 = _http_error(429, "Too Many Requests")
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch",
            side_effect=[error_429, error_429, mock_handle],
        ), patch("tta_primer_design.modules.sequence_fetcher.time.sleep"):
            record = fetcher.fetch_fasta("NM_001101")

        assert isinstance(record, SeqRecord)

    def test_raises_after_exhausting_retries(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        fetcher.config.ncbi.retries = 2

        error_429 = _http_error(429, "Too Many Requests")

        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch",
            side_effect=[error_429, error_429, error_429],
        ), patch("tta_primer_design.modules.sequence_fetcher.time.sleep"):
            with pytest.raises(RuntimeError, match="Failed to fetch"):
                fetcher.fetch_fasta("NM_001101")

    def test_raises_on_non_retryable_http_error(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        error_404 = _http_error(404, "Not Found")

        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch",
            side_effect=error_404,
        ), patch("tta_primer_design.modules.sequence_fetcher.time.sleep"):
            with pytest.raises(RuntimeError, match="NCBI HTTP error 404"):
                fetcher.fetch_fasta("INVALID")

    def test_retries_on_network_error(self, tmp_path: Path) -> None:
        fetcher = _make_fetcher(tmp_path)
        fetcher.config.ncbi.retries = 3

        url_error = URLError("Network unreachable")
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch",
            side_effect=[url_error, mock_handle],
        ), patch("tta_primer_design.modules.sequence_fetcher.time.sleep"):
            record = fetcher.fetch_fasta("NM_001101")

        assert isinstance(record, SeqRecord)

    def test_rate_limit_sleep_called(self, tmp_path: Path) -> None:
        """Kiểm tra sleep được gọi giữa các request để tuân rate limit."""
        fetcher = _make_fetcher(tmp_path)
        mock_handle = io.StringIO(_FASTA_TEXT)
        mock_handle.close = lambda: None  # type: ignore[assignment]

        with patch(
            "tta_primer_design.modules.sequence_fetcher.Entrez.efetch",
            return_value=mock_handle,
        ), patch("tta_primer_design.modules.sequence_fetcher.time.sleep") as mock_sleep:
            fetcher.fetch_fasta("NM_001101")

        mock_sleep.assert_called()
