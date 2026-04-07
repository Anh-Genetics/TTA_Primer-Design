"""Tests cho modules/sequence_fetcher.py.

Tất cả các lời gọi mạng tới NCBI được mock hoàn toàn.
"""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import patch

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.modules.sequence_fetcher import SequenceFetcher

# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

FASTA_TEXT = ">NM_001101.5 Homo sapiens ACTB mRNA\nATCGATCGATCGATCGATCG\n"

FIXTURES = Path(__file__).parent / "fixtures"


def _make_fasta_record(acc: str = "NM_001101", seq: str = "ATCGATCGATCGATCGATCG") -> SeqRecord:
    """Build a minimal FASTA SeqRecord."""
    record = SeqRecord(Seq(seq), id=acc, name=acc, description="test fasta")
    record.annotations["molecule_type"] = "DNA"
    return record


def _make_genbank_record(
    acc: str = "NM_001101",
    seq: str = "A" * 600,
    exon_ranges: list[tuple[int, int]] | None = None,
) -> SeqRecord:
    """Build a minimal GenBank SeqRecord with optional exon features."""
    record = SeqRecord(Seq(seq), id=acc, name=acc, description="test genbank")
    record.annotations["molecule_type"] = "mRNA"

    if exon_ranges:
        for start, end in exon_ranges:
            feat = SeqFeature(FeatureLocation(start, end, strand=1), type="exon")
            record.features.append(feat)

    return record


def _fasta_handle(record: SeqRecord) -> io.StringIO:
    """Serialise a SeqRecord as FASTA and return a StringIO handle."""
    from Bio import SeqIO

    buf = io.StringIO()
    SeqIO.write(record, buf, "fasta")
    buf.seek(0)
    return buf


def _genbank_handle(record: SeqRecord) -> io.StringIO:
    """Serialise a SeqRecord as GenBank and return a StringIO handle."""
    from Bio import SeqIO

    buf = io.StringIO()
    SeqIO.write(record, buf, "genbank")
    buf.seek(0)
    return buf


@pytest.fixture
def config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def fetcher(config: AppConfig) -> SequenceFetcher:
    with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
        return SequenceFetcher(config)


@pytest.fixture
def fetcher_with_cache(config: AppConfig, tmp_path: Path) -> SequenceFetcher:
    with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
        f = SequenceFetcher(config, cache_dir=tmp_path / "cache")
    return f


# ---------------------------------------------------------------------------
# Tests: fetch_fasta
# ---------------------------------------------------------------------------


class TestFetchFasta:
    """Test SequenceFetcher.fetch_fasta."""

    def test_returns_seq_record(self, fetcher: SequenceFetcher) -> None:
        record = _make_fasta_record()
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_fasta_handle(record),
        ):
            result = fetcher.fetch_fasta("NM_001101")
        assert isinstance(result, SeqRecord)

    def test_record_id_matches(self, fetcher: SequenceFetcher) -> None:
        record = _make_fasta_record("NM_001101")
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_fasta_handle(record),
        ):
            result = fetcher.fetch_fasta("NM_001101")
        assert result.id == "NM_001101"

    def test_calls_entrez_efetch(self, fetcher: SequenceFetcher) -> None:
        record = _make_fasta_record()
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_fasta_handle(record),
        ) as mock_fetch:
            fetcher.fetch_fasta("NM_001101")
        mock_fetch.assert_called_once()
        call_kwargs = mock_fetch.call_args.kwargs
        assert call_kwargs["rettype"] == "fasta"
        assert call_kwargs["db"] == "nucleotide"
        assert call_kwargs["id"] == "NM_001101"

    def test_second_call_uses_memory_cache(self, fetcher: SequenceFetcher) -> None:
        record = _make_fasta_record()
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_fasta_handle(record),
        ) as mock_fetch:
            fetcher.fetch_fasta("NM_001101")
            fetcher.fetch_fasta("NM_001101")  # second call — should hit cache
        # Entrez called only once
        mock_fetch.assert_called_once()

    def test_sequence_content(self, fetcher: SequenceFetcher) -> None:
        seq = "ATCGATCGATCGATCGATCG"
        record = _make_fasta_record(seq=seq)
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_fasta_handle(record),
        ):
            result = fetcher.fetch_fasta("NM_001101")
        assert str(result.seq) == seq


# ---------------------------------------------------------------------------
# Tests: fetch_genbank
# ---------------------------------------------------------------------------


class TestFetchGenbank:
    """Test SequenceFetcher.fetch_genbank."""

    def test_returns_seq_record(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record()
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_genbank_handle(record),
        ):
            result = fetcher.fetch_genbank("NM_001101")
        assert isinstance(result, SeqRecord)

    def test_calls_entrez_efetch_with_gb(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record()
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_genbank_handle(record),
        ) as mock_fetch:
            fetcher.fetch_genbank("NM_001101")
        call_kwargs = mock_fetch.call_args.kwargs
        assert call_kwargs["rettype"] == "gb"

    def test_second_call_uses_memory_cache(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record()
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_genbank_handle(record),
        ) as mock_fetch:
            fetcher.fetch_genbank("NM_001101")
            fetcher.fetch_genbank("NM_001101")
        mock_fetch.assert_called_once()

    def test_features_preserved(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record(exon_ranges=[(0, 100), (200, 400)])
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_genbank_handle(record),
        ):
            result = fetcher.fetch_genbank("NM_001101")
        exon_feats = [f for f in result.features if f.type == "exon"]
        assert len(exon_feats) == 2


# ---------------------------------------------------------------------------
# Tests: extract_exon_coords
# ---------------------------------------------------------------------------


class TestExtractExonCoords:
    """Test SequenceFetcher.extract_exon_coords."""

    def test_extracts_exon_features(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record(exon_ranges=[(0, 200), (300, 600)])
        coords = fetcher.extract_exon_coords(record)
        assert coords == [(0, 200), (300, 600)]

    def test_coords_sorted_by_start(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record(exon_ranges=[(500, 600), (0, 100), (200, 400)])
        coords = fetcher.extract_exon_coords(record)
        assert coords == [(0, 100), (200, 400), (500, 600)]

    def test_fallback_to_mrna_feature(self, fetcher: SequenceFetcher) -> None:
        from Bio.SeqFeature import CompoundLocation

        record = SeqRecord(Seq("A" * 600), id="NM_X", name="NM_X")
        record.annotations["molecule_type"] = "mRNA"
        # Joined mRNA feature
        loc = CompoundLocation(
            [FeatureLocation(0, 100, strand=1), FeatureLocation(200, 400, strand=1)]
        )
        record.features.append(SeqFeature(loc, type="mRNA"))
        coords = fetcher.extract_exon_coords(record)
        assert (0, 100) in coords
        assert (200, 400) in coords

    def test_fallback_to_cds_feature(self, fetcher: SequenceFetcher) -> None:
        from Bio.SeqFeature import CompoundLocation

        record = SeqRecord(Seq("A" * 600), id="NM_X", name="NM_X")
        record.annotations["molecule_type"] = "mRNA"
        loc = CompoundLocation(
            [FeatureLocation(50, 150, strand=1), FeatureLocation(250, 500, strand=1)]
        )
        record.features.append(SeqFeature(loc, type="CDS"))
        coords = fetcher.extract_exon_coords(record)
        assert (50, 150) in coords
        assert (250, 500) in coords

    def test_fallback_full_sequence(self, fetcher: SequenceFetcher) -> None:
        record = SeqRecord(Seq("A" * 300), id="NM_X", name="NM_X")
        record.annotations["molecule_type"] = "mRNA"
        coords = fetcher.extract_exon_coords(record)
        assert coords == [(0, 300)]

    def test_single_exon_feature(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record(exon_ranges=[(0, 500)])
        coords = fetcher.extract_exon_coords(record)
        assert coords == [(0, 500)]


# ---------------------------------------------------------------------------
# Tests: get_mrna_sequence
# ---------------------------------------------------------------------------


class TestGetMrnaSequence:
    """Test SequenceFetcher.get_mrna_sequence."""

    def test_returns_seq_record(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record(exon_ranges=[(0, 100), (200, 300)])
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_genbank_handle(record),
        ):
            mrna = fetcher.get_mrna_sequence("NM_001101")
        assert isinstance(mrna, SeqRecord)

    def test_mrna_sequence_is_exon_concatenation(self, fetcher: SequenceFetcher) -> None:
        seq = "AAAAACCCCC" * 60  # 600 bp
        record = _make_genbank_record(seq=seq, exon_ranges=[(0, 100), (200, 300)])
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_genbank_handle(record),
        ):
            mrna = fetcher.get_mrna_sequence("NM_001101")
        expected = seq[0:100] + seq[200:300]
        assert str(mrna.seq) == expected

    def test_mrna_id_matches_accession(self, fetcher: SequenceFetcher) -> None:
        record = _make_genbank_record(exon_ranges=[(0, 200)])
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_genbank_handle(record),
        ):
            mrna = fetcher.get_mrna_sequence("NM_001101")
        assert mrna.id == "NM_001101"


# ---------------------------------------------------------------------------
# Tests: in-memory cache (cache_sequence / load_from_cache)
# ---------------------------------------------------------------------------


class TestMemoryCache:
    """Test in-memory cache_sequence / load_from_cache."""

    def test_cache_and_load(self, fetcher: SequenceFetcher) -> None:
        record = _make_fasta_record()
        fetcher.cache_sequence("NM_001101", record)
        loaded = fetcher.load_from_cache("NM_001101")
        assert loaded is record

    def test_miss_returns_none(self, fetcher: SequenceFetcher) -> None:
        assert fetcher.load_from_cache("NM_NOTEXIST") is None

    def test_cache_overwrite(self, fetcher: SequenceFetcher) -> None:
        record1 = _make_fasta_record(acc="NM_001101", seq="AAAA")
        record2 = _make_fasta_record(acc="NM_001101", seq="CCCC")
        fetcher.cache_sequence("NM_001101", record1)
        fetcher.cache_sequence("NM_001101", record2)
        loaded = fetcher.load_from_cache("NM_001101")
        assert str(loaded.seq) == "CCCC"  # type: ignore[union-attr]


# ---------------------------------------------------------------------------
# Tests: SQLite cache
# ---------------------------------------------------------------------------


class TestSQLiteCache:
    """Test SQLite-backed cache persistence."""

    def test_sqlite_db_created(self, config: AppConfig, tmp_path: Path) -> None:
        with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
            SequenceFetcher(config, cache_dir=tmp_path / "cache")
        assert (tmp_path / "cache" / "sequences.db").exists()

    def test_store_and_reload_fasta(self, config: AppConfig, tmp_path: Path) -> None:
        cache_dir = tmp_path / "cache"
        record = _make_fasta_record()

        with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
            f1 = SequenceFetcher(config, cache_dir=cache_dir)
        f1._store_record("NM_001101", record, fmt="fasta")
        f1.close()

        # New fetcher instance, same cache dir — should hit SQLite.
        with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
            f2 = SequenceFetcher(config, cache_dir=cache_dir)
        loaded = f2._load_record("NM_001101")
        assert loaded is not None
        assert loaded.id == "NM_001101"
        f2.close()

    def test_store_and_reload_genbank(self, config: AppConfig, tmp_path: Path) -> None:
        cache_dir = tmp_path / "cache"
        record = _make_genbank_record(exon_ranges=[(0, 100), (200, 400)])

        with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
            f1 = SequenceFetcher(config, cache_dir=cache_dir)
        f1._store_record("NM_001101", record, fmt="genbank")
        f1.close()

        with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
            f2 = SequenceFetcher(config, cache_dir=cache_dir)
        loaded = f2._load_record("NM_001101")
        assert loaded is not None
        exon_feats = [feat for feat in loaded.features if feat.type == "exon"]
        assert len(exon_feats) == 2
        f2.close()

    def test_fetch_fasta_stores_in_sqlite(self, fetcher_with_cache: SequenceFetcher) -> None:
        record = _make_fasta_record()
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_fasta_handle(record),
        ):
            fetcher_with_cache.fetch_fasta("NM_001101")

        # Direct SQLite check
        assert fetcher_with_cache._db_conn is not None
        row = fetcher_with_cache._db_conn.execute(
            "SELECT accession FROM sequence_cache WHERE accession = ?",
            ("NM_001101",),
        ).fetchone()
        assert row is not None
        fetcher_with_cache.close()

    def test_fetch_uses_sqlite_cache_on_second_instance(
        self, config: AppConfig, tmp_path: Path
    ) -> None:
        cache_dir = tmp_path / "cache"
        record = _make_fasta_record()

        # First fetch — populates SQLite.
        with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
            f1 = SequenceFetcher(config, cache_dir=cache_dir)
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            return_value=_fasta_handle(record),
        ):
            f1.fetch_fasta("NM_001101")
        f1.close()

        # Second fetch — should hit SQLite without calling Entrez.
        with patch("tta_primer_design.modules.sequence_fetcher.setup_entrez"):
            f2 = SequenceFetcher(config, cache_dir=cache_dir)
        with patch("tta_primer_design.modules.sequence_fetcher.entrez_efetch") as mock_fetch:
            f2.fetch_fasta("NM_001101")
        mock_fetch.assert_not_called()
        f2.close()

    def test_close_releases_connection(self, fetcher_with_cache: SequenceFetcher) -> None:
        fetcher_with_cache.close()
        assert fetcher_with_cache._db_conn is None

    def test_close_idempotent(self, fetcher_with_cache: SequenceFetcher) -> None:
        fetcher_with_cache.close()
        fetcher_with_cache.close()  # second close must not raise


# ---------------------------------------------------------------------------
# Tests: retry behaviour
# ---------------------------------------------------------------------------


class TestRetryBehaviour:
    """Verify that RuntimeError propagates when all retries are exhausted."""

    def test_runtime_error_on_all_retries_exhausted(self, fetcher: SequenceFetcher) -> None:
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            side_effect=RuntimeError("all retries failed"),
        ):
            with pytest.raises(RuntimeError, match="all retries failed"):
                fetcher.fetch_fasta("NM_FAIL")

    def test_runtime_error_propagates_for_genbank(self, fetcher: SequenceFetcher) -> None:
        with patch(
            "tta_primer_design.modules.sequence_fetcher.entrez_efetch",
            side_effect=RuntimeError("network error"),
        ):
            with pytest.raises(RuntimeError):
                fetcher.fetch_genbank("NM_FAIL")
