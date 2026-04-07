"""Tests cho InputParser — parse JSON, CSV, FASTA, TXT."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.models import DesignTarget
from tta_primer_design.modules.input_parser import InputParser

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def default_config() -> AppConfig:
    return load_config(None)


@pytest.fixture
def test_config() -> AppConfig:
    return load_config(FIXTURES / "test_config.yaml")


@pytest.fixture
def parser(default_config: AppConfig) -> InputParser:
    return InputParser(default_config)


class TestParseJSON:
    """Test parse JSON input file."""

    def test_parse_fixture_json(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        assert len(targets) == 2

    def test_targets_are_design_target(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        for t in targets:
            assert isinstance(t, DesignTarget)

    def test_first_target_id(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        assert targets[0].target_id == "ACTB"

    def test_accession_parsed(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        assert targets[0].accession == "NM_001101"

    def test_design_mode_parsed(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        assert targets[0].design_mode == "qpcr"

    def test_exon_junction_parsed(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        assert targets[0].exon_junction is True

    def test_organism_parsed(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        assert targets[0].organism == "Homo sapiens"

    def test_input_type_is_accession(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_input.json")
        assert targets[0].input_type == "accession"

    def test_inline_json(self, parser: InputParser, tmp_path: Path) -> None:
        data = {"targets": [{"id": "GENE1", "accession": "NM_000001", "design_mode": "pcr"}]}
        path = tmp_path / "inline.json"
        path.write_text(json.dumps(data))
        targets = parser.parse(path)
        assert len(targets) == 1
        assert targets[0].target_id == "GENE1"

    def test_empty_targets_list(self, parser: InputParser, tmp_path: Path) -> None:
        data = {"targets": []}
        path = tmp_path / "empty.json"
        path.write_text(json.dumps(data))
        targets = parser.parse(path)
        assert targets == []

    def test_sequence_input_type(self, parser: InputParser, tmp_path: Path) -> None:
        data = {"targets": [{"id": "SEQ1", "sequence": "ATCGATCGATCG", "design_mode": "pcr"}]}
        path = tmp_path / "seq.json"
        path.write_text(json.dumps(data))
        targets = parser.parse(path)
        assert targets[0].input_type == "sequence"
        assert targets[0].sequence == "ATCGATCGATCG"


class TestParseCSV:
    """Test parse CSV input file."""

    def test_parse_fixture_csv(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_targets.csv")
        assert len(targets) == 3

    def test_csv_ids(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_targets.csv")
        ids = [t.target_id for t in targets]
        assert "ACTB" in ids
        assert "GAPDH" in ids
        assert "TP53" in ids

    def test_csv_organism(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_targets.csv")
        for t in targets:
            assert t.organism == "Homo sapiens"

    def test_csv_exon_junction_bool(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_targets.csv")
        actb = next(t for t in targets if t.target_id == "ACTB")
        gapdh = next(t for t in targets if t.target_id == "GAPDH")
        assert actb.exon_junction is True
        assert gapdh.exon_junction is False

    def test_csv_missing_id_raises(self, parser: InputParser, tmp_path: Path) -> None:
        bad_csv = tmp_path / "bad.csv"
        bad_csv.write_text("accession,organism\nNM_001101,Homo sapiens\n")
        with pytest.raises(ValueError, match="id"):
            parser.parse(bad_csv)


class TestParseFASTA:
    """Test parse FASTA input file."""

    def test_parse_fixture_fasta(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_fasta.fa")
        assert len(targets) == 2

    def test_fasta_ids(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_fasta.fa")
        ids = [t.target_id for t in targets]
        assert "ACTB_fragment" in ids
        assert "TP53_fragment" in ids

    def test_fasta_input_type(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_fasta.fa")
        for t in targets:
            assert t.input_type == "sequence"

    def test_fasta_sequence_not_empty(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_fasta.fa")
        for t in targets:
            assert t.sequence and len(t.sequence) > 0

    def test_fasta_sequence_uppercase(self, parser: InputParser) -> None:
        targets = parser.parse(FIXTURES / "sample_fasta.fa")
        for t in targets:
            assert t.sequence == t.sequence.upper()


class TestParseTXT:
    """Test parse accession list TXT file."""

    def test_parse_accession_list(self, parser: InputParser, tmp_path: Path) -> None:
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("NM_001101\nNM_000546\nNM_002046\n")
        targets = parser.parse(acc_file)
        assert len(targets) == 3

    def test_accession_ids(self, parser: InputParser, tmp_path: Path) -> None:
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("NM_001101\nNM_000546\n")
        targets = parser.parse(acc_file)
        assert targets[0].target_id == "NM_001101"
        assert targets[0].accession == "NM_001101"

    def test_txt_skips_blank_lines(self, parser: InputParser, tmp_path: Path) -> None:
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("NM_001101\n\nNM_000546\n\n")
        targets = parser.parse(acc_file)
        assert len(targets) == 2

    def test_txt_skips_comment_lines(self, parser: InputParser, tmp_path: Path) -> None:
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("# Human genes\nNM_001101\n# another comment\nNM_000546\n")
        targets = parser.parse(acc_file)
        assert len(targets) == 2

    def test_txt_input_type_is_accession(self, parser: InputParser, tmp_path: Path) -> None:
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("NM_001101\n")
        targets = parser.parse(acc_file)
        assert targets[0].input_type == "accession"


class TestParseErrors:
    """Test parse error cases."""

    def test_file_not_found(self, parser: InputParser) -> None:
        with pytest.raises(FileNotFoundError):
            parser.parse("/nonexistent/file.json")

    def test_unsupported_extension(self, parser: InputParser, tmp_path: Path) -> None:
        bad_file = tmp_path / "data.xyz"
        bad_file.write_text("data")
        with pytest.raises(ValueError, match="Unsupported"):
            parser.parse(bad_file)
