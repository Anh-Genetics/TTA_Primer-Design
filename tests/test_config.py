"""Tests cho config.py — load và validate AppConfig từ YAML."""

from __future__ import annotations

from pathlib import Path

import pytest

from tta_primer_design.config import AppConfig, load_config

FIXTURES = Path(__file__).parent / "fixtures"


class TestLoadConfigDefault:
    """Test load_config với config mặc định (không có file)."""

    def test_returns_app_config(self) -> None:
        cfg = load_config(None)
        assert isinstance(cfg, AppConfig)

    def test_default_mode_is_qpcr(self) -> None:
        cfg = load_config(None)
        assert cfg.pipeline.mode == "qpcr"

    def test_default_use_local_primer3(self) -> None:
        cfg = load_config(None)
        assert cfg.pipeline.use_local_primer3 is True

    def test_default_top_n_pairs(self) -> None:
        cfg = load_config(None)
        assert cfg.pipeline.top_n_pairs == 5

    def test_default_ncbi_email(self) -> None:
        cfg = load_config(None)
        assert "@" in cfg.ncbi.email

    def test_default_blast_database(self) -> None:
        cfg = load_config(None)
        assert cfg.blast.database == "refseq_rna"


class TestLoadConfigFromFile:
    """Test load_config từ fixture YAML."""

    def test_load_test_config(self) -> None:
        cfg = load_config(FIXTURES / "test_config.yaml")
        assert isinstance(cfg, AppConfig)

    def test_name_from_file(self) -> None:
        cfg = load_config(FIXTURES / "test_config.yaml")
        assert cfg.pipeline.name == "Test_Pipeline"

    def test_top_n_from_file(self) -> None:
        cfg = load_config(FIXTURES / "test_config.yaml")
        assert cfg.pipeline.top_n_pairs == 3

    def test_ncbi_email_from_file(self) -> None:
        cfg = load_config(FIXTURES / "test_config.yaml")
        assert cfg.ncbi.email == "test@example.com"

    def test_output_format_from_file(self) -> None:
        cfg = load_config(FIXTURES / "test_config.yaml")
        assert "json" in cfg.output.formats

    def test_filters_from_file(self) -> None:
        cfg = load_config(FIXTURES / "test_config.yaml")
        assert cfg.filters.min_specificity_score == 80.0
        assert cfg.filters.avoid_snp_in_primer is True


class TestLoadConfigErrors:
    """Test load_config với input không hợp lệ."""

    def test_file_not_found(self) -> None:
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/path/config.yaml")

    def test_invalid_mode_raises_value_error(self, tmp_path: Path) -> None:
        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text("pipeline:\n  mode: invalid_mode\n")
        with pytest.raises(ValueError, match="invalid_mode"):
            load_config(bad_config)

    def test_invalid_top_n_raises_value_error(self, tmp_path: Path) -> None:
        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text("pipeline:\n  top_n_pairs: 0\n")
        with pytest.raises(ValueError, match="top_n_pairs"):
            load_config(bad_config)

    def test_invalid_specificity_score_raises_value_error(self, tmp_path: Path) -> None:
        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text("filters:\n  min_specificity_score: 150\n")
        with pytest.raises(ValueError, match="min_specificity_score"):
            load_config(bad_config)


class TestLoadConfigAllModes:
    """Test tất cả pipeline modes hợp lệ."""

    @pytest.mark.parametrize("mode", ["pcr", "qpcr", "taqman", "sybr", "multiplex"])
    def test_valid_modes(self, tmp_path: Path, mode: str) -> None:
        config_file = tmp_path / "cfg.yaml"
        config_file.write_text(f"pipeline:\n  mode: {mode}\n")
        cfg = load_config(config_file)
        assert cfg.pipeline.mode == mode
