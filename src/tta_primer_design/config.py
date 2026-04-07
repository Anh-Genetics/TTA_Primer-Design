"""Load và validate cấu hình pipeline từ file YAML."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# Sub-config dataclasses
# ---------------------------------------------------------------------------


@dataclass
class NCBIConfig:
    email: str = "user@example.com"
    api_key: str = ""
    retries: int = 3
    timeout: int = 120
    rate_limit: int = 3


@dataclass
class BLASTConfig:
    database: str = "refseq_rna"
    organism: str = "Homo sapiens"
    max_targets: int = 500
    word_size: int = 7
    evalue: float = 1000.0
    num_threads: int = 4


@dataclass
class OutputConfig:
    formats: list[str] = field(default_factory=lambda: ["csv", "json"])
    include_plots: bool = False


@dataclass
class FiltersConfig:
    min_specificity_score: float = 80.0
    max_off_targets: int = 0
    require_exon_junction: bool = False
    avoid_snp_in_primer: bool = True
    avoid_snp_in_probe: bool = True


@dataclass
class PipelineConfig:
    name: str = "Primer_Probe_Design_Pipeline"
    version: str = "0.1.0"
    mode: str = "qpcr"
    use_local_primer3: bool = True
    top_n_pairs: int = 5


# ---------------------------------------------------------------------------
# Top-level config
# ---------------------------------------------------------------------------


@dataclass
class AppConfig:
    """Cấu hình toàn bộ ứng dụng.

    Được load từ YAML và validate tại runtime.
    """

    pipeline: PipelineConfig = field(default_factory=PipelineConfig)
    ncbi: NCBIConfig = field(default_factory=NCBIConfig)
    blast: BLASTConfig = field(default_factory=BLASTConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    filters: FiltersConfig = field(default_factory=FiltersConfig)


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

_VALID_MODES = {"pcr", "qpcr", "taqman", "sybr", "multiplex"}


def load_config(config_path: str | Path | None = None) -> AppConfig:
    """Load AppConfig từ file YAML.

    Args:
        config_path: Đường dẫn tới file YAML (tuỳ chọn).
                     Nếu None, trả về config mặc định.

    Returns:
        AppConfig đã được populate.

    Raises:
        FileNotFoundError: Nếu file không tồn tại.
        ValueError: Nếu giá trị config không hợp lệ.
    """
    if config_path is None:
        return AppConfig()

    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r", encoding="utf-8") as fh:
        raw: dict = yaml.safe_load(fh) or {}

    cfg = AppConfig()

    # --- pipeline ---
    if pipeline_raw := raw.get("pipeline"):
        cfg.pipeline = PipelineConfig(
            name=pipeline_raw.get("name", cfg.pipeline.name),
            version=pipeline_raw.get("version", cfg.pipeline.version),
            mode=pipeline_raw.get("mode", cfg.pipeline.mode),
            use_local_primer3=pipeline_raw.get("use_local_primer3", cfg.pipeline.use_local_primer3),
            top_n_pairs=pipeline_raw.get("top_n_pairs", cfg.pipeline.top_n_pairs),
        )

    # --- ncbi ---
    if ncbi_raw := raw.get("ncbi"):
        cfg.ncbi = NCBIConfig(
            email=ncbi_raw.get("email", cfg.ncbi.email),
            api_key=ncbi_raw.get("api_key", cfg.ncbi.api_key),
            retries=ncbi_raw.get("retries", cfg.ncbi.retries),
            timeout=ncbi_raw.get("timeout", cfg.ncbi.timeout),
            rate_limit=ncbi_raw.get("rate_limit", cfg.ncbi.rate_limit),
        )

    # --- blast ---
    if blast_raw := raw.get("blast"):
        cfg.blast = BLASTConfig(
            database=blast_raw.get("database", cfg.blast.database),
            organism=blast_raw.get("organism", cfg.blast.organism),
            max_targets=blast_raw.get("max_targets", cfg.blast.max_targets),
            word_size=blast_raw.get("word_size", cfg.blast.word_size),
            evalue=blast_raw.get("evalue", cfg.blast.evalue),
            num_threads=blast_raw.get("num_threads", cfg.blast.num_threads),
        )

    # --- output ---
    if output_raw := raw.get("output"):
        cfg.output = OutputConfig(
            formats=output_raw.get("formats", cfg.output.formats),
            include_plots=output_raw.get("include_plots", cfg.output.include_plots),
        )

    # --- filters ---
    if filters_raw := raw.get("filters"):
        cfg.filters = FiltersConfig(
            min_specificity_score=filters_raw.get(
                "min_specificity_score", cfg.filters.min_specificity_score
            ),
            max_off_targets=filters_raw.get("max_off_targets", cfg.filters.max_off_targets),
            require_exon_junction=filters_raw.get(
                "require_exon_junction", cfg.filters.require_exon_junction
            ),
            avoid_snp_in_primer=filters_raw.get(
                "avoid_snp_in_primer", cfg.filters.avoid_snp_in_primer
            ),
            avoid_snp_in_probe=filters_raw.get(
                "avoid_snp_in_probe", cfg.filters.avoid_snp_in_probe
            ),
        )

    # Validate
    _validate(cfg)
    return cfg


def _validate(cfg: AppConfig) -> None:
    """Kiểm tra tính hợp lệ của config.

    Raises:
        ValueError: Nếu có trường không hợp lệ.
    """
    if cfg.pipeline.mode not in _VALID_MODES:
        raise ValueError(
            f"Invalid pipeline.mode '{cfg.pipeline.mode}'. "
            f"Must be one of: {sorted(_VALID_MODES)}"
        )
    if cfg.pipeline.top_n_pairs < 1:
        raise ValueError("pipeline.top_n_pairs must be >= 1")
    if cfg.ncbi.rate_limit < 1:
        raise ValueError("ncbi.rate_limit must be >= 1")
    if cfg.filters.min_specificity_score < 0 or cfg.filters.min_specificity_score > 100:
        raise ValueError("filters.min_specificity_score must be between 0 and 100")
