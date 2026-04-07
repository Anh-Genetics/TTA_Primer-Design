"""Entry point chính của pipeline TTA Primer Design."""

from __future__ import annotations

import logging
from pathlib import Path

from tta_primer_design.config import AppConfig, load_config
from tta_primer_design.logging_setup import setup_logging
from tta_primer_design.models import DesignResult
from tta_primer_design.modules.input_parser import InputParser

logger = logging.getLogger("tta_primer_design.main")


def run_pipeline(
    config_path: str | Path | None,
    input_path: str | Path,
    output_dir: str | Path,
) -> list[DesignResult]:
    """Chạy toàn bộ pipeline thiết kế primer/probe.

    Luồng xử lý:
        1. Load config
        2. Parse input → List[DesignTarget]
        3. Với mỗi target:
           a. Fetch sequence (SequenceFetcher)
           b. Preprocess (SequencePreprocessor)
           c. Thiết kế primer (Primer3Runner hoặc NCBIPrimerBlast)
           d. BLAST specificity check
           e. Thiết kế probe nếu mode = taqman/qpcr
           f. SNP check
           g. Filter & rank
        4. Tạo báo cáo (ReportGenerator)

    Args:
        config_path: Đường dẫn tới file YAML config (None = dùng default).
        input_path: Đường dẫn file input (JSON/CSV/FASTA/txt).
        output_dir: Thư mục để lưu kết quả.

    Returns:
        Danh sách DesignResult cho mỗi target.
    """
    # --- 1. Load config & setup logging ---
    cfg: AppConfig = load_config(config_path)
    setup_logging(log_file=Path(output_dir) / "pipeline.log")
    logger.info("Pipeline '%s' v%s started", cfg.pipeline.name, cfg.pipeline.version)
    logger.info("Mode: %s | Local Primer3: %s", cfg.pipeline.mode, cfg.pipeline.use_local_primer3)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # --- 2. Parse input ---
    parser = InputParser(cfg)
    targets = parser.parse(input_path)
    logger.info("Loaded %d design target(s)", len(targets))

    results: list[DesignResult] = []

    for target in targets:
        logger.info("Processing target: %s", target.target_id)
        try:
            result = _process_target(target, cfg)
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed target '%s': %s", target.target_id, exc, exc_info=True)
            result = DesignResult(target=target, status="failed", error=str(exc))

        results.append(result)

    # --- 4. Report (stub) ---
    logger.info(
        "Pipeline completed. %d/%d targets succeeded.",
        sum(1 for r in results if r.status == "success"),
        len(results),
    )
    return results


def _process_target(target, cfg: AppConfig) -> DesignResult:
    """Xử lý một DesignTarget qua toàn bộ pipeline."""
    from tta_primer_design.modules.filter_ranker import FilterRanker
    from tta_primer_design.modules.primer3_runner import Primer3Runner
    from tta_primer_design.modules.probe_designer import ProbeDesigner
    from tta_primer_design.modules.sequence_fetcher import SequenceFetcher
    from tta_primer_design.modules.sequence_preprocessor import SequencePreprocessor

    preprocessor = SequencePreprocessor(cfg)
    primer3_runner = Primer3Runner(cfg)
    filter_ranker = FilterRanker(cfg)

    # Step 1: Get sequence
    seq_input = None
    if target.input_type == "sequence" and target.sequence:
        seq_input = target.sequence
    elif target.input_type in ("accession", "gene_name"):
        fetcher = SequenceFetcher(cfg)
        acc = target.accession or target.gene_name
        try:
            seq_input = fetcher.fetch_fasta(acc)
        except Exception as exc:
            logger.error("Failed to fetch sequence for '%s': %s", target.target_id, exc)
            return DesignResult(target=target, status="failed", error=str(exc))

    if seq_input is None:
        return DesignResult(target=target, status="failed", error="No sequence available")

    # Step 2: Preprocess
    try:
        processed_seq = preprocessor.process(seq_input, target)
    except ValueError as exc:
        return DesignResult(
            target=target, status="failed", error=f"Sequence validation failed: {exc}"
        )

    # Step 3: Design primers
    if cfg.pipeline.use_local_primer3:
        try:
            pairs = primer3_runner.run(processed_seq, target)
        except Exception as exc:
            logger.error("Primer3 failed for '%s': %s", target.target_id, exc)
            return DesignResult(target=target, status="failed", error=str(exc))
    else:
        from tta_primer_design.modules.ncbi_primer_blast import NCBIPrimerBlast

        api = NCBIPrimerBlast(cfg)
        try:
            pairs = api.design_primers(target, processed_seq)
        except NotImplementedError:
            logger.warning("NCBIPrimerBlast not implemented — falling back to local Primer3")
            try:
                pairs = primer3_runner.run(processed_seq, target)
            except Exception as exc:
                return DesignResult(target=target, status="failed", error=str(exc))
        except Exception as exc:
            return DesignResult(target=target, status="failed", error=str(exc))

    if not pairs:
        return DesignResult(target=target, primer_pairs=[], status="no_primers")

    # Step 4: Design probes if mode requires it
    if cfg.pipeline.mode in ("taqman", "qpcr"):
        probe_designer = ProbeDesigner(cfg)
        try:
            pairs = probe_designer.design(pairs, processed_seq)
        except Exception as exc:
            logger.warning("Probe design failed for '%s': %s", target.target_id, exc)

    # Step 5: Filter & rank
    final_pairs = filter_ranker.process(pairs)

    status = "success" if final_pairs else "no_primers"
    return DesignResult(target=target, primer_pairs=final_pairs, status=status)
