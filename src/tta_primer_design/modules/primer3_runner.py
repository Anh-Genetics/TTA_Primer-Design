"""Module 4 — Primer3Runner: chạy Primer3 để thiết kế primer.

Sử dụng ``primer3-py`` (Python binding cho Primer3 C library).
Tài liệu: https://primer3-py.readthedocs.io/

Nếu ``primer3-py`` chưa cài, raise ImportError rõ ràng với hướng dẫn.

Luồng:
    1. Load tham số từ YAML config (primer3_qpcr_params.yaml, v.v.)
    2. Merge với custom_params từ DesignTarget
    3. Gọi ``primer3.bindings.design_primers()``
    4. Parse output thành List[PrimerPair]
    5. Xử lý lỗi: không tìm được primer, constraint quá chặt

TODO (Sprint 2):
    - Implement ``run()``
    - Implement ``_parse_primer3_output()``
    - Implement auto-relax constraints nếu không tìm được primer
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget, PrimerPair, ProcessedSequence

logger = logging.getLogger("tta_primer_design.modules.primer3_runner")

# Map từ pipeline.mode → file config primer3 tương ứng
_MODE_TO_CONFIG_FILE: dict[str, str] = {
    "qpcr": "primer3_qpcr_params.yaml",
    "taqman": "primer3_qpcr_params.yaml",
    "pcr": "primer3_pcr_params.yaml",
    "sybr": "primer3_pcr_params.yaml",
    "multiplex": "primer3_pcr_params.yaml",
}


class Primer3Runner:
    """Chạy Primer3 thông qua primer3-py.

    Args:
        config: AppConfig đã load.
        config_dir: Thư mục chứa các file YAML tham số Primer3.

    Example::

        runner = Primer3Runner(config, config_dir="config/")
        pairs = runner.run(processed_seq, target)
    """

    def __init__(self, config: AppConfig, config_dir: str | Path = "config") -> None:
        self.config = config
        self.config_dir = Path(config_dir)
        self._p3_params: dict[str, Any] = {}

    def load_primer3_params(self) -> dict[str, Any]:
        """Load tham số Primer3 từ YAML tương ứng với pipeline.mode.

        Returns:
            Dict tham số Primer3.

        Raises:
            FileNotFoundError: Nếu file YAML không tồn tại.
        """
        mode = self.config.pipeline.mode
        filename = _MODE_TO_CONFIG_FILE.get(mode, "primer3_qpcr_params.yaml")
        path = self.config_dir / filename

        if not path.exists():
            logger.warning("Primer3 config file not found: %s — using empty params", path)
            return {}

        with path.open("r", encoding="utf-8") as fh:
            params: dict[str, Any] = yaml.safe_load(fh) or {}

        logger.debug("Loaded %d Primer3 parameters from %s", len(params), filename)
        return params

    def run(
        self,
        processed_seq: ProcessedSequence,
        target: DesignTarget,
    ) -> list[PrimerPair]:
        """Chạy Primer3 và trả về danh sách PrimerPair.

        Args:
            processed_seq: Sequence đã tiền xử lý.
            target: DesignTarget (chứa custom_params nếu có).

        Returns:
            Danh sách PrimerPair từ Primer3.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
            ImportError: Nếu primer3-py chưa được cài.
        """
        try:
            import primer3  # noqa: F401
        except ImportError as exc:
            raise ImportError("primer3-py chưa được cài. Chạy: pip install primer3-py") from exc

        raise NotImplementedError("Primer3Runner.run() chưa được implement (Sprint 2)")

    def _build_seq_args(
        self,
        processed_seq: ProcessedSequence,
        target: DesignTarget,
    ) -> dict[str, Any]:
        """Tạo dict sequence arguments cho primer3.

        Args:
            processed_seq: ProcessedSequence đã xử lý.
            target: DesignTarget.

        Returns:
            Dict ``SEQUENCE_*`` arguments cho primer3-py.
        """
        seq_args: dict[str, Any] = {
            "SEQUENCE_ID": target.target_id,
            "SEQUENCE_TEMPLATE": processed_seq.sequence,
        }
        if processed_seq.included_region:
            start, length = processed_seq.included_region
            seq_args["SEQUENCE_INCLUDED_REGION"] = [start, length]
        if processed_seq.excluded_regions:
            seq_args["SEQUENCE_EXCLUDED_REGION"] = [
                [s, e - s] for s, e in processed_seq.excluded_regions
            ]
        if processed_seq.target_regions:
            seq_args["SEQUENCE_TARGET"] = [[s, e - s] for s, e in processed_seq.target_regions]
        return seq_args

    def _parse_primer3_output(self, result: dict[str, Any], target_id: str) -> list[PrimerPair]:
        """Parse output từ primer3-py thành List[PrimerPair].

        Args:
            result: Dict trả về từ ``primer3.bindings.design_primers()``.
            target_id: ID target để đặt tên pair.

        Returns:
            Danh sách PrimerPair.

        Raises:
            NotImplementedError: Chưa implement (Sprint 2).
        """
        raise NotImplementedError("_parse_primer3_output chưa được implement (Sprint 2)")
