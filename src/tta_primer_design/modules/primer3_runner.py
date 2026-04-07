"""Module 4 — Primer3Runner: chạy Primer3 để thiết kế primer.

Sử dụng ``primer3-py`` (Python binding cho Primer3 C library).
Tài liệu: https://primer3-py.readthedocs.io/

Nếu ``primer3-py`` chưa cài, raise ImportError rõ ràng với hướng dẫn.

Luồng:
    1. Load tham số từ YAML config (primer3_qpcr_params.yaml, v.v.)
    2. Merge với custom_params từ DesignTarget
    3. Gọi ``primer3.bindings.design_primers()``
    4. Parse output thành List[PrimerPair]
    5. Xử lý lỗi: không tìm được primer, constraint quá chặt → auto-relax
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget, Oligo, PrimerPair, ProcessedSequence

logger = logging.getLogger("tta_primer_design.modules.primer3_runner")

# Map từ pipeline.mode → file config primer3 tương ứng
_MODE_TO_CONFIG_FILE: dict[str, str] = {
    "qpcr": "primer3_qpcr_params.yaml",
    "taqman": "primer3_qpcr_params.yaml",
    "pcr": "primer3_pcr_params.yaml",
    "sybr": "primer3_pcr_params.yaml",
    "multiplex": "primer3_pcr_params.yaml",
}

# Các tham số được nới lỏng khi không tìm được primer (auto-relax)
_RELAX_STEPS: list[dict[str, Any]] = [
    {"PRIMER_MIN_TM": 2.0, "PRIMER_MAX_TM": 2.0, "PRIMER_PAIR_MAX_DIFF_TM": 1.0},
    {"PRIMER_MIN_SIZE": -2, "PRIMER_MAX_SIZE": 2},
    {"PRIMER_MIN_GC": -5.0, "PRIMER_MAX_GC": 5.0},
]


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

        Luồng:
            1. Load global params từ YAML.
            2. Merge với target.custom_params (nếu có).
            3. Build seq_args từ processed_seq.
            4. Gọi primer3.design_primers().
            5. Nếu không có kết quả → thử auto-relax constraints (tối đa len(_RELAX_STEPS) lần).
            6. Parse và trả về danh sách PrimerPair.

        Args:
            processed_seq: Sequence đã tiền xử lý.
            target: DesignTarget (chứa custom_params nếu có).

        Returns:
            Danh sách PrimerPair từ Primer3 (có thể rỗng nếu không tìm được).

        Raises:
            ImportError: Nếu primer3-py chưa được cài.
        """
        try:
            import primer3
        except ImportError as exc:
            raise ImportError("primer3-py chưa được cài. Chạy: pip install primer3-py") from exc

        global_args = self.load_primer3_params()

        # Merge với custom_params của target
        if target.custom_params:
            global_args.update(target.custom_params)

        # Giới hạn số lượng pair trả về theo config
        global_args.setdefault("PRIMER_NUM_RETURN", self.config.pipeline.top_n_pairs)

        seq_args = self._build_seq_args(processed_seq, target)

        try:
            result = primer3.design_primers(seq_args=seq_args, global_args=global_args)
        except OSError as exc:
            logger.warning(
                "Primer3 lỗi cho target '%s': %s — trả về danh sách rỗng.",
                target.target_id,
                exc,
            )
            return []
        num_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)

        # Auto-relax nếu không tìm được primer
        if num_returned == 0:
            relaxed_args = dict(global_args)
            for step in _RELAX_STEPS:
                relaxed_args = _apply_relax_step(relaxed_args, step)
                try:
                    result = primer3.design_primers(seq_args=seq_args, global_args=relaxed_args)
                except OSError as exc:
                    logger.warning(
                        "Primer3 lỗi khi auto-relax cho target '%s': %s",
                        target.target_id,
                        exc,
                    )
                    break
                num_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
                if num_returned > 0:
                    logger.info(
                        "Auto-relax thành công sau khi nới lỏng constraints cho target '%s': "
                        "%d pairs",
                        target.target_id,
                        num_returned,
                    )
                    break
            else:
                logger.warning(
                    "Không tìm được primer cho target '%s' sau khi auto-relax.",
                    target.target_id,
                )

        pairs = self._parse_primer3_output(result, target.target_id)
        logger.debug("Primer3 trả về %d pair(s) cho target '%s'", len(pairs), target.target_id)
        return pairs

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
            Danh sách PrimerPair (rỗng nếu PRIMER_PAIR_NUM_RETURNED == 0).
        """
        num_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if not num_returned:
            return []

        pairs: list[PrimerPair] = []
        for i in range(num_returned):
            left = _parse_oligo(result, "LEFT", i)
            right = _parse_oligo(result, "RIGHT", i)
            if left is None or right is None:
                continue

            probe: Oligo | None = None
            if f"PRIMER_INTERNAL_{i}_SEQUENCE" in result or f"PRIMER_INTERNAL_{i}" in result:
                probe = _parse_oligo(result, "INTERNAL", i)

            amplicon_size = result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)
            pair_penalty = result.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0)

            pair = PrimerPair(
                pair_id=f"{target_id}_pair_{i}",
                left_primer=left,
                right_primer=right,
                probe=probe,
                amplicon_size=int(amplicon_size),
                pair_penalty=float(pair_penalty),
            )
            pairs.append(pair)

        return pairs


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _parse_oligo(result: dict[str, Any], side: str, idx: int) -> Oligo | None:
    """Parse một oligo từ kết quả primer3-py.

    Args:
        result: Dict kết quả từ primer3.
        side: "LEFT", "RIGHT", hoặc "INTERNAL".
        idx: Index của pair.

    Returns:
        Oligo hoặc None nếu không có dữ liệu.
    """
    prefix = f"PRIMER_{side}_{idx}"
    seq_key = f"{prefix}_SEQUENCE"
    pos_key = f"PRIMER_{side}_{idx}"  # same key, value = [start, length]

    sequence = result.get(seq_key, "")
    if not sequence:
        return None

    position = result.get(pos_key, [0, len(sequence)])
    start = int(position[0]) if isinstance(position, (list, tuple)) else 0

    return Oligo(
        sequence=str(sequence),
        start=start,
        tm=float(result.get(f"{prefix}_TM", 0.0)),
        gc_percent=float(result.get(f"{prefix}_GC_PERCENT", 0.0)),
        self_any_th=float(result.get(f"{prefix}_SELF_ANY_TH", 0.0)),
        self_end_th=float(result.get(f"{prefix}_SELF_END_TH", 0.0)),
        hairpin_th=float(result.get(f"{prefix}_HAIRPIN_TH", 0.0)),
        end_stability=float(result.get(f"{prefix}_END_STABILITY", 0.0)),
        penalty=float(result.get(f"{prefix}_PENALTY", 0.0)),
    )


def _apply_relax_step(params: dict[str, Any], step: dict[str, Any]) -> dict[str, Any]:
    """Nới lỏng constraints theo một bước relax.

    Giá trị dương trong step → trừ vào MIN / cộng vào MAX.
    Tên key dùng quy ước: MIN_* giảm, MAX_* tăng, DIFF_* tăng.
    """
    relaxed = dict(params)
    for key, delta in step.items():
        if key.startswith("PRIMER_MIN_") or key.startswith("PRIMER_INTERNAL_MIN_"):
            if key in relaxed:
                relaxed[key] = relaxed[key] - abs(delta)
        elif key.startswith("PRIMER_MAX_") or key.startswith("PRIMER_INTERNAL_MAX_"):
            if key in relaxed:
                relaxed[key] = relaxed[key] + abs(delta)
        elif "DIFF" in key or "PAIR_MAX" in key:
            if key in relaxed:
                relaxed[key] = relaxed[key] + abs(delta)
    return relaxed
