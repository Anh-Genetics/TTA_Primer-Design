"""Thiết lập logging cho TTA Primer Design pipeline."""

from __future__ import annotations

import logging
import sys
from pathlib import Path


def setup_logging(
    level: str = "INFO",
    log_file: str | Path | None = None,
    fmt: str = "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt: str = "%Y-%m-%d %H:%M:%S",
) -> logging.Logger:
    """Cấu hình root logger cho pipeline.

    Args:
        level: Mức log — "DEBUG" | "INFO" | "WARNING" | "ERROR".
        log_file: Đường dẫn file log (tuỳ chọn). None = chỉ log ra console.
        fmt: Format string cho log messages.
        datefmt: Format cho timestamp.

    Returns:
        Root logger đã được cấu hình.

    Example::

        from tta_primer_design.logging_setup import setup_logging
        logger = setup_logging(level="DEBUG", log_file="logs/run.log")
        logger.info("Pipeline started")
    """
    numeric_level = getattr(logging, level.upper(), logging.INFO)

    handlers: list[logging.Handler] = [
        logging.StreamHandler(sys.stdout),
    ]

    if log_file is not None:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path, encoding="utf-8"))

    logging.basicConfig(
        level=numeric_level,
        format=fmt,
        datefmt=datefmt,
        handlers=handlers,
        force=True,  # override bất kỳ config logging nào trước đó
    )

    logger = logging.getLogger("tta_primer_design")
    logger.setLevel(numeric_level)
    return logger


def get_logger(name: str) -> logging.Logger:
    """Lấy logger con theo tên module.

    Args:
        name: Tên module (thường dùng ``__name__``).

    Returns:
        Logger đã được namespace hoá dưới "tta_primer_design".

    Example::

        logger = get_logger(__name__)
        logger.debug("Processing target: %s", target_id)
    """
    return logging.getLogger(f"tta_primer_design.{name}")
