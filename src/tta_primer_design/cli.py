"""CLI cho TTA Primer Design — dùng Click."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from tta_primer_design import __version__
from tta_primer_design.main import run_pipeline


@click.group()
@click.version_option(version=__version__, prog_name="tta-primer-design")
def main() -> None:
    """TTA Primer Design — Pipeline thiết kế primer/probe.

    Sử dụng lệnh ``run`` để bắt đầu pipeline:

        tta-primer-design run --config config/pipeline_config.yaml \\
            --input targets.csv --output results/
    """


@main.command("run")
@click.option(
    "--config",
    "-c",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    show_default=True,
    help="Đường dẫn file YAML config (tuỳ chọn — dùng default nếu bỏ qua).",
)
@click.option(
    "--input",
    "-i",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="File input: JSON, CSV, FASTA, hoặc danh sách accession (.txt).",
)
@click.option(
    "--output",
    "-o",
    "output_dir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Thư mục lưu kết quả.",
)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    help="Mức độ log.",
)
def run_command(
    config_path: Path | None,
    input_path: Path,
    output_dir: Path,
    log_level: str,
) -> None:
    """Chạy pipeline thiết kế primer/probe.

    Ví dụ:

    \b
        # PCR cơ bản
        tta-primer-design run -i targets.csv -o results/

        # qPCR + TaqMan với config tuỳ chỉnh
        tta-primer-design run -c config/pipeline_config.yaml \\
            -i targets.json -o results/run_001/
    """
    from tta_primer_design.logging_setup import setup_logging

    setup_logging(level=log_level)

    try:
        results = run_pipeline(
            config_path=config_path,
            input_path=input_path,
            output_dir=output_dir,
        )
        success = sum(1 for r in results if r.status == "success")
        click.echo(f"✅ Pipeline hoàn thành: {success}/{len(results)} target(s) thành công.")
        sys.exit(0 if success == len(results) else 1)
    except FileNotFoundError as exc:
        click.echo(f"❌ File không tìm thấy: {exc}", err=True)
        sys.exit(2)
    except ValueError as exc:
        click.echo(f"❌ Lỗi cấu hình: {exc}", err=True)
        sys.exit(2)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"❌ Lỗi không xác định: {exc}", err=True)
        sys.exit(3)
