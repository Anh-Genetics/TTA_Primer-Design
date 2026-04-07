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


@main.command("evaluate")
@click.option(
    "--left-primer",
    "-l",
    "left_seq",
    required=True,
    help="Trình tự primer xuôi (5'→3', ký tự IUPAC DNA).",
)
@click.option(
    "--right-primer",
    "-r",
    "right_seq",
    required=True,
    help="Trình tự primer ngược (5'→3', ký tự IUPAC DNA).",
)
@click.option(
    "--probe",
    "-p",
    "probe_seq",
    default=None,
    show_default=True,
    help="Trình tự probe (tuỳ chọn).",
)
@click.option(
    "--pair-id",
    default="user_pair",
    show_default=True,
    help="Tên/ID cho cặp primer.",
)
@click.option(
    "--organism",
    "-org",
    default="Homo sapiens",
    show_default=True,
    help="Tên loài NCBI dùng để lọc BLAST (ví dụ: 'Homo sapiens').",
)
@click.option(
    "--output",
    "-o",
    "output_dir",
    default=None,
    type=click.Path(file_okay=False, path_type=Path),
    help="Thư mục lưu báo cáo JSON (tuỳ chọn).",
)
@click.option(
    "--no-blast",
    "skip_blast",
    is_flag=True,
    default=False,
    help="Bỏ qua bước BLAST (chỉ tính nhiệt động học).",
)
@click.option(
    "--max-amplicon-size",
    "max_amplicon_size",
    default=4000,
    show_default=True,
    type=int,
    help="Kích thước amplicon off-target tối đa cần flag (bp).",
)
@click.option(
    "--config",
    "-c",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Đường dẫn file YAML config (tuỳ chọn).",
)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    help="Mức độ log.",
)
def evaluate_command(
    left_seq: str,
    right_seq: str,
    probe_seq: str | None,
    pair_id: str,
    organism: str,
    output_dir: Path | None,
    skip_blast: bool,
    max_amplicon_size: int,
    config_path: Path | None,
    log_level: str,
) -> None:
    """Đánh giá một cặp mồi sẵn có: nhiệt động học + BLAST specificity.

    Ví dụ:

    \b
        # Chỉ tính nhiệt động học (không cần internet)
        tta-primer-design evaluate \\
            --left-primer GCACTGACCTCCCACTTCAA \\
            --right-primer TTGCTGATCCACATCTGCTG \\
            --no-blast

        # Đánh giá đầy đủ với BLAST (cần internet)
        tta-primer-design evaluate \\
            --left-primer GCACTGACCTCCCACTTCAA \\
            --right-primer TTGCTGATCCACATCTGCTG \\
            --organism "Homo sapiens" \\
            --output results/
    """
    import json
    from dataclasses import asdict

    from tta_primer_design.config import load_config
    from tta_primer_design.logging_setup import setup_logging
    from tta_primer_design.modules.primer_evaluator import PrimerEvaluator

    setup_logging(level=log_level)

    try:
        cfg = load_config(config_path)
    except (FileNotFoundError, ValueError) as exc:
        click.echo(f"❌ Lỗi config: {exc}", err=True)
        sys.exit(2)

    evaluator = PrimerEvaluator(cfg, max_amplicon_size=max_amplicon_size)

    try:
        report = evaluator.evaluate(
            pair_id=pair_id,
            left_seq=left_seq,
            right_seq=right_seq,
            probe_seq=probe_seq,
            organism=organism,
            run_blast=not skip_blast,
        )
    except (ValueError, ImportError) as exc:
        click.echo(f"❌ Lỗi đánh giá: {exc}", err=True)
        sys.exit(2)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"❌ Lỗi không xác định: {exc}", err=True)
        sys.exit(3)

    # --- Print summary ---
    icon = {"PASS": "✅", "WARNING": "⚠️ ", "FAIL": "❌"}.get(
        report.overall_recommendation, "❓"
    )
    click.echo(f"\n{icon} Kết quả đánh giá cặp mồi: {report.pair_id}")
    click.echo("=" * 60)
    click.echo(f"  Left primer : {report.left_sequence}")
    click.echo(f"  Right primer: {report.right_sequence}")
    if report.probe_sequence:
        click.echo(f"  Probe       : {report.probe_sequence}")
    click.echo("")
    click.echo("📊 Nhiệt động học:")
    lp = report.left_thermo
    rp = report.right_thermo
    click.echo(
        f"  Left  — Tm={lp.tm:.1f}°C | GC={lp.gc_percent:.1f}% | "
        f"Hairpin ΔG={lp.hairpin_dg:.2f} | Self-dimer ΔG={lp.homodimer_dg:.2f} kcal/mol"
    )
    click.echo(
        f"  Right — Tm={rp.tm:.1f}°C | GC={rp.gc_percent:.1f}% | "
        f"Hairpin ΔG={rp.hairpin_dg:.2f} | Self-dimer ΔG={rp.homodimer_dg:.2f} kcal/mol"
    )
    click.echo(
        f"  Hetero-dimer (left×right) ΔG={report.heterodimer_dg:.2f} kcal/mol"
    )
    if report.probe_thermo:
        pp = report.probe_thermo
        click.echo(
            f"  Probe — Tm={pp.tm:.1f}°C | GC={pp.gc_percent:.1f}% | "
            f"Hairpin ΔG={pp.hairpin_dg:.2f} | Self-dimer ΔG={pp.homodimer_dg:.2f} kcal/mol"
        )
    click.echo("")

    if report.blast_performed and report.specificity:
        sp = report.specificity
        click.echo(
            f"🔍 BLAST Specificity: score={sp.specificity_score:.0f}/100 | "
            f"Off-target amplicons={len(sp.off_target_amplicons)}"
        )
    elif not report.blast_performed:
        click.echo("🔍 BLAST: bỏ qua (--no-blast)")
    click.echo("")

    if report.summary_warnings:
        click.echo("⚠️  Cảnh báo:")
        for w in report.summary_warnings:
            click.echo(f"   • {w}")
        click.echo("")
    elif not report.blast_performed:
        click.echo("ℹ️  Ghi chú: BLAST bị bỏ qua — không có dữ liệu specificity.")
        click.echo("")

    click.echo(f"📋 Khuyến nghị: {icon} {report.overall_recommendation}")

    # --- Save JSON if output_dir given ---
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_file = output_dir / f"evaluation_{pair_id}.json"
        try:
            data = asdict(report)
        except Exception:
            data = {"pair_id": pair_id, "recommendation": report.overall_recommendation}
        out_file.write_text(
            json.dumps(data, indent=2, ensure_ascii=False, default=str),
            encoding="utf-8",
        )
        click.echo(f"\n💾 Báo cáo JSON đã lưu: {out_file}")

    sys.exit(0 if report.overall_recommendation != "FAIL" else 1)
