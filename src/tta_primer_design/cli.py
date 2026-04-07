"""CLI cho TTA Primer Design — dùng Click."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import click

from tta_primer_design import __version__
from tta_primer_design.main import run_pipeline

if TYPE_CHECKING:
    from tta_primer_design.models import PrimerPair
    from tta_primer_design.modules.blast_specificity import SpecificityResult

_VALID_BASES = re.compile(r"^[ATCGNatcgn]+$")


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


# ---------------------------------------------------------------------------
# evaluate command
# ---------------------------------------------------------------------------


@main.command("evaluate")
@click.option(
    "--left",
    "-l",
    "left_seq",
    required=True,
    help="Chuỗi primer trái (5'→3').",
)
@click.option(
    "--right",
    "-r",
    "right_seq",
    required=True,
    help="Chuỗi primer phải (5'→3').",
)
@click.option(
    "--probe",
    "-p",
    "probe_seq",
    default=None,
    help="Chuỗi TaqMan probe (tuỳ chọn).",
)
@click.option(
    "--name",
    "-n",
    "pair_name",
    default="pair_001",
    show_default=True,
    help="Tên/ID cho cặp mồi.",
)
@click.option(
    "--organism",
    default="Homo sapiens",
    show_default=True,
    help="Tên loài dùng để lọc BLAST (ví dụ: 'Mus musculus').",
)
@click.option(
    "--database",
    default="nt",
    show_default=True,
    help="NCBI BLAST database: nt | refseq_rna | refseq_genomic.",
)
@click.option(
    "--skip-blast",
    is_flag=True,
    default=False,
    help="Bỏ qua bước BLAST (chỉ đánh giá nhiệt động học).",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    default=None,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Lưu kết quả ra file JSON.",
)
@click.option(
    "--config",
    "-c",
    "config_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="File YAML config (tuỳ chọn).",
)
@click.option(
    "--log-level",
    default="WARNING",
    show_default=True,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    help="Mức độ log.",
)
def evaluate_command(
    left_seq: str,
    right_seq: str,
    probe_seq: str | None,
    pair_name: str,
    organism: str,
    database: str,
    skip_blast: bool,
    output_file: Path | None,
    config_path: Path | None,
    log_level: str,
) -> None:
    """Đánh giá nhiệt động học và kiểm tra BLAST cho một cặp mồi sẵn có.

    Lệnh này nhận trực tiếp trình tự primer/probe, tính toán các thông số
    nhiệt động học (Tm, ΔG, GC%, GC clamp, repeat runs) và tuỳ chọn kiểm tra
    tính đặc hiệu bằng NCBI BLAST (cần kết nối internet).

    \b
    Ví dụ:

        # Chỉ đánh giá nhiệt động học (không cần mạng, nhanh)
        tta-primer-design evaluate \\
            --left  GCAAGGAATGGTTTCAGAAATCCA \\
            --right CAGGACTCCATGTCGTCCA \\
            --skip-blast

        # Đánh giá đầy đủ kèm BLAST (cần internet, ~1-2 phút)
        tta-primer-design evaluate \\
            --left  GCAAGGAATGGTTTCAGAAATCCA \\
            --right CAGGACTCCATGTCGTCCA \\
            --organism "Homo sapiens" --database nt

        # TaqMan — có probe, lưu kết quả ra file JSON
        tta-primer-design evaluate \\
            --left  GCAAGGAATGGTTTCAGAAATCCA \\
            --right CAGGACTCCATGTCGTCCA \\
            --probe TGCAGCCACACTTTCTACAATGAGC \\
            --name  ACTB_pair1 \\
            --output results/eval_ACTB.json
    """
    from tta_primer_design.config import load_config
    from tta_primer_design.logging_setup import setup_logging
    from tta_primer_design.models import Oligo, PrimerPair
    from tta_primer_design.modules.filter_ranker import FilterRanker
    from tta_primer_design.modules.thermodynamics import Thermodynamics

    setup_logging(level=log_level)

    # ---------- Validate sequences ----------
    for label, seq in [("Left", left_seq), ("Right", right_seq)]:
        seq = seq.strip()
        if not seq:
            click.echo(f"❌ {label} primer không được để trống.", err=True)
            sys.exit(2)
        if not _VALID_BASES.match(seq):
            bad = sorted(set(seq.upper()) - set("ATCGN"))
            click.echo(
                f"❌ {label} primer chứa ký tự không hợp lệ: {bad}",
                err=True,
            )
            sys.exit(2)
    if probe_seq is not None:
        probe_seq = probe_seq.strip()
        if probe_seq and not _VALID_BASES.match(probe_seq):
            bad = sorted(set(probe_seq.upper()) - set("ATCGN"))
            click.echo(f"❌ Probe chứa ký tự không hợp lệ: {bad}", err=True)
            sys.exit(2)

    left_seq = left_seq.strip().upper()
    right_seq = right_seq.strip().upper()
    if probe_seq:
        probe_seq = probe_seq.strip().upper()

    # ---------- Load config ----------
    try:
        cfg = load_config(config_path)
    except (FileNotFoundError, ValueError) as exc:
        click.echo(f"❌ Lỗi config: {exc}", err=True)
        sys.exit(2)

    cfg.blast.database = database

    # ---------- Build model objects ----------
    left_oligo = Oligo(sequence=left_seq)
    right_oligo = Oligo(sequence=right_seq)
    probe_oligo = Oligo(sequence=probe_seq) if probe_seq else None

    pair = PrimerPair(
        pair_id=pair_name,
        left_primer=left_oligo,
        right_primer=right_oligo,
        probe=probe_oligo,
    )

    _W = 60  # display width

    click.echo()
    click.echo("=" * _W)
    click.echo("  TTA Primer Design — Evaluate Existing Primer Pair")
    click.echo("=" * _W)
    click.echo(f"  Pair: {pair_name}")
    click.echo()

    # ==================== THERMODYNAMICS ====================
    click.echo("─── THERMODYNAMICS " + "─" * (_W - 19))
    try:
        thermo = Thermodynamics(cfg)

        def _print_oligo_profile(label: str, seq: str, oligo: Oligo) -> None:
            profile = thermo.full_thermodynamic_profile(oligo)
            oligo.tm = profile.tm
            oligo.gc_percent = profile.gc_percent
            oligo.hairpin_th = profile.hairpin_dg
            oligo.self_any_th = profile.self_dimer_dg
            oligo.end_stability = profile.end_stability_dg

            gc_icon = "✅" if profile.gc_clamp_ok else "⚠️ "
            rep_icon = "✅" if profile.repeat_ok else "⚠️ "
            hp_icon = "✅" if profile.hairpin_dg > -9.0 else "⚠️ "
            sd_icon = "✅" if profile.self_dimer_dg > -9.0 else "⚠️ "

            click.echo(f"\n  {label}: {seq}")
            click.echo(f"    Tm             : {profile.tm:.1f} °C")
            click.echo(f"    GC%            : {profile.gc_percent:.1f}%")
            click.echo(f"    Hairpin ΔG     : {profile.hairpin_dg:.2f} kcal/mol  {hp_icon}")
            click.echo(f"    Self-dimer ΔG  : {profile.self_dimer_dg:.2f} kcal/mol  {sd_icon}")
            click.echo(f"    3'-end ΔG      : {profile.end_stability_dg:.2f} kcal/mol")
            click.echo(f"    GC Clamp       : {gc_icon}")
            click.echo(f"    Repeat runs    : {rep_icon}")

        _print_oligo_profile("Left  primer", left_seq, left_oligo)
        _print_oligo_profile("Right primer", right_seq, right_oligo)

        heterodimer_dg = thermo.calculate_dimer_dg(left_seq, right_seq)
        hd_icon = "✅" if heterodimer_dg > -9.0 else "⚠️ "
        click.echo(f"\n  Hetero-dimer ΔG  : {heterodimer_dg:.2f} kcal/mol  {hd_icon}")

        if probe_oligo and probe_seq:
            _print_oligo_profile("Probe       ", probe_seq, probe_oligo)

    except ImportError:
        click.echo("  ⚠️  primer3-py chưa được cài — bỏ qua nhiệt động học.")
        click.echo("      Chạy: pip install primer3-py")
    except Exception as exc:  # noqa: BLE001
        click.echo(f"  ⚠️  Lỗi thermodynamics: {exc}", err=True)

    click.echo()

    # ==================== BLAST SPECIFICITY ====================
    spec_result = None
    if not skip_blast:
        click.echo("─── BLAST SPECIFICITY " + "─" * (_W - 22))
        click.echo(f"  Database : {database} | Organism : {organism}")
        click.echo("  ⏳ Đang chạy BLAST (30–120 giây)…")
        try:
            from tta_primer_design.modules.blast_specificity import BlastSpecificity

            blast_checker = BlastSpecificity(cfg)
            spec_result = blast_checker.check_pair(pair, organism=organism)
            pair.specificity_result = spec_result

            left_perfect = sum(1 for h in spec_result.blast_hits_left if h.mismatches == 0)
            right_perfect = sum(1 for h in spec_result.blast_hits_right if h.mismatches == 0)
            left_3end_ok = sum(1 for h in spec_result.blast_hits_left if h.mismatches_3prime > 0)
            right_3end_ok = sum(1 for h in spec_result.blast_hits_right if h.mismatches_3prime > 0)

            click.echo(
                f"\n  Left primer  : {len(spec_result.blast_hits_left):3d} hits"
                f"  | perfect: {left_perfect}"
                f"  | 3'-mismatch: {left_3end_ok}"
            )
            click.echo(
                f"  Right primer : {len(spec_result.blast_hits_right):3d} hits"
                f"  | perfect: {right_perfect}"
                f"  | 3'-mismatch: {right_3end_ok}"
            )

            if spec_result.off_target_amplicons:
                click.echo(f"\n  ⚠️  Off-target amplicons: {len(spec_result.off_target_amplicons)}")
                for amp in spec_result.off_target_amplicons[:5]:
                    click.echo(f"     • {amp.subject_id}: ~{amp.amplicon_size} bp")
                if len(spec_result.off_target_amplicons) > 5:
                    click.echo(
                        f"     … và {len(spec_result.off_target_amplicons) - 5} amplicon khác"
                    )
            else:
                click.echo("\n  Off-target amplicons : 0  ✅")

            click.echo(f"  Specificity score    : {spec_result.specificity_score:.1f} / 100")

        except ImportError:
            click.echo("  ⚠️  biopython chưa được cài — bỏ qua BLAST.")
            click.echo("      Chạy: pip install biopython")
        except RuntimeError as exc:
            click.echo(f"  ⚠️  BLAST API lỗi: {exc}", err=True)
        except Exception as exc:  # noqa: BLE001
            click.echo(f"  ⚠️  Lỗi BLAST không xác định: {exc}", err=True)

        click.echo()
    else:
        click.echo("  (BLAST bị bỏ qua — bỏ cờ --skip-blast để kiểm tra tính đặc hiệu)")
        click.echo()

    # ==================== SUMMARY ====================
    click.echo("─── SUMMARY " + "─" * (_W - 12))
    try:
        ranker = FilterRanker(cfg)
        score = ranker.calculate_score(pair)
        pair.score = score

        if score >= 70:
            verdict = "✅  PASS"
        elif score >= 50:
            verdict = "⚠️   MARGINAL"
        else:
            verdict = "❌  FAIL"

        click.echo(f"\n  Overall Score : {score:.1f} / 100  →  {verdict}")
        if pair.snp_flags:
            click.echo(f"  Flags         : {', '.join(pair.snp_flags)}")

    except Exception as exc:  # noqa: BLE001
        click.echo(f"  ⚠️  Lỗi tính điểm: {exc}", err=True)

    click.echo()

    # ==================== SAVE JSON ====================
    if output_file is not None:
        _save_evaluate_json(
            output_file,
            pair,
            left_seq,
            right_seq,
            probe_seq,
            spec_result,
        )


def _save_evaluate_json(
    output_file: Path,
    pair: PrimerPair,
    left_seq: str,
    right_seq: str,
    probe_seq: str | None,
    spec_result: SpecificityResult | None,
) -> None:
    """Lưu kết quả evaluate ra file JSON."""
    data: dict = {
        "pair_id": pair.pair_id,
        "left_primer": {
            "sequence": left_seq,
            "tm": round(pair.left_primer.tm, 2),
            "gc_percent": round(pair.left_primer.gc_percent, 2),
            "hairpin_th": round(pair.left_primer.hairpin_th, 2),
            "self_any_th": round(pair.left_primer.self_any_th, 2),
            "end_stability": round(pair.left_primer.end_stability, 2),
        },
        "right_primer": {
            "sequence": right_seq,
            "tm": round(pair.right_primer.tm, 2),
            "gc_percent": round(pair.right_primer.gc_percent, 2),
            "hairpin_th": round(pair.right_primer.hairpin_th, 2),
            "self_any_th": round(pair.right_primer.self_any_th, 2),
            "end_stability": round(pair.right_primer.end_stability, 2),
        },
        "score": round(pair.score, 2),
    }
    if probe_seq and pair.probe is not None:
        data["probe"] = {
            "sequence": probe_seq,
            "tm": round(pair.probe.tm, 2),
            "gc_percent": round(pair.probe.gc_percent, 2),
        }
    if spec_result is not None:
        data["blast"] = {
            "left_hits": len(spec_result.blast_hits_left),
            "right_hits": len(spec_result.blast_hits_right),
            "off_target_amplicons": len(spec_result.off_target_amplicons),
            "specificity_score": round(spec_result.specificity_score, 2),
            "is_specific": spec_result.is_specific,
        }

    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(json.dumps(data, indent=2, ensure_ascii=False))
        click.echo(f"💾 Kết quả đã lưu: {output_file}")
    except OSError as exc:
        click.echo(f"⚠️  Không thể lưu file: {exc}", err=True)
