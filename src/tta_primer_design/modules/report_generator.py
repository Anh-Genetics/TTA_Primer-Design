"""Module 11 — ReportGenerator: tạo báo cáo đầu ra đa định dạng.

Output files:
    1. results_summary.csv   — Tổng hợp tất cả primer pairs
    2. results_detailed.xlsx — Excel: Summary, Primers, Probes, BLAST_hits, SNPs
    3. primers.fasta         — Sequences để đặt hàng tổng hợp
    4. results.json          — Full data cho downstream tools
    5. report.html           — Báo cáo trực quan (Jinja2 template)

Visualizations (theo kế hoạch, Sprint 4):
    - Amplicon position map trên transcript/gene
    - Tm distribution chart
    - GC% scatter plot
    - Specificity score bar chart

TODO (Sprint 4):
    - Implement generate_csv()
    - Implement generate_excel() với nhiều sheet
    - Implement generate_json()
    - Implement generate_fasta()
    - Implement generate_html() với Jinja2 template
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from pathlib import Path

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignResult

logger = logging.getLogger("tta_primer_design.modules.report_generator")

_SUPPORTED_FORMATS = {"csv", "xlsx", "json", "fasta", "html"}


class ReportGenerator:
    """Tạo báo cáo kết quả thiết kế primer.

    Args:
        config: AppConfig đã load.

    Example::

        reporter = ReportGenerator(config)
        reporter.generate(results, output_dir="results/run_001/")
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def generate(
        self,
        results: list[DesignResult],
        output_dir: str | Path,
    ) -> list[Path]:
        """Tạo tất cả báo cáo theo config.output.formats.

        Args:
            results: Danh sách DesignResult.
            output_dir: Thư mục lưu output.

        Returns:
            Danh sách đường dẫn file đã tạo.
        """
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        generated: list[Path] = []
        for fmt in self.config.output.formats:
            fmt = fmt.lower()
            if fmt not in _SUPPORTED_FORMATS:
                logger.warning("Unsupported output format: %s — skipping", fmt)
                continue

            try:
                if fmt == "csv":
                    path = self.generate_csv(results, out_path / "results_summary.csv")
                elif fmt == "json":
                    path = self.generate_json(results, out_path / "results.json")
                elif fmt == "xlsx":
                    path = self.generate_excel(results, out_path / "results_detailed.xlsx")
                elif fmt == "fasta":
                    path = self.generate_fasta(results, out_path / "primers.fasta")
                elif fmt == "html":
                    path = self.generate_html(results, out_path / "report.html")
                else:
                    continue

                generated.append(path)
                logger.info("Generated report: %s", path.name)

            except NotImplementedError:
                logger.debug("Report format '%s' not yet implemented", fmt)
            except Exception as exc:  # noqa: BLE001
                logger.error("Failed to generate '%s' report: %s", fmt, exc)

        return generated

    def generate_csv(self, results: list[DesignResult], output_path: Path) -> Path:
        """Tạo CSV tổng hợp.

        Args:
            results: Danh sách DesignResult.
            output_path: Đường dẫn file output.

        Returns:
            Đường dẫn file đã tạo.
        """
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError("pandas is required for CSV generation") from exc
        rows = []
        for r in results:
            for pair in r.primer_pairs:
                rows.append(
                    {
                        "target_id": r.target.target_id,
                        "pair_id": pair.pair_id,
                        "status": r.status,
                        "left_primer": pair.left_primer.sequence,
                        "right_primer": pair.right_primer.sequence,
                        "probe": pair.probe.sequence if pair.probe else "",
                        "amplicon_size": pair.amplicon_size,
                        "left_tm": round(pair.left_primer.tm, 2),
                        "right_tm": round(pair.right_primer.tm, 2),
                        "probe_tm": round(pair.probe.tm, 2) if pair.probe else "",
                        "left_gc_pct": round(pair.left_primer.gc_percent, 2),
                        "right_gc_pct": round(pair.right_primer.gc_percent, 2),
                        "pair_penalty": round(pair.pair_penalty, 4),
                        "score": round(pair.score, 2),
                    }
                )
            if not r.primer_pairs:
                rows.append(
                    {
                        "target_id": r.target.target_id,
                        "pair_id": "",
                        "status": r.status,
                        "left_primer": "",
                        "right_primer": "",
                        "probe": "",
                        "amplicon_size": "",
                        "left_tm": "",
                        "right_tm": "",
                        "probe_tm": "",
                        "left_gc_pct": "",
                        "right_gc_pct": "",
                        "pair_penalty": "",
                        "score": "",
                    }
                )
        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False, encoding="utf-8")
        return output_path

    def generate_excel(self, results: list[DesignResult], output_path: Path) -> Path:
        """Tạo Excel với nhiều sheet.

        Args:
            results: Danh sách DesignResult.
            output_path: Đường dẫn file .xlsx.

        Returns:
            Đường dẫn file đã tạo.
        """
        try:
            from openpyxl import Workbook
        except ImportError as exc:
            raise ImportError("openpyxl is required for Excel generation") from exc

        wb = Workbook()

        # Sheet 1: Summary
        ws_summary = wb.active
        ws_summary.title = "Summary"
        ws_summary.append(["target_id", "status", "num_pairs", "best_score"])
        for r in results:
            best_score = max((p.score for p in r.primer_pairs), default=0.0)
            ws_summary.append([r.target.target_id, r.status, len(r.primer_pairs), best_score])

        # Sheet 2: Primers
        ws_primers = wb.create_sheet("Primers")
        ws_primers.append(
            [
                "target_id",
                "pair_id",
                "left_primer",
                "right_primer",
                "amplicon_size",
                "left_tm",
                "right_tm",
                "left_gc_pct",
                "right_gc_pct",
                "pair_penalty",
                "score",
            ]
        )
        for r in results:
            for pair in r.primer_pairs:
                ws_primers.append(
                    [
                        r.target.target_id,
                        pair.pair_id,
                        pair.left_primer.sequence,
                        pair.right_primer.sequence,
                        pair.amplicon_size,
                        round(pair.left_primer.tm, 2),
                        round(pair.right_primer.tm, 2),
                        round(pair.left_primer.gc_percent, 2),
                        round(pair.right_primer.gc_percent, 2),
                        round(pair.pair_penalty, 4),
                        round(pair.score, 2),
                    ]
                )

        # Sheet 3: Probes
        ws_probes = wb.create_sheet("Probes")
        ws_probes.append(["target_id", "pair_id", "probe_sequence", "probe_tm", "probe_gc_pct"])
        for r in results:
            for pair in r.primer_pairs:
                if pair.probe:
                    ws_probes.append(
                        [
                            r.target.target_id,
                            pair.pair_id,
                            pair.probe.sequence,
                            round(pair.probe.tm, 2),
                            round(pair.probe.gc_percent, 2),
                        ]
                    )

        wb.save(output_path)
        return output_path

    def generate_json(self, results: list[DesignResult], output_path: Path) -> Path:
        """Tạo JSON đầy đủ.

        Args:
            results: Danh sách DesignResult.
            output_path: Đường dẫn file .json.

        Returns:
            Đường dẫn file đã tạo.
        """
        serializable = []
        for r in results:
            try:
                serializable.append(asdict(r))
            except Exception as exc:  # noqa: BLE001
                logger.warning("Could not serialize result for %s: %s", r.target.target_id, exc)
                serializable.append({"target_id": r.target.target_id, "status": r.status})

        with output_path.open("w", encoding="utf-8") as fh:
            json.dump(serializable, fh, indent=2, ensure_ascii=False, default=str)

        return output_path

    def generate_fasta(self, results: list[DesignResult], output_path: Path) -> Path:
        """Tạo FASTA file chứa tất cả oligo sequences.

        Args:
            results: Danh sách DesignResult.
            output_path: Đường dẫn file .fasta.

        Returns:
            Đường dẫn file đã tạo.
        """
        lines = []
        for r in results:
            for pair in r.primer_pairs:
                lines.append(f">{pair.pair_id}_LEFT")
                lines.append(pair.left_primer.sequence)
                lines.append(f">{pair.pair_id}_RIGHT")
                lines.append(pair.right_primer.sequence)
                if pair.probe:
                    lines.append(f">{pair.pair_id}_PROBE")
                    lines.append(pair.probe.sequence)
        with output_path.open("w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))
            if lines:
                fh.write("\n")
        return output_path

    def generate_html(self, results: list[DesignResult], output_path: Path) -> Path:
        """Tạo báo cáo HTML với Jinja2.

        Args:
            results: Danh sách DesignResult.
            output_path: Đường dẫn file .html.

        Returns:
            Đường dẫn file đã tạo.
        """
        try:
            from jinja2 import Template
        except ImportError as exc:
            raise ImportError("jinja2 is required for HTML generation") from exc

        _TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>TTA Primer Design Report</title>
<style>
body{font-family:sans-serif;margin:2em}
table{border-collapse:collapse;width:100%}
th,td{border:1px solid #ccc;padding:6px 10px;text-align:left}
th{background:#4a7ebf;color:#fff}
tr:nth-child(even){background:#f2f2f2}
h1{color:#2c3e50}
.status-success{color:green} .status-failed{color:red} .status-no_primers{color:orange}
</style>
</head>
<body>
<h1>TTA Primer Design Report</h1>
<p>Total targets: {{ results|length }}</p>
{% for r in results %}
<h2>Target: {{ r.target.target_id }} — <span class="status-{{ r.status }}">{{ r.status }}</span></h2>
{% if r.primer_pairs %}
<table>
<tr><th>Pair ID</th><th>Left Primer</th><th>Right Primer</th><th>Probe</th>
    <th>Amplicon (bp)</th><th>Left Tm</th><th>Right Tm</th><th>Score</th></tr>
{% for pair in r.primer_pairs %}
<tr>
  <td>{{ pair.pair_id }}</td>
  <td>{{ pair.left_primer.sequence }}</td>
  <td>{{ pair.right_primer.sequence }}</td>
  <td>{{ pair.probe.sequence if pair.probe else '' }}</td>
  <td>{{ pair.amplicon_size }}</td>
  <td>{{ '%.1f'|format(pair.left_primer.tm) }}</td>
  <td>{{ '%.1f'|format(pair.right_primer.tm) }}</td>
  <td>{{ '%.2f'|format(pair.score) }}</td>
</tr>
{% endfor %}
</table>
{% else %}
<p><em>No primer pairs found.</em></p>
{% endif %}
{% endfor %}
</body></html>"""

        template = Template(_TEMPLATE)
        html = template.render(results=results)
        with output_path.open("w", encoding="utf-8") as fh:
            fh.write(html)
        return output_path
