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

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("generate_csv chưa được implement (Sprint 4)")

    def generate_excel(self, results: list[DesignResult], output_path: Path) -> Path:
        """Tạo Excel với nhiều sheet.

        Args:
            results: Danh sách DesignResult.
            output_path: Đường dẫn file .xlsx.

        Returns:
            Đường dẫn file đã tạo.

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("generate_excel chưa được implement (Sprint 4)")

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

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("generate_fasta chưa được implement (Sprint 4)")

    def generate_html(self, results: list[DesignResult], output_path: Path) -> Path:
        """Tạo báo cáo HTML với Jinja2.

        Args:
            results: Danh sách DesignResult.
            output_path: Đường dẫn file .html.

        Returns:
            Đường dẫn file đã tạo.

        Raises:
            NotImplementedError: Chưa implement (Sprint 4).
        """
        raise NotImplementedError("generate_html chưa được implement (Sprint 4)")
