"""Module 1 — InputParser: parse các định dạng input khác nhau.

Định dạng được hỗ trợ:
    - JSON  : ``{"targets": [...]}`` — xem ``tests/fixtures/sample_input.json``
    - CSV   : cột ``id``, ``accession`` (hoặc ``sequence``), các cột tuỳ chọn
    - FASTA : mỗi record là một target (ID = header)
    - TXT   : mỗi dòng là một accession NCBI

Output:
    List[DesignTarget]
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget

logger = logging.getLogger("tta_primer_design.modules.input_parser")

_SUPPORTED_EXTENSIONS = {".json", ".csv", ".fasta", ".fa", ".txt"}


class InputParser:
    """Parse file input thành danh sách DesignTarget.

    Args:
        config: AppConfig đã load.

    Example::

        parser = InputParser(config)
        targets = parser.parse("data/input/targets.csv")
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def parse(self, input_path: str | Path) -> list[DesignTarget]:
        """Parse file input.

        Args:
            input_path: Đường dẫn tới file input.

        Returns:
            Danh sách DesignTarget.

        Raises:
            FileNotFoundError: Nếu file không tồn tại.
            ValueError: Nếu định dạng file không hỗ trợ hoặc bị lỗi.
        """
        path = Path(input_path)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {path}")

        ext = path.suffix.lower()
        logger.info("Parsing input file: %s (format: %s)", path.name, ext)

        if ext == ".json":
            targets = self._parse_json(path)
        elif ext == ".csv":
            targets = self._parse_csv(path)
        elif ext in {".fasta", ".fa"}:
            targets = self._parse_fasta(path)
        elif ext == ".txt":
            targets = self._parse_accession_list(path)
        else:
            raise ValueError(
                f"Unsupported input format '{ext}'. " f"Supported: {sorted(_SUPPORTED_EXTENSIONS)}"
            )

        logger.info("Parsed %d target(s)", len(targets))
        return targets

    # ------------------------------------------------------------------
    # Private parsers
    # ------------------------------------------------------------------

    def _parse_json(self, path: Path) -> list[DesignTarget]:
        """Parse JSON file.

        Expected schema::

            {
              "targets": [
                {
                  "id": "ACTB",
                  "accession": "NM_001101",
                  "organism": "Homo sapiens",
                  "design_mode": "qpcr"
                }
              ]
            }
        """
        with path.open("r", encoding="utf-8") as fh:
            data: dict[str, Any] = json.load(fh)

        raw_targets = data.get("targets", [])
        if not isinstance(raw_targets, list):
            raise ValueError("JSON file must have a top-level 'targets' list")

        targets: list[DesignTarget] = []
        for i, item in enumerate(raw_targets):
            if not isinstance(item, dict):
                raise ValueError(f"Target #{i} is not a JSON object")
            targets.append(self._dict_to_target(item))
        return targets

    def _parse_csv(self, path: Path) -> list[DesignTarget]:
        """Parse CSV file.

        Required columns: ``id``
        Optional columns: ``accession``, ``sequence``, ``gene_name``,
                          ``organism``, ``design_mode``, ``exon_junction``
        """
        targets: list[DesignTarget] = []
        with path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for i, row in enumerate(reader):
                if "id" not in row:
                    raise ValueError(f"CSV row #{i + 1} is missing required column 'id'")
                targets.append(self._dict_to_target(dict(row)))
        return targets

    def _parse_fasta(self, path: Path) -> list[DesignTarget]:
        """Parse FASTA file — mỗi record là một DesignTarget.

        Dùng BioPython nếu có; fallback sang parser đơn giản.
        """
        try:
            from Bio import SeqIO  # type: ignore[import-untyped]

            targets: list[DesignTarget] = []
            for record in SeqIO.parse(str(path), "fasta"):
                targets.append(
                    DesignTarget(
                        target_id=record.id,
                        input_type="sequence",
                        sequence=str(record.seq).upper(),
                        organism=self.config.blast.organism,
                        design_mode=self.config.pipeline.mode,
                    )
                )
            return targets
        except ImportError:
            return self._parse_fasta_simple(path)

    def _parse_fasta_simple(self, path: Path) -> list[DesignTarget]:
        """FASTA parser đơn giản (không cần BioPython)."""
        targets: list[DesignTarget] = []
        current_id: str | None = None
        seq_lines: list[str] = []

        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if current_id is not None:
                        targets.append(
                            DesignTarget(
                                target_id=current_id,
                                input_type="sequence",
                                sequence="".join(seq_lines).upper(),
                                organism=self.config.blast.organism,
                                design_mode=self.config.pipeline.mode,
                            )
                        )
                    current_id = line[1:].split()[0]
                    seq_lines = []
                else:
                    seq_lines.append(line)

        if current_id is not None:
            targets.append(
                DesignTarget(
                    target_id=current_id,
                    input_type="sequence",
                    sequence="".join(seq_lines).upper(),
                    organism=self.config.blast.organism,
                    design_mode=self.config.pipeline.mode,
                )
            )

        return targets

    def _parse_accession_list(self, path: Path) -> list[DesignTarget]:
        """Parse danh sách accession (mỗi dòng một accession)."""
        targets: list[DesignTarget] = []
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                acc = line.strip()
                if acc and not acc.startswith("#"):
                    targets.append(
                        DesignTarget(
                            target_id=acc,
                            input_type="accession",
                            accession=acc,
                            organism=self.config.blast.organism,
                            design_mode=self.config.pipeline.mode,
                        )
                    )
        return targets

    # ------------------------------------------------------------------
    # Helper
    # ------------------------------------------------------------------

    def _dict_to_target(self, d: dict[str, Any]) -> DesignTarget:
        """Chuyển dict (từ JSON hoặc CSV) thành DesignTarget."""
        target_id = str(d.get("id", d.get("target_id", "unknown")))
        accession = d.get("accession") or None
        sequence = d.get("sequence") or None
        gene_name = d.get("gene_name") or None

        # Xác định input_type tự động nếu chưa có
        input_type = d.get("input_type", None)
        if input_type is None:
            if accession:
                input_type = "accession"
            elif sequence:
                input_type = "sequence"
            elif gene_name:
                input_type = "gene_name"
            else:
                input_type = "accession"

        exon_junction_raw = d.get("exon_junction", False)
        exon_junction: bool
        if isinstance(exon_junction_raw, str):
            exon_junction = exon_junction_raw.strip().lower() in {"true", "1", "yes"}
        else:
            exon_junction = bool(exon_junction_raw)

        return DesignTarget(
            target_id=target_id,
            input_type=input_type,
            accession=accession,
            sequence=sequence.upper() if sequence else None,
            gene_name=gene_name,
            organism=d.get("organism", self.config.blast.organism),
            design_mode=d.get("design_mode", self.config.pipeline.mode),
            exon_junction=exon_junction,
        )
