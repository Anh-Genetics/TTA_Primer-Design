"""Module 5 — NCBIPrimerBlast: tương tác với NCBI Primer-BLAST Web API.

Endpoint: https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi

Luồng:
    Step 1: POST request với sequence + params → nhận Job ID (CMD=request)
    Step 2: Polling GET với Job ID → check status (CMD=Get)
    Step 3: Parse HTML response khi complete

Tham số quan trọng:
    - INPUT_TYPE / SEQUENCE_TEMPLATE
    - ORGANISM, DATABASE, SEARCH_SPECIFIC_PRIMER
    - PRIMER_OPT_SIZE, PRIMER_OPT_TM, PRIMER_PRODUCT_SIZE_RANGE, v.v.
"""

from __future__ import annotations

import logging
import re
import time
from typing import Any

import requests
from bs4 import BeautifulSoup

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget, Oligo, PrimerPair, ProcessedSequence

logger = logging.getLogger("tta_primer_design.modules.ncbi_primer_blast")

_BASE_URL = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"
_POLL_INTERVAL = 10  # giây
_MAX_WAIT = 300  # giây tối đa chờ kết quả

# Keyword patterns in the NCBI Primer-BLAST HTML response
_STATUS_RUNNING = [
    "Your job is still running",
    "Job still running",
    "still running",
    "Please wait",
]
_STATUS_COMPLETE = [
    "Primer pairs found",
    "primers found",
    "Graphic view",
]
_STATUS_NOPRIMERS = [
    "No primers found",
    "no primers were found",
]
_STATUS_ERROR = [
    "Error",
    "error occurred",
    "Internal Server Error",
]

# Default primer design parameters for NCBI Primer-BLAST
_DEFAULT_PARAMS: dict[str, Any] = {
    "INPUT_TYPE": "0",  # 0 = template sequence
    "PRIMER_SPECIFICITY_DATABASE": "7",  # refseq_rna
    "SEARCH_SPECIFIC_PRIMER": "1",
    "PRIMER_OPT_SIZE": "20",
    "PRIMER_MIN_SIZE": "18",
    "PRIMER_MAX_SIZE": "25",
    "PRIMER_OPT_TM": "60.0",
    "PRIMER_MIN_TM": "57.0",
    "PRIMER_MAX_TM": "63.0",
    "PRIMER_OPT_GC_PERCENT": "50",
    "PRIMER_MIN_GC": "40",
    "PRIMER_MAX_GC": "60",
    "PRIMER_PRODUCT_MIN": "70",
    "PRIMER_PRODUCT_MAX": "200",
    "NUM_TARGETS_WITH_PRIMERS": "20",
    "NUM_TARGETS_WITH_PRIMERS_USER": "10",
    "PRIMER_NUM_RETURN": "10",
    "PRIMER_MAX_DIFF_TM": "3",
}


class NCBIPrimerBlast:
    """Tương tác với NCBI Primer-BLAST Web API.

    Args:
        config: AppConfig đã load.

    Example::

        api = NCBIPrimerBlast(config)
        pairs = api.design_primers(target, processed_seq)
    """

    BASE_URL = _BASE_URL
    POLL_INTERVAL = _POLL_INTERVAL
    MAX_WAIT = _MAX_WAIT

    def __init__(self, config: AppConfig) -> None:
        self.config = config
        self._session = requests.Session()

    def submit_job(self, params: dict[str, Any]) -> str:
        """Submit job lên NCBI Primer-BLAST.

        POST request tới NCBI Primer-BLAST với các tham số thiết kế primer.
        Trích xuất JOB_KEY từ HTML response.

        Args:
            params: Dict các tham số gửi lên API.

        Returns:
            Job key (string) dùng để poll kết quả.

        Raises:
            ValueError: Nếu không tìm thấy JOB_KEY trong response.
            requests.RequestException: Nếu gặp lỗi mạng.
        """
        if self.config.ncbi.email:
            params.setdefault("EMAIL", self.config.ncbi.email)

        resp = self._session.post(
            self.BASE_URL,
            data=params,
            timeout=self.config.ncbi.timeout,
        )
        resp.raise_for_status()

        job_key = self._extract_job_key(resp.text)
        if not job_key:
            raise ValueError("Cannot extract JOB_KEY from NCBI Primer-BLAST response")

        logger.info("Submitted NCBI Primer-BLAST job: %s", job_key)
        return job_key

    def poll_status(self, job_key: str) -> str:
        """Kiểm tra trạng thái job.

        GET request với JOB_KEY; phân tích text HTML để xác định status.

        Args:
            job_key: Job key nhận được từ submit_job.

        Returns:
            "running" | "complete" | "no_primers" | "error"
        """
        try:
            resp = self._session.get(
                self.BASE_URL,
                params={"CMD": "Get", "JOB_KEY": job_key},
                timeout=self.config.ncbi.timeout,
            )
            resp.raise_for_status()
        except requests.RequestException as exc:
            logger.warning("Poll error for job %s: %s", job_key, exc)
            return "error"

        html = resp.text
        return self._classify_status(html)

    def fetch_results(self, job_key: str) -> str:
        """Lấy kết quả HTML khi job hoàn thành.

        Args:
            job_key: Job key.

        Returns:
            HTML string của trang kết quả.

        Raises:
            requests.RequestException: Nếu gặp lỗi mạng.
        """
        resp = self._session.get(
            self.BASE_URL,
            params={"CMD": "Get", "JOB_KEY": job_key},
            timeout=self.config.ncbi.timeout,
        )
        resp.raise_for_status()
        return resp.text

    def parse_html_results(self, html: str) -> list[PrimerPair]:
        """Parse HTML kết quả thành List[PrimerPair].

        Sử dụng BeautifulSoup để trích xuất thông tin primer pairs từ
        HTML của NCBI Primer-BLAST kết quả.

        Args:
            html: HTML string từ NCBI Primer-BLAST.

        Returns:
            Danh sách PrimerPair.
        """
        soup = BeautifulSoup(html, "lxml")
        primer_pairs: list[PrimerPair] = []

        # NCBI Primer-BLAST result tables have class "prPrimer" rows
        # Each primer pair section starts with "Primer pair <N>"
        pair_sections = self._find_primer_pair_sections(soup)

        for idx, section in enumerate(pair_sections):
            pair = self._parse_primer_pair_section(section, idx)
            if pair is not None:
                primer_pairs.append(pair)

        if not primer_pairs:
            # Fallback: try regex-based extraction
            primer_pairs = self._parse_html_fallback(html)

        logger.info("Parsed %d primer pairs from NCBI Primer-BLAST HTML", len(primer_pairs))
        return primer_pairs

    def design_primers(
        self,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> list[PrimerPair]:
        """Thiết kế primer thông qua NCBI Primer-BLAST API.

        Wrapper tổng hợp: submit → poll với exponential backoff → parse.

        Args:
            target: DesignTarget.
            processed_seq: ProcessedSequence.

        Returns:
            Danh sách PrimerPair từ NCBI.

        Raises:
            TimeoutError: Nếu job không hoàn thành trong MAX_WAIT giây.
            RuntimeError: Nếu NCBI trả về lỗi.
        """
        params = self._build_params(target, processed_seq)
        job_key = self.submit_job(params)

        # Poll with exponential backoff
        elapsed = 0
        interval = self.POLL_INTERVAL

        while elapsed < self.MAX_WAIT:
            time.sleep(interval)
            elapsed += interval

            status = self.poll_status(job_key)
            logger.debug("Job %s status: %s (%ds elapsed)", job_key, status, elapsed)

            if status == "complete":
                html = self.fetch_results(job_key)
                return self.parse_html_results(html)

            if status == "no_primers":
                logger.info("NCBI Primer-BLAST found no primers for job %s", job_key)
                return []

            if status == "error":
                raise RuntimeError(f"NCBI Primer-BLAST job {job_key} reported an error")

            # Exponential backoff, cap at 60 seconds
            interval = min(interval * 1.5, 60)

        raise TimeoutError(
            f"NCBI Primer-BLAST job {job_key} did not complete within {self.MAX_WAIT}s"
        )

    # ------------------------------------------------------------------
    # Private helpers — HTTP / parsing
    # ------------------------------------------------------------------

    def _build_params(
        self, target: DesignTarget, processed_seq: ProcessedSequence
    ) -> dict[str, Any]:
        """Xây dựng dict tham số từ target và processed_seq.

        Args:
            target: DesignTarget.
            processed_seq: ProcessedSequence.

        Returns:
            Dict tham số cho NCBI Primer-BLAST POST request.
        """
        params: dict[str, Any] = dict(_DEFAULT_PARAMS)
        params["SEQUENCE_TEMPLATE"] = processed_seq.sequence

        # Map organism name to NCBI taxid where common
        organism_map = {
            "Homo sapiens": "9606",
            "Mus musculus": "10090",
            "Rattus norvegicus": "10116",
        }
        organism_id = organism_map.get(target.organism, target.organism)
        params["ORGANISM"] = organism_id

        # Apply custom params from target
        if target.custom_params:
            params.update(target.custom_params)

        # Product size range
        if target.region_include:
            start, end = target.region_include
            params["SEQUENCE_TARGET"] = f"{start},{end - start}"

        return params

    @staticmethod
    def _extract_job_key(html: str) -> str | None:
        """Trích xuất JOB_KEY từ HTML response.

        Args:
            html: HTML string từ NCBI Primer-BLAST POST response.

        Returns:
            JOB_KEY string hoặc None nếu không tìm thấy.
        """
        # Try hidden input field first
        match = re.search(r'<input[^>]+name="JOB_KEY"[^>]+value="([^"]+)"', html, re.IGNORECASE)
        if match:
            return match.group(1)

        # Try value before name ordering
        match = re.search(r'<input[^>]+value="([^"]+)"[^>]+name="JOB_KEY"', html, re.IGNORECASE)
        if match:
            return match.group(1)

        # Try URL-style parameter
        match = re.search(r"JOB_KEY=([A-Za-z0-9_\-]+)", html)
        if match:
            return match.group(1)

        return None

    @staticmethod
    def _classify_status(html: str) -> str:
        """Phân loại trạng thái job từ HTML.

        Args:
            html: HTML của poll response.

        Returns:
            "running" | "complete" | "no_primers" | "error"
        """
        for pattern in _STATUS_NOPRIMERS:
            if pattern.lower() in html.lower():
                return "no_primers"
        for pattern in _STATUS_COMPLETE:
            if pattern.lower() in html.lower():
                return "complete"
        for pattern in _STATUS_ERROR:
            if pattern.lower() in html.lower():
                return "error"
        for pattern in _STATUS_RUNNING:
            if pattern.lower() in html.lower():
                return "running"
        # Default: still running if can't determine
        return "running"

    @staticmethod
    def _find_primer_pair_sections(soup: BeautifulSoup) -> list[Any]:
        """Tìm các section chứa primer pair trong HTML.

        Args:
            soup: BeautifulSoup object.

        Returns:
            Danh sách các table/section elements chứa primer pair.
        """
        sections: list[Any] = []

        # Strategy 1: Find tables with primer-related headers
        for table in soup.find_all("table"):
            text = table.get_text()
            if "Forward primer" in text or "forward primer" in text.lower():
                sections.append(table)

        # Strategy 2: Find divs/sections with "Primer pair" heading
        if not sections:
            for elem in soup.find_all(["div", "section", "p"]):
                if re.search(r"Primer pair\s+\d+", elem.get_text()):
                    sections.append(elem)

        return sections

    def _parse_primer_pair_section(self, section: Any, idx: int) -> PrimerPair | None:
        """Parse một section HTML thành PrimerPair.

        Xử lý cấu trúc bảng của NCBI Primer-BLAST kết quả:
        mỗi bảng có hàng Forward/Reverse primer với Sequence, Tm, GC%.

        Args:
            section: BeautifulSoup element chứa primer pair info.
            idx: Index (0-based) của primer pair.

        Returns:
            PrimerPair hoặc None nếu parse thất bại.
        """
        rows = section.find_all("tr")
        left_data: dict[str, str] = {}
        right_data: dict[str, str] = {}
        amplicon_size = 0

        for row in rows:
            cells = [td.get_text(strip=True) for td in row.find_all(["td", "th"])]
            if not cells:
                continue

            row_text = " ".join(cells).lower()

            if "forward" in row_text and len(cells) >= 2:
                left_data = _extract_oligo_data(cells)
            elif "reverse" in row_text and len(cells) >= 2:
                right_data = _extract_oligo_data(cells)
            elif "product length" in row_text or "amplicon" in row_text:
                for cell in cells:
                    if cell.isdigit():
                        amplicon_size = int(cell)
                        break

        if not left_data.get("sequence") or not right_data.get("sequence"):
            return None

        left_primer = Oligo(
            sequence=left_data["sequence"],
            tm=float(left_data.get("tm", 0.0)),
            gc_percent=float(left_data.get("gc", 0.0)),
            start=int(left_data.get("start", 0)),
        )
        right_primer = Oligo(
            sequence=right_data["sequence"],
            tm=float(right_data.get("tm", 0.0)),
            gc_percent=float(right_data.get("gc", 0.0)),
            start=int(right_data.get("start", 0)),
        )

        return PrimerPair(
            pair_id=f"ncbi_{idx}",
            left_primer=left_primer,
            right_primer=right_primer,
            amplicon_size=amplicon_size,
        )

    @staticmethod
    def _parse_html_fallback(html: str) -> list[PrimerPair]:
        """Fallback: trích xuất primer sequences từ HTML bằng regex.

        Args:
            html: HTML string.

        Returns:
            Danh sách PrimerPair (tối thiểu — chỉ có sequences).
        """
        pairs: list[PrimerPair] = []

        # Look for patterns like "Forward primer: ATCGATCG..."
        fwd_pattern = re.compile(r"[Ff]orward\s+primer[^:]*:\s*([ACGTacgt]{10,40})", re.IGNORECASE)
        rev_pattern = re.compile(r"[Rr]everse\s+primer[^:]*:\s*([ACGTacgt]{10,40})", re.IGNORECASE)

        fwd_matches = fwd_pattern.findall(html)
        rev_matches = rev_pattern.findall(html)

        for idx, (fwd_seq, rev_seq) in enumerate(zip(fwd_matches, rev_matches, strict=False)):
            left = Oligo(sequence=fwd_seq.upper())
            right = Oligo(sequence=rev_seq.upper())
            pairs.append(PrimerPair(pair_id=f"ncbi_{idx}", left_primer=left, right_primer=right))

        return pairs


# ------------------------------------------------------------------
# Module-level helpers
# ------------------------------------------------------------------


def _extract_oligo_data(cells: list[str]) -> dict[str, str]:
    """Trích xuất dữ liệu oligo từ danh sách cells của một hàng bảng.

    Expected cell order (NCBI Primer-BLAST):
    Label | Sequence | Template strand | Length | Start | Stop | Tm | GC% | ...

    Args:
        cells: Danh sách text từ các cells trong một hàng.

    Returns:
        Dict với keys: sequence, tm, gc, start, stop.
    """
    data: dict[str, str] = {}

    for cell in cells:
        # DNA sequence: only ACGTacgt and length 10-40
        if re.match(r"^[ACGTacgt]{10,40}$", cell):
            data["sequence"] = cell.upper()
        # Temperature: number with optional decimal
        elif re.match(r"^\d{2}(\.\d+)?$", cell) and "tm" not in data:
            val = float(cell)
            if 40.0 <= val <= 90.0:
                data["tm"] = cell
        # GC%: number with % sign or 20-80 range
        elif cell.endswith("%") and re.match(r"^\d+", cell):
            data["gc"] = cell.rstrip("%")
        elif re.match(r"^\d+$", cell) and "start" not in data and int(cell) < 10000:
            data["start"] = cell

    return data
