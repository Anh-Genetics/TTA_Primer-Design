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

# Cụm từ xuất hiện trong HTML khi job vẫn đang chạy
_RUNNING_PATTERNS = [
    "job is still running",
    "being processed",
    "please wait",
    "job id#",
]
# Cụm từ xuất hiện khi có kết quả
_COMPLETE_PATTERNS = [
    "primer pair",
    "left primer",
    "right primer",
    "no primers were found",
    "no primer pairs",
]
# Cụm từ xuất hiện khi có lỗi
_ERROR_PATTERNS = [
    "error occurred",
    "server error",
    "invalid sequence",
    "too short",
    "too long",
]

# Tham số mặc định gửi lên Primer-BLAST
_DEFAULT_PARAMS: dict[str, str] = {
    "CMD": "request",
    "INPUT_TYPE": "0",  # 0 = sequence template
    "PRIMER_SPECIFICITY_DATABASE": "refseq_mrna",
    "SEARCH_SPECIFIC_PRIMER": "on",
    "PRIMER_MIN_SIZE": "18",
    "PRIMER_OPT_SIZE": "20",
    "PRIMER_MAX_SIZE": "25",
    "PRIMER_MIN_TM": "57.0",
    "PRIMER_OPT_TM": "60.0",
    "PRIMER_MAX_TM": "63.0",
    "PRIMER_MIN_GC": "40.0",
    "PRIMER_MAX_GC": "60.0",
    "PRIMER_MAX_SELF_ANY_TH": "45.0",
    "PRIMER_MAX_SELF_END_TH": "35.0",
    "PRIMER_MAX_HAIRPIN_TH": "24.0",
    "PRIMER_NUM_RETURN": "10",
    "PRIMER_PRODUCT_SIZE_RANGE": "70-200",
    "SHOW_SVIEWER": "true",
}


def _extract_job_key(html: str) -> str:
    """Trích xuất job key từ HTML response của NCBI Primer-BLAST.

    Tìm theo thứ tự ưu tiên:
        1. Hidden input field ``name="JOB_KEY"``
        2. Pattern ``JOB_KEY=<value>`` trong URL/text
        3. Pattern ``job id# <value>`` trong text

    Args:
        html: HTML response từ submit request.

    Returns:
        Job key string.

    Raises:
        ValueError: Nếu không tìm thấy job key.
    """
    soup = BeautifulSoup(html, "lxml")

    # 1. Hidden input
    job_input = soup.find("input", {"name": "JOB_KEY"})
    if job_input and job_input.get("value"):
        return str(job_input["value"])

    # 2. URL pattern JOB_KEY=...
    match = re.search(r"JOB_KEY=([A-Za-z0-9_\-]+)", html)
    if match:
        return match.group(1)

    # 3. "job id# XXXXX" in text
    match = re.search(r"[Jj]ob\s+[Ii][Dd]#?\s*([A-Za-z0-9_\-]+)", html)
    if match:
        return match.group(1)

    raise ValueError("Không tìm thấy job key trong HTML response từ NCBI Primer-BLAST")


def _detect_status(html: str) -> str:
    """Phân tích HTML để xác định trạng thái job.

    Args:
        html: HTML response từ poll request.

    Returns:
        ``"running"`` | ``"complete"`` | ``"error"``
    """
    html_lower = html.lower()

    for pattern in _ERROR_PATTERNS:
        if pattern in html_lower:
            logger.debug("NCBI Primer-BLAST: error pattern '%s' detected", pattern)
            return "error"

    for pattern in _COMPLETE_PATTERNS:
        if pattern in html_lower:
            logger.debug("NCBI Primer-BLAST: complete pattern '%s' detected", pattern)
            return "complete"

    for pattern in _RUNNING_PATTERNS:
        if pattern in html_lower:
            logger.debug("NCBI Primer-BLAST: running pattern '%s' detected", pattern)
            return "running"

    # Mặc định: coi như đang chạy nếu không rõ
    return "running"


def _parse_primer_table(soup: BeautifulSoup, pair_index: int) -> PrimerPair | None:
    """Parse một cặp primer từ HTML NCBI Primer-BLAST.

    Cấu trúc HTML điển hình:
        <pre>
          Primer pair 1
            Sequence (5'->3')     Start  Stop  Len   Tm    GC%
          Forward primer  ATCG...    10    30   21  60.1   52.4
          Reverse primer  TAGC...   200   220   21  59.8   47.6
          Product length  191
        </pre>

    Args:
        soup: BeautifulSoup object của toàn bộ trang.
        pair_index: Chỉ số cặp primer (0-based) để đặt pair_id.

    Returns:
        :class:`PrimerPair` nếu parse thành công, ``None`` nếu không.
    """
    # Tìm tất cả <pre> blocks chứa primer data
    pre_blocks = soup.find_all("pre")
    target_pre = None
    for pre in pre_blocks:
        text = pre.get_text()
        if re.search(r"[Ff]orward primer", text) and re.search(r"[Rr]everse primer", text):
            # Đếm để chọn đúng block theo pair_index
            if pair_index == 0 or target_pre is None:
                target_pre = pre
                if pair_index == 0:
                    break

    if target_pre is None:
        return None

    text = target_pre.get_text()
    return _parse_primer_text_block(text, pair_index)


def _parse_primer_text_block(text: str, pair_index: int) -> PrimerPair | None:
    """Parse text block chứa thông tin một cặp primer.

    Args:
        text: Text từ <pre> block.
        pair_index: Chỉ số cặp (0-based) để đặt pair_id.

    Returns:
        :class:`PrimerPair` nếu parse đủ dữ liệu, ``None`` nếu thiếu.
    """
    # Pattern: Forward/Left primer   SEQUENCE   start  stop  len   tm   gc%
    fwd_match = re.search(
        r"(?:Forward|Left)\s+primer\s+([ACGTacgt]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)",
        text,
    )
    rev_match = re.search(
        r"(?:Reverse|Right)\s+primer\s+([ACGTacgt]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)",
        text,
    )

    if not fwd_match or not rev_match:
        return None

    left = Oligo(
        sequence=fwd_match.group(1).upper(),
        start=int(fwd_match.group(2)) - 1,  # chuyển về 0-based
        length=int(fwd_match.group(4)),
        tm=float(fwd_match.group(5)),
        gc_percent=float(fwd_match.group(6)),
    )
    right = Oligo(
        sequence=rev_match.group(1).upper(),
        start=int(rev_match.group(2)) - 1,
        length=int(rev_match.group(4)),
        tm=float(rev_match.group(5)),
        gc_percent=float(rev_match.group(6)),
    )

    # Product length
    prod_match = re.search(r"[Pp]roduct\s+(?:length|size)\s+(\d+)", text)
    amplicon_size = int(prod_match.group(1)) if prod_match else 0

    return PrimerPair(
        pair_id=f"ncbi_pair_{pair_index + 1:02d}",
        left_primer=left,
        right_primer=right,
        amplicon_size=amplicon_size,
    )


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

    def submit_job(self, params: dict[str, Any]) -> str:
        """Submit job lên NCBI Primer-BLAST.

        Gộp ``params`` vào tập tham số mặc định rồi POST lên API.
        Trích xuất job key từ HTML response.

        Args:
            params: Dict các tham số gửi lên API (ghi đè default).

        Returns:
            Job key (string) dùng để poll kết quả.

        Raises:
            requests.HTTPError: Nếu HTTP request thất bại.
            ValueError: Nếu không tìm thấy job key trong response.
        """
        merged: dict[str, Any] = {**_DEFAULT_PARAMS, **params}
        merged["CMD"] = "request"

        timeout = self.config.ncbi.timeout
        logger.info("Submitting NCBI Primer-BLAST job (organism=%s)", merged.get("ORGANISM", ""))
        response = requests.post(self.BASE_URL, data=merged, timeout=timeout)
        response.raise_for_status()

        job_key = _extract_job_key(response.text)
        logger.info("NCBI Primer-BLAST job submitted: key=%s", job_key)
        return job_key

    def poll_status(self, job_key: str) -> str:
        """Kiểm tra trạng thái job một lần.

        Args:
            job_key: Job key nhận được từ :meth:`submit_job`.

        Returns:
            ``"running"`` | ``"complete"`` | ``"error"``

        Raises:
            requests.HTTPError: Nếu HTTP request thất bại.
        """
        params = {"CMD": "Get", "JOB_KEY": job_key}
        timeout = self.config.ncbi.timeout
        response = requests.get(self.BASE_URL, params=params, timeout=timeout)
        response.raise_for_status()

        status = _detect_status(response.text)
        logger.debug("NCBI Primer-BLAST poll: key=%s status=%s", job_key, status)
        return status

    def fetch_results(self, job_key: str) -> str:
        """Lấy HTML kết quả khi job đã hoàn thành.

        Args:
            job_key: Job key.

        Returns:
            HTML string của trang kết quả.

        Raises:
            requests.HTTPError: Nếu HTTP request thất bại.
        """
        params = {"CMD": "Get", "JOB_KEY": job_key}
        timeout = self.config.ncbi.timeout
        response = requests.get(self.BASE_URL, params=params, timeout=timeout)
        response.raise_for_status()
        return response.text

    def parse_html_results(self, html: str) -> list[PrimerPair]:
        """Parse HTML kết quả NCBI Primer-BLAST thành List[PrimerPair].

        Tìm tất cả ``<pre>`` block chứa dữ liệu primer pair và parse
        sequence, vị trí, Tm, GC%.

        Args:
            html: HTML string từ NCBI Primer-BLAST kết quả.

        Returns:
            Danh sách :class:`PrimerPair` (có thể rỗng nếu không tìm thấy primer).
        """
        soup = BeautifulSoup(html, "lxml")
        pairs: list[PrimerPair] = []

        # Tìm tất cả <pre> block chứa primer data
        pre_blocks = soup.find_all("pre")
        pair_index = 0
        for pre in pre_blocks:
            text = pre.get_text()
            if re.search(r"(?:Forward|Left)\s+primer", text) and re.search(
                r"(?:Reverse|Right)\s+primer", text
            ):
                # Có thể chứa nhiều cặp primer trong một <pre> block
                # Tách theo "Primer pair N"
                sub_blocks = re.split(r"Primer pair\s+\d+", text)
                if len(sub_blocks) > 1:
                    for block in sub_blocks[1:]:
                        pair = _parse_primer_text_block(block, pair_index)
                        if pair is not None:
                            pairs.append(pair)
                            pair_index += 1
                else:
                    pair = _parse_primer_text_block(text, pair_index)
                    if pair is not None:
                        pairs.append(pair)
                        pair_index += 1

        if not pairs:
            logger.info("NCBI Primer-BLAST: không tìm thấy primer pair trong HTML response")
        else:
            logger.info("NCBI Primer-BLAST: đã parse %d primer pair(s)", len(pairs))
        return pairs

    def _build_params(
        self,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> dict[str, Any]:
        """Xây dựng dict tham số từ DesignTarget và ProcessedSequence.

        Args:
            target: DesignTarget.
            processed_seq: ProcessedSequence.

        Returns:
            Dict tham số sẵn sàng gửi lên API.
        """
        params: dict[str, Any] = {
            "ORGANISM": target.organism,
            "SEQUENCE_TEMPLATE": processed_seq.sequence,
        }

        if target.accession:
            params["INPUT_TYPE"] = "0"
        if target.region_include:
            start, end = target.region_include
            params["SEQUENCE_TARGET"] = f"{start},{end - start}"

        if target.custom_params:
            params.update(target.custom_params)

        return params

    def design_primers(
        self,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> list[PrimerPair]:
        """Thiết kế primer thông qua NCBI Primer-BLAST API.

        Luồng: submit → poll (với exponential backoff) → parse HTML.

        Args:
            target: DesignTarget.
            processed_seq: ProcessedSequence.

        Returns:
            Danh sách :class:`PrimerPair` từ NCBI.

        Raises:
            RuntimeError: Nếu job kết thúc với trạng thái ``"error"``.
            TimeoutError: Nếu job không hoàn thành trong :attr:`MAX_WAIT` giây.
        """
        params = self._build_params(target, processed_seq)
        job_key = self.submit_job(params)

        elapsed = 0.0
        interval = float(self.POLL_INTERVAL)
        while elapsed < self.MAX_WAIT:
            status = self.poll_status(job_key)
            if status == "complete":
                html = self.fetch_results(job_key)
                return self.parse_html_results(html)
            if status == "error":
                raise RuntimeError(
                    f"NCBI Primer-BLAST job {job_key!r} kết thúc với lỗi. "
                    "Kiểm tra sequence hoặc tham số đầu vào."
                )
            logger.info("NCBI Primer-BLAST job %s: đang chờ (elapsed=%.0fs)...", job_key, elapsed)
            time.sleep(interval)
            elapsed += interval
            # Exponential backoff tối đa 60 giây
            interval = min(interval * 1.5, 60.0)

        raise TimeoutError(
            f"NCBI Primer-BLAST job {job_key!r} không hoàn thành sau {self.MAX_WAIT}s"
        )
