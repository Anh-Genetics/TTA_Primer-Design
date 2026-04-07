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

TODO (Sprint 3):
    - Implement ``submit_job()`` với requests
    - Implement ``poll_status()`` với exponential backoff
    - Implement ``parse_html_results()`` với BeautifulSoup
"""

from __future__ import annotations

import logging
from typing import Any

from tta_primer_design.config import AppConfig
from tta_primer_design.models import DesignTarget, PrimerPair, ProcessedSequence

logger = logging.getLogger("tta_primer_design.modules.ncbi_primer_blast")

_BASE_URL = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"
_POLL_INTERVAL = 10  # giây
_MAX_WAIT = 300  # giây tối đa chờ kết quả


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

        Args:
            params: Dict các tham số gửi lên API.

        Returns:
            Job key (string) dùng để poll kết quả.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("submit_job chưa được implement (Sprint 3)")

    def poll_status(self, job_key: str) -> str:
        """Kiểm tra trạng thái job.

        Args:
            job_key: Job key nhận được từ submit_job.

        Returns:
            "running" | "complete" | "error"

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("poll_status chưa được implement (Sprint 3)")

    def fetch_results(self, job_key: str) -> str:
        """Lấy kết quả HTML khi job hoàn thành.

        Args:
            job_key: Job key.

        Returns:
            HTML string của trang kết quả.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("fetch_results chưa được implement (Sprint 3)")

    def parse_html_results(self, html: str) -> list[PrimerPair]:
        """Parse HTML kết quả thành List[PrimerPair].

        Args:
            html: HTML string từ NCBI Primer-BLAST.

        Returns:
            Danh sách PrimerPair.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("parse_html_results chưa được implement (Sprint 3)")

    def design_primers(
        self,
        target: DesignTarget,
        processed_seq: ProcessedSequence,
    ) -> list[PrimerPair]:
        """Thiết kế primer thông qua NCBI Primer-BLAST API.

        Wrapper tổng hợp: submit → poll → parse.

        Args:
            target: DesignTarget.
            processed_seq: ProcessedSequence.

        Returns:
            Danh sách PrimerPair từ NCBI.

        Raises:
            NotImplementedError: Chưa implement (Sprint 3).
        """
        raise NotImplementedError("design_primers chưa được implement (Sprint 3)")
