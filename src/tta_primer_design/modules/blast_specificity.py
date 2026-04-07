"""Module 6 — BlastSpecificity: kiểm tra tính đặc hiệu của primer/probe.

Ba phương án:
    A) Dùng kết quả sẵn từ NCBI Primer-BLAST (ưu tiên)
    B) Chạy BLAST+ local cho batch lớn
    C) Gọi NCBI BLAST API cho từng oligo

Tiêu chí đánh giá:
    - Số off-target hits có thể tạo amplicon
    - Max amplicon size từ off-target (nếu > cutoff → FAIL)
    - Perfect match (mismatch = 0) trên non-target sequences
    - 3' end mismatch (≥1 mismatch ở 3 bp cuối = acceptable)
"""

from __future__ import annotations

import logging
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import requests

from tta_primer_design.config import AppConfig
from tta_primer_design.models import PrimerPair

logger = logging.getLogger("tta_primer_design.modules.blast_specificity")

_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
_BLAST_POLL_INTERVAL = 10  # seconds between status polls
_BLAST_MAX_WAIT = 180  # max seconds to wait per BLAST job
_MAX_OFFTARGET_AMPLICON_SIZE = 3000  # bp

# Specificity scoring constants
_OFFTARGET_PENALTY_PER_AMPLICON = 25.0  # score penalty per off-target amplicon
_MIN_SCORE_WITH_HITS = 80.0  # minimum score floor when many non-specific hits
_HIT_THRESHOLD = 10  # total hit count above which minor penalty applies
_HIT_PENALTY = 0.5  # score penalty per extra hit above threshold


@dataclass
class BlastHit:
    """Một BLAST hit cho một oligo.

    Attributes:
        subject_id: ID sequence đích trong database.
        subject_title: Tên/mô tả sequence đích.
        identity: % identity.
        alignment_length: Độ dài alignment.
        mismatches: Số mismatch.
        gaps: Số gap.
        query_start: Vị trí bắt đầu trên query.
        query_end: Vị trí kết thúc trên query.
        subject_start: Vị trí bắt đầu trên subject.
        subject_end: Vị trí kết thúc trên subject.
        evalue: E-value.
        bit_score: Bit score.
        mismatches_3prime: Số mismatch ở 3 bp cuối primer.
    """

    subject_id: str
    subject_title: str = ""
    identity: float = 0.0
    alignment_length: int = 0
    mismatches: int = 0
    gaps: int = 0
    query_start: int = 0
    query_end: int = 0
    subject_start: int = 0
    subject_end: int = 0
    evalue: float = 0.0
    bit_score: float = 0.0
    mismatches_3prime: int = 0


@dataclass
class OffTargetAmplicon:
    """Một amplicon off-target tiềm năng (2 primer hit gần nhau).

    Attributes:
        subject_id: ID sequence đích.
        amplicon_size: Kích thước amplicon ước tính.
        left_hit: BLAST hit của primer trái.
        right_hit: BLAST hit của primer phải.
    """

    subject_id: str
    amplicon_size: int
    left_hit: BlastHit
    right_hit: BlastHit


@dataclass
class SpecificityResult:
    """Kết quả kiểm tra specificity cho một PrimerPair.

    Attributes:
        primer_pair_id: ID cặp primer.
        is_specific: True nếu pass specificity check.
        off_target_amplicons: Danh sách amplicon off-target.
        specificity_score: Điểm specificity (0–100).
        blast_hits_left: BLAST hits của primer trái.
        blast_hits_right: BLAST hits của primer phải.
        blast_hits_probe: BLAST hits của probe (nếu có).
    """

    primer_pair_id: str
    is_specific: bool = True
    off_target_amplicons: list[OffTargetAmplicon] = field(default_factory=list)
    specificity_score: float = 100.0
    blast_hits_left: list[BlastHit] = field(default_factory=list)
    blast_hits_right: list[BlastHit] = field(default_factory=list)
    blast_hits_probe: list[BlastHit] | None = None


class BlastSpecificity:
    """Kiểm tra specificity của primer/probe bằng NCBI BLAST API.

    Args:
        config: AppConfig đã load.

    Example::

        checker = BlastSpecificity(config)
        results = checker.check_all(primer_pairs, organism="Homo sapiens")
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config
        self._session = requests.Session()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def check_all(
        self,
        primer_pairs: list[PrimerPair],
        organism: str = "Homo sapiens",
    ) -> list[PrimerPair]:
        """Kiểm tra specificity cho tất cả primer pairs.

        Gán ``specificity_result`` vào từng PrimerPair.

        Args:
            primer_pairs: Danh sách PrimerPair cần kiểm tra.
            organism: Tên loài.

        Returns:
            Danh sách PrimerPair đã có specificity_result.
        """
        for pair in primer_pairs:
            try:
                result = self.check_pair(pair, organism=organism)
                pair.specificity_result = result
            except Exception as exc:
                logger.warning("Specificity check failed for pair %s: %s", pair.pair_id, exc)
                pair.specificity_result = SpecificityResult(
                    primer_pair_id=pair.pair_id,
                    is_specific=False,
                    specificity_score=0.0,
                )
        return primer_pairs

    def check_pair(
        self,
        pair: PrimerPair,
        organism: str = "Homo sapiens",
    ) -> SpecificityResult:
        """Kiểm tra specificity cho một PrimerPair.

        BLAST left và right primer; tìm off-target amplicons; tính score.

        Args:
            pair: PrimerPair cần kiểm tra.
            organism: Tên loài.

        Returns:
            SpecificityResult.
        """
        database = self.config.blast.database
        evalue = self.config.blast.evalue
        word_size = self.config.blast.word_size

        logger.info("BLASTing pair %s", pair.pair_id)

        left_hits = self._blast_sequence(
            pair.left_primer.sequence,
            database=database,
            organism=organism,
            evalue=evalue,
            word_size=word_size,
        )
        right_hits = self._blast_sequence(
            pair.right_primer.sequence,
            database=database,
            organism=organism,
            evalue=evalue,
            word_size=word_size,
        )

        probe_hits: list[BlastHit] | None = None
        if pair.probe is not None:
            probe_hits = self._blast_sequence(
                pair.probe.sequence,
                database=database,
                organism=organism,
                evalue=evalue,
                word_size=word_size,
            )

        off_targets = self._find_off_target_amplicons(left_hits, right_hits)
        score = self._calculate_specificity_score(off_targets, len(left_hits), len(right_hits))

        return SpecificityResult(
            primer_pair_id=pair.pair_id,
            is_specific=len(off_targets) == 0,
            off_target_amplicons=off_targets,
            specificity_score=score,
            blast_hits_left=left_hits,
            blast_hits_right=right_hits,
            blast_hits_probe=probe_hits,
        )

    # ------------------------------------------------------------------
    # Private helpers — NCBI BLAST API
    # ------------------------------------------------------------------

    def _blast_sequence(
        self,
        sequence: str,
        database: str = "refseq_rna",
        organism: str = "Homo sapiens",
        evalue: float = 1000.0,
        word_size: int = 7,
    ) -> list[BlastHit]:
        """Gửi sequence lên NCBI BLAST và trả về danh sách BlastHit.

        Args:
            sequence: Chuỗi oligo cần BLAST.
            database: BLAST database.
            organism: Tên loài lọc.
            evalue: E-value threshold.
            word_size: Word size cho blastn.

        Returns:
            Danh sách BlastHit (có thể rỗng nếu gặp lỗi).
        """
        rid = self._submit_blast(sequence, database, organism, evalue, word_size)
        if rid is None:
            return []

        status = self._poll_blast(rid)
        if status != "READY":
            logger.warning("BLAST job %s did not complete (status=%s)", rid, status)
            return []

        xml_data = self._fetch_blast_results(rid)
        if not xml_data:
            return []

        return self._parse_blast_xml(xml_data, len(sequence))

    def _submit_blast(
        self,
        sequence: str,
        database: str,
        organism: str,
        evalue: float,
        word_size: int,
    ) -> str | None:
        """Submit BLAST job; trả về RID (Request ID) hoặc None nếu lỗi."""
        params: dict[str, str | float | int] = {
            "CMD": "Put",
            "PROGRAM": "blastn",
            "DATABASE": database,
            "QUERY": sequence,
            "ENTREZ_QUERY": f'"{organism}"[organism]',
            "EXPECT": evalue,
            "WORD_SIZE": word_size,
            "FORMAT_TYPE": "XML",
            "HITLIST_SIZE": self.config.blast.max_targets,
            "FILTER": "L",  # low-complexity filter
        }
        if self.config.ncbi.email:
            params["email"] = self.config.ncbi.email

        try:
            resp = self._session.post(_BLAST_URL, data=params, timeout=self.config.ncbi.timeout)
            resp.raise_for_status()
        except requests.RequestException as exc:
            logger.error("BLAST submit failed: %s", exc)
            return None

        # Extract RID from response HTML
        for line in resp.text.splitlines():
            if "RID" in line and "=" in line:
                parts = line.strip().split("=")
                if len(parts) == 2 and parts[0].strip() == "RID":
                    return parts[1].strip()

        logger.warning("Could not extract RID from BLAST response")
        return None

    def _poll_blast(self, rid: str) -> str:
        """Poll BLAST job until READY or FAILED.

        Args:
            rid: BLAST Request ID.

        Returns:
            "READY" | "FAILED" | "WAITING"
        """
        elapsed = 0
        interval = _BLAST_POLL_INTERVAL

        while elapsed < _BLAST_MAX_WAIT:
            time.sleep(interval)
            elapsed += interval

            try:
                resp = self._session.get(
                    _BLAST_URL,
                    params={
                        "CMD": "Get",
                        "FORMAT_OBJECT": "SearchInfo",
                        "RID": rid,
                    },
                    timeout=self.config.ncbi.timeout,
                )
                resp.raise_for_status()
            except requests.RequestException as exc:
                logger.warning("BLAST poll error: %s", exc)
                continue

            text = resp.text
            if "Status=READY" in text:
                return "READY"
            if "Status=FAILED" in text or "Status=UNKNOWN" in text:
                return "FAILED"

            logger.debug("BLAST job %s still waiting (%ds elapsed)", rid, elapsed)

        return "WAITING"

    def _fetch_blast_results(self, rid: str) -> str:
        """Lấy kết quả BLAST XML cho RID.

        Args:
            rid: BLAST Request ID.

        Returns:
            XML string hoặc rỗng nếu lỗi.
        """
        try:
            resp = self._session.get(
                _BLAST_URL,
                params={
                    "CMD": "Get",
                    "FORMAT_TYPE": "XML",
                    "RID": rid,
                },
                timeout=self.config.ncbi.timeout,
            )
            resp.raise_for_status()
            return resp.text
        except requests.RequestException as exc:
            logger.error("BLAST fetch results failed: %s", exc)
            return ""

    def _parse_blast_xml(self, xml_data: str, query_len: int) -> list[BlastHit]:
        """Parse BLAST XML output thành danh sách BlastHit.

        Args:
            xml_data: BLAST XML string.
            query_len: Độ dài query (để tính 3' mismatch).

        Returns:
            Danh sách BlastHit.
        """
        hits: list[BlastHit] = []
        try:
            root = ET.fromstring(xml_data)
        except ET.ParseError as exc:
            logger.error("Failed to parse BLAST XML: %s", exc)
            return hits

        # Support both old BlastOutput and new BlastXML2 formats
        # Old format: BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_hsps/Hsp
        for hit_elem in root.iter("Hit"):
            subject_id = _xml_text(hit_elem, "Hit_id") or _xml_text(hit_elem, "Hit_accession", "")
            subject_title = _xml_text(hit_elem, "Hit_def", "")

            for hsp in hit_elem.iter("Hsp"):
                try:
                    ident = int(_xml_text(hsp, "Hsp_identity", "0"))
                    align_len = int(_xml_text(hsp, "Hsp_align-len", "0"))
                    mismatches = int(_xml_text(hsp, "Hsp_mismatch", "0"))
                    gaps = int(_xml_text(hsp, "Hsp_gaps", "0"))
                    q_start = int(_xml_text(hsp, "Hsp_query-from", "0"))
                    q_end = int(_xml_text(hsp, "Hsp_query-to", "0"))
                    s_start = int(_xml_text(hsp, "Hsp_hit-from", "0"))
                    s_end = int(_xml_text(hsp, "Hsp_hit-to", "0"))
                    evalue = float(_xml_text(hsp, "Hsp_evalue", "1.0"))
                    bit_score = float(_xml_text(hsp, "Hsp_bit-score", "0.0"))
                    identity_pct = (ident / align_len * 100.0) if align_len else 0.0

                    # Estimate 3' mismatches from midline string
                    midline = _xml_text(hsp, "Hsp_midline", "")
                    mm_3prime = _count_3prime_mismatches(midline, s_start, s_end)

                    blast_hit = BlastHit(
                        subject_id=subject_id or "",
                        subject_title=subject_title,
                        identity=identity_pct,
                        alignment_length=align_len,
                        mismatches=mismatches,
                        gaps=gaps,
                        query_start=q_start,
                        query_end=q_end,
                        subject_start=s_start,
                        subject_end=s_end,
                        evalue=evalue,
                        bit_score=bit_score,
                        mismatches_3prime=mm_3prime,
                    )
                    hits.append(blast_hit)
                except (ValueError, TypeError) as exc:
                    logger.debug("Skipping malformed HSP: %s", exc)

        return hits

    # ------------------------------------------------------------------
    # Private helpers — off-target amplicon detection
    # ------------------------------------------------------------------

    def _find_off_target_amplicons(
        self,
        left_hits: list[BlastHit],
        right_hits: list[BlastHit],
        max_amplicon_size: int = _MAX_OFFTARGET_AMPLICON_SIZE,
    ) -> list[OffTargetAmplicon]:
        """Phát hiện amplicon off-target từ các BLAST hits.

        Một amplicon off-target xảy ra khi left primer và right primer đều
        hit cùng một subject sequence trong khoảng cách hợp lý và đúng hướng:
            - Left primer hit plus strand (s_start < s_end)
            - Right primer hit minus strand (s_start > s_end)
            - Khoảng cách giữa hai hit < max_amplicon_size

        Args:
            left_hits: BLAST hits của primer trái.
            right_hits: BLAST hits của primer phải.
            max_amplicon_size: Kích thước amplicon off-target tối đa.

        Returns:
            Danh sách OffTargetAmplicon tiềm năng.
        """
        off_targets: list[OffTargetAmplicon] = []

        # Index right hits by subject_id for fast lookup
        right_by_subject: dict[str, list[BlastHit]] = {}
        for rh in right_hits:
            right_by_subject.setdefault(rh.subject_id, []).append(rh)

        for lh in left_hits:
            # Left primer should hit on plus strand
            if lh.subject_start >= lh.subject_end:
                continue
            # Ignore hits with 3' mismatches (primer won't extend)
            if lh.mismatches_3prime > 0:
                continue

            rh_list = right_by_subject.get(lh.subject_id, [])
            for rh in rh_list:
                # Right primer should hit on minus strand (coordinates reversed)
                if rh.subject_start <= rh.subject_end:
                    continue
                if rh.mismatches_3prime > 0:
                    continue

                # Amplicon: left_start to right_start (minus strand end = 5' of right primer)
                amp_start = lh.subject_start
                amp_end = rh.subject_start  # highest coordinate of right hit
                amplicon_size = amp_end - amp_start

                if 0 < amplicon_size <= max_amplicon_size:
                    off_targets.append(
                        OffTargetAmplicon(
                            subject_id=lh.subject_id,
                            amplicon_size=amplicon_size,
                            left_hit=lh,
                            right_hit=rh,
                        )
                    )

        return off_targets

    def _calculate_specificity_score(
        self,
        off_targets: list[OffTargetAmplicon],
        n_left_hits: int,
        n_right_hits: int,
    ) -> float:
        """Tính điểm specificity (0–100).

        Score = 100 nếu không có off-target amplicon.
        Giảm dần theo số off-target và số BLAST hits.

        Args:
            off_targets: Danh sách off-target amplicons.
            n_left_hits: Số hits của primer trái.
            n_right_hits: Số hits của primer phải.

        Returns:
            Điểm specificity (0.0–100.0).
        """
        if off_targets:
            # Each off-target amplicon reduces score significantly
            penalty = min(len(off_targets) * _OFFTARGET_PENALTY_PER_AMPLICON, 100.0)
            return max(0.0, 100.0 - penalty)

        # Minor penalty for many non-specific hits
        total_hits = n_left_hits + n_right_hits
        if total_hits > _HIT_THRESHOLD:
            return max(_MIN_SCORE_WITH_HITS, 100.0 - (total_hits - _HIT_THRESHOLD) * _HIT_PENALTY)

        return 100.0


# ------------------------------------------------------------------
# Module-level helpers
# ------------------------------------------------------------------


def _xml_text(elem: ET.Element, tag: str, default: str = "") -> str:
    """Extract text from a child XML element."""
    child = elem.find(tag)
    return child.text.strip() if child is not None and child.text else default


def _count_3prime_mismatches(midline: str, s_start: int, s_end: int) -> int:
    """Đếm mismatch ở 3 bp cuối primer từ midline string.

    Args:
        midline: Alignment midline (spaces = mismatch/gap, | or letter = match).
        s_start: Subject start coordinate.
        s_end: Subject end coordinate.

    Returns:
        Số mismatch ở 3 bp cuối.
    """
    if not midline:
        return 0
    # Plus strand: 3' end is at right side of alignment
    # Minus strand (s_start > s_end): 3' end is at left side
    tail = midline[-3:] if s_start < s_end else midline[:3]
    return sum(1 for c in tail if c == " ")
