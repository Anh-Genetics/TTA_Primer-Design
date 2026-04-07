"""Module 6 — BlastSpecificity: kiểm tra tính đặc hiệu của primer/probe.

Ba phương án:
    A) Dùng kết quả sẵn từ NCBI Primer-BLAST (ưu tiên)
    B) Chạy BLAST+ local cho batch lớn
    C) Gọi NCBI BLAST API cho từng oligo ← ĐÃ IMPLEMENT

Tiêu chí đánh giá:
    - Số off-target hits có thể tạo amplicon
    - Max amplicon size từ off-target (nếu > cutoff → FAIL)
    - Perfect match (mismatch = 0) trên non-target sequences
    - 3' end mismatch (≥1 mismatch ở 3 bp cuối = acceptable)

Sử dụng biopython NCBIWWW để gọi NCBI BLAST Web API.
Yêu cầu kết nối internet và email NCBI hợp lệ.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field

from tta_primer_design.config import AppConfig
from tta_primer_design.models import PrimerPair

logger = logging.getLogger("tta_primer_design.modules.blast_specificity")


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


def _calc_3prime_mismatches(hsp: object, query_len: int, prime_len: int = 3) -> int:
    """Đếm số mismatch trong ``prime_len`` base cuối của query trong một HSP.

    Args:
        hsp: Bio.Blast HSP object (có .query, .sbjct, .query_end).
        query_len: Độ dài đầy đủ của oligo query.
        prime_len: Số base 3'-end cần kiểm tra (default 3).

    Returns:
        Số mismatch ở đuôi 3' của query trong alignment này.
    """
    if hsp.query_end < query_len:
        return 0  # 3' end không nằm trong alignment này

    q_align: str = hsp.query
    s_align: str = hsp.sbjct
    count_from_3prime = 0
    mismatches = 0

    for i in range(len(q_align) - 1, -1, -1):
        q_base = q_align[i].upper()
        s_base = s_align[i].upper()
        if q_base != "-":
            count_from_3prime += 1
            if count_from_3prime <= prime_len:
                if s_base != "-" and q_base != s_base:
                    mismatches += 1
        if count_from_3prime >= prime_len:
            break

    return mismatches


def _blast_oligo_ncbi(
    sequence: str,
    database: str,
    organism: str,
    email: str,
    evalue: float = 1000.0,
    word_size: int = 7,
    hitlist_size: int = 50,
) -> list[BlastHit]:
    """BLAST một oligo đơn bằng NCBI BLAST Web API (biopython).

    Sử dụng ``Bio.Blast.NCBIWWW.qblast()`` với tham số tối ưu cho primer ngắn:
    ``megablast=False``, ``filter=None`` (tắt DUST), ``word_size=7``.

    Args:
        sequence: Chuỗi oligo (5'→3').
        database: NCBI BLAST database (ví dụ: ``"nt"``, ``"refseq_rna"``).
        organism: Tên loài để lọc kết quả (ví dụ: ``"Homo sapiens"``).
        email: Email NCBI (bắt buộc theo NCBI policy).
        evalue: E-value threshold (default 1000 để bắt các hit yếu).
        word_size: Kích thước word BLAST (7 phù hợp cho primer 18–30 bp).
        hitlist_size: Số hits tối đa trả về.

    Returns:
        Danh sách :class:`BlastHit`.

    Raises:
        ImportError: Nếu biopython chưa được cài.
        RuntimeError: Nếu gặp lỗi khi gọi API hoặc parse kết quả.
    """
    try:
        from Bio import Entrez
        from Bio.Blast import NCBIWWW, NCBIXML
    except ImportError as exc:
        raise ImportError(
            "biopython is required for NCBI BLAST. Run: pip install biopython"
        ) from exc

    Entrez.email = email
    entrez_query = f"{organism}[Organism]" if organism else ""

    logger.info("BLASTing oligo %s... against %s (organism=%s)", sequence[:12], database, organism)

    try:
        handle = NCBIWWW.qblast(
            "blastn",
            database,
            sequence,
            entrez_query=entrez_query,
            expect=evalue,
            word_size=word_size,
            hitlist_size=hitlist_size,
            megablast=False,
            filter=None,  # tắt DUST filter để giữ lại hits ngắn
        )
        blast_record = NCBIXML.read(handle)
        handle.close()
    except Exception as exc:
        raise RuntimeError(f"NCBI BLAST API error for sequence {sequence[:12]!r}: {exc}") from exc

    query_len = len(sequence)
    hits: list[BlastHit] = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            mismatches = hsp.align_length - hsp.identities - hsp.gaps
            mismatches_3p = _calc_3prime_mismatches(hsp, query_len)
            hits.append(
                BlastHit(
                    subject_id=alignment.accession,
                    subject_title=alignment.title[:120],
                    identity=hsp.identities / hsp.align_length * 100.0,
                    alignment_length=hsp.align_length,
                    mismatches=mismatches,
                    gaps=hsp.gaps,
                    query_start=hsp.query_start,
                    query_end=hsp.query_end,
                    subject_start=hsp.sbjct_start,
                    subject_end=hsp.sbjct_end,
                    evalue=hsp.expect,
                    bit_score=hsp.bits,
                    mismatches_3prime=mismatches_3p,
                )
            )

    logger.info("BLAST returned %d hit(s) for %s...", len(hits), sequence[:12])
    return hits


def _find_off_target_amplicons(
    left_hits: list[BlastHit],
    right_hits: list[BlastHit],
    max_size: int = 2000,
    max_mismatches: int = 2,
) -> list[OffTargetAmplicon]:
    """Tìm các off-target amplicon tiềm năng từ BLAST hits của 2 primer.

    Logic:
        - Left primer: hit trên forward strand (``sbjct_start < sbjct_end``)
        - Right primer: hit trên reverse strand (``sbjct_start > sbjct_end``)
        - Cùng subject ID, khoảng cách ≤ ``max_size``

    Args:
        left_hits: BLAST hits của primer trái.
        right_hits: BLAST hits của primer phải.
        max_size: Kích thước amplicon tối đa (bp) để coi là off-target.
        max_mismatches: Số mismatch tối đa cho một hit để xét vào phân tích.

    Returns:
        Danh sách :class:`OffTargetAmplicon` tiềm năng.
    """
    left_fwd: dict[str, list[BlastHit]] = defaultdict(list)
    right_rev: dict[str, list[BlastHit]] = defaultdict(list)

    for hit in left_hits:
        # Forward strand: sbjct_start < sbjct_end
        if hit.subject_start < hit.subject_end and hit.mismatches <= max_mismatches:
            left_fwd[hit.subject_id].append(hit)

    for hit in right_hits:
        # Reverse strand: BLAST reports sbjct_start > sbjct_end cho minus-strand hit
        if hit.subject_start > hit.subject_end and hit.mismatches <= max_mismatches:
            right_rev[hit.subject_id].append(hit)

    off_targets: list[OffTargetAmplicon] = []
    for subject_id in set(left_fwd) & set(right_rev):
        for lhit in left_fwd[subject_id]:
            for rhit in right_rev[subject_id]:
                left_pos = lhit.subject_start  # 5' end of left primer (lower coord)
                # On reverse strand, subject_start > subject_end; subject_start is
                # the 5' end of the right primer in genomic coordinates (higher coord).
                right_pos = rhit.subject_start
                if right_pos > left_pos:
                    amp_size = right_pos - left_pos + 1
                    if amp_size <= max_size:
                        off_targets.append(
                            OffTargetAmplicon(
                                subject_id=subject_id,
                                amplicon_size=amp_size,
                                left_hit=lhit,
                                right_hit=rhit,
                            )
                        )

    return off_targets


class BlastSpecificity:
    """Kiểm tra specificity của primer/probe bằng NCBI BLAST Web API.

    Sử dụng :func:`_blast_oligo_ncbi` để BLAST từng oligo, sau đó
    :func:`_find_off_target_amplicons` để phát hiện amplicon off-target tiềm năng.

    Args:
        config: AppConfig đã load.

    Example::

        checker = BlastSpecificity(config)
        result = checker.check_pair(pair, organism="Homo sapiens")
        pairs  = checker.check_all(primer_pairs, organism="Homo sapiens")
    """

    def __init__(self, config: AppConfig) -> None:
        self.config = config

    def check_all(
        self,
        primer_pairs: list[PrimerPair],
        organism: str = "Homo sapiens",
    ) -> list[PrimerPair]:
        """Kiểm tra specificity cho tất cả primer pairs.

        Gán ``specificity_result`` vào từng :class:`PrimerPair`.
        Các pair bị lỗi BLAST sẽ được log warning và bỏ qua.

        Args:
            primer_pairs: Danh sách PrimerPair cần kiểm tra.
            organism: Tên loài (ví dụ: ``"Homo sapiens"``).

        Returns:
            Danh sách PrimerPair đã có ``specificity_result`` (nếu thành công).
        """
        for pair in primer_pairs:
            try:
                result = self.check_pair(pair, organism=organism)
                pair.specificity_result = result
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "BLAST failed for pair '%s' — skipping specificity check: %s",
                    pair.pair_id,
                    exc,
                )
        return primer_pairs

    def check_pair(
        self,
        pair: PrimerPair,
        organism: str = "Homo sapiens",
    ) -> SpecificityResult:
        """Kiểm tra specificity cho một PrimerPair bằng NCBI BLAST.

        Thực hiện:
            1. BLAST left primer
            2. BLAST right primer
            3. (Nếu có probe) BLAST probe
            4. Tìm off-target amplicons từ hits của left + right
            5. Tính specificity_score

        Args:
            pair: PrimerPair cần kiểm tra.
            organism: Tên loài lọc BLAST.

        Returns:
            :class:`SpecificityResult` đầy đủ.

        Raises:
            RuntimeError: Nếu NCBI BLAST API thất bại.
        """
        db = self.config.blast.database
        email = self.config.ncbi.email
        evalue = self.config.blast.evalue
        word_size = self.config.blast.word_size

        left_hits = _blast_oligo_ncbi(
            pair.left_primer.sequence,
            db,
            organism,
            email,
            evalue=evalue,
            word_size=word_size,
        )
        right_hits = _blast_oligo_ncbi(
            pair.right_primer.sequence,
            db,
            organism,
            email,
            evalue=evalue,
            word_size=word_size,
        )
        probe_hits: list[BlastHit] | None = None
        if pair.probe is not None:
            probe_hits = _blast_oligo_ncbi(
                pair.probe.sequence,
                db,
                organism,
                email,
                evalue=evalue,
                word_size=word_size,
            )

        off_targets = _find_off_target_amplicons(left_hits, right_hits)

        is_specific = len(off_targets) <= self.config.filters.max_off_targets
        # Trừ 20 điểm mỗi off-target amplicon, tối thiểu 0
        specificity_score = max(0.0, 100.0 - len(off_targets) * 20.0)

        return SpecificityResult(
            primer_pair_id=pair.pair_id,
            is_specific=is_specific,
            off_target_amplicons=off_targets,
            specificity_score=specificity_score,
            blast_hits_left=left_hits,
            blast_hits_right=right_hits,
            blast_hits_probe=probe_hits,
        )
