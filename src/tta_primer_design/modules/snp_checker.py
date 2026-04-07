"""Module 8 — SNPChecker: kiểm tra SNP trong vùng primer/probe.

Nguồn dữ liệu:
    NCBI Variation Services API (online)
    Endpoint: https://api.ncbi.nlm.nih.gov/variation/v0/

Logic:
    1. Lấy genomic coordinates của primer/probe từ SeqRecord (nếu có) hoặc
       dùng vị trí tương đối trên template sequence.
    2. Query NCBI dbSNP cho region đó qua eUtils esearch + esummary.
    3. Flag SNPs trong primer với MAF > threshold (default: 0.01).
    4. Cảnh báo đặc biệt SNP ở 3' end primer (3 bp cuối).
    5. Cảnh báo tất cả SNP trong probe.

Tiêu chí đánh giá:
    - "PASS"    : không có SNP nào có MAF >= threshold
    - "WARNING" : có SNP nhưng MAF < threshold hoặc ở giữa primer (không ở 3')
    - "FAIL"    : SNP có MAF >= threshold ở 3' end primer, hoặc bất kỳ SNP
                  nào trong probe
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import requests

from tta_primer_design.config import AppConfig
from tta_primer_design.models import PrimerPair

logger = logging.getLogger("tta_primer_design.modules.snp_checker")

_DEFAULT_MAF_THRESHOLD = 0.01
_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_ESEARCH_URL = f"{_EUTILS_BASE}/esearch.fcgi"
_ESUMMARY_URL = f"{_EUTILS_BASE}/esummary.fcgi"
_REQUEST_TIMEOUT = 30  # giây


@dataclass
class SNPInfo:
    """Thông tin về một SNP trong vùng oligo.

    Attributes:
        rsid: ID dbSNP (ví dụ: rs12345).
        position_in_oligo: Vị trí SNP trong oligo (0-based từ 5').
        maf: Minor Allele Frequency (0.0–1.0).
        alleles: Allele string (ví dụ: "A/G").
        is_3prime: True nếu SNP nằm ở 3 bp cuối primer.
    """

    rsid: str
    position_in_oligo: int = 0
    maf: float = 0.0
    alleles: str = ""
    is_3prime: bool = False


@dataclass
class SNPResult:
    """Kết quả kiểm tra SNP cho một oligo.

    Attributes:
        oligo_id: ID oligo (primer/probe).
        has_snp: True nếu có bất kỳ SNP nào.
        snp_list: Danh sách tất cả SNP.
        critical_snps: SNP nghiêm trọng (3' end hoặc trong probe).
        recommendation: "PASS" | "WARNING" | "FAIL".
    """

    oligo_id: str
    has_snp: bool = False
    snp_list: list[SNPInfo] = field(default_factory=list)
    critical_snps: list[SNPInfo] = field(default_factory=list)
    recommendation: str = "PASS"


def _search_snps_by_accession_range(
    accession: str,
    start: int,
    end: int,
    email: str,
    api_key: str = "",
    retmax: int = 100,
) -> list[str]:
    """Tìm SNP IDs (rsIDs) trong vùng [start, end] của một accession.

    Dùng NCBI eSearch với truy vấn:
    ``{accession}[SeqID] AND {start}:{end}[Base Position]``

    Args:
        accession: GenBank/RefSeq accession (ví dụ: ``"NM_001101"``).
        start: Vị trí bắt đầu (1-based, inclusive).
        end: Vị trí kết thúc (1-based, inclusive).
        email: Email NCBI (bắt buộc theo NCBI policy).
        api_key: NCBI API key (tuỳ chọn, tăng rate limit).
        retmax: Số kết quả tối đa trả về.

    Returns:
        Danh sách rsID strings (ví dụ: ``["12345", "67890"]``).

    Raises:
        requests.RequestException: Nếu HTTP request thất bại.
        ValueError: Nếu response không parse được.
    """
    query = f"{accession}[SeqID] AND {start}:{end}[Base Position]"
    params: dict[str, Any] = {
        "db": "snp",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key

    logger.debug("dbSNP esearch: %s [%d-%d]", accession, start, end)
    response = requests.get(_ESEARCH_URL, params=params, timeout=_REQUEST_TIMEOUT)
    response.raise_for_status()

    data = response.json()
    id_list: list[str] = data.get("esearchresult", {}).get("idlist", [])
    logger.debug("dbSNP esearch returned %d SNP(s)", len(id_list))
    return id_list


def _fetch_snp_summaries(
    snp_ids: list[str],
    email: str,
    api_key: str = "",
) -> list[dict[str, Any]]:
    """Lấy thông tin chi tiết của các SNP từ esummary.

    Args:
        snp_ids: Danh sách rsID (dạng string, không có tiền tố ``"rs"``).
        email: Email NCBI.
        api_key: NCBI API key.

    Returns:
        Danh sách dict, mỗi dict chứa thông tin một SNP (từ esummary).
        Trả về ``[]`` nếu ``snp_ids`` rỗng.

    Raises:
        requests.RequestException: Nếu HTTP request thất bại.
    """
    if not snp_ids:
        return []

    params: dict[str, Any] = {
        "db": "snp",
        "id": ",".join(snp_ids),
        "retmode": "json",
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key

    response = requests.get(_ESUMMARY_URL, params=params, timeout=_REQUEST_TIMEOUT)
    response.raise_for_status()

    data = response.json()
    result_dict: dict[str, Any] = data.get("result", {})
    uids: list[str] = result_dict.get("uids", [])
    return [result_dict[uid] for uid in uids if uid in result_dict]


def _extract_maf(summary: dict[str, Any]) -> float:
    """Trích xuất MAF từ esummary dict của một SNP.

    Thử lấy từ field ``"maf"`` hoặc ``"global_mafs"`` (tuỳ phiên bản API).

    Args:
        summary: Dict esummary của một SNP.

    Returns:
        MAF value (0.0–1.0). Trả về 0.0 nếu không có dữ liệu.
    """
    # Trường hợp 1: field "maf" trực tiếp
    maf_val = summary.get("maf")
    if maf_val is not None:
        try:
            return float(maf_val)
        except (TypeError, ValueError):
            pass

    # Trường hợp 2: "global_mafs" list
    global_mafs = summary.get("global_mafs")
    if isinstance(global_mafs, list) and global_mafs:
        for entry in global_mafs:
            freq_str = entry.get("freq")
            if freq_str:
                try:
                    return float(freq_str)
                except (TypeError, ValueError):
                    continue

    return 0.0


def _extract_alleles(summary: dict[str, Any]) -> str:
    """Trích xuất allele string từ esummary dict.

    Args:
        summary: Dict esummary của một SNP.

    Returns:
        Allele string (ví dụ: ``"A/G"``). Trả về ``""`` nếu không có.
    """
    return str(summary.get("allele_origin", "") or summary.get("docsum", ""))


def _build_snp_result(
    oligo_id: str,
    oligo_length: int,
    snp_summaries: list[dict[str, Any]],
    genomic_start: int,
    oligo_genomic_start: int,
    maf_threshold: float,
    is_probe: bool,
) -> SNPResult:
    """Xây dựng SNPResult từ danh sách SNP summaries.

    Args:
        oligo_id: ID oligo.
        oligo_length: Độ dài oligo.
        snp_summaries: Danh sách dict từ esummary.
        genomic_start: Vị trí genomic bắt đầu của oligo (1-based).
        oligo_genomic_start: Vị trí genomic bắt đầu của oligo (1-based, dùng để tính offset).
        maf_threshold: Ngưỡng MAF.
        is_probe: True nếu là probe.

    Returns:
        :class:`SNPResult` đã được phân loại.
    """
    snp_list: list[SNPInfo] = []
    critical_snps: list[SNPInfo] = []

    for summary in snp_summaries:
        rsid = f"rs{summary.get('snp_id', summary.get('uid', '?'))}"
        maf = _extract_maf(summary)
        alleles = _extract_alleles(summary)

        # Tính vị trí trong oligo (0-based)
        chrpos = summary.get("chrpos") or summary.get("chrpos_prev_assm")
        if chrpos:
            try:
                genomic_pos = int(str(chrpos).split(":")[1]) if ":" in str(chrpos) else int(chrpos)
                position_in_oligo = genomic_pos - oligo_genomic_start
            except (ValueError, IndexError):
                position_in_oligo = 0
        else:
            position_in_oligo = 0

        # Clamp về [0, oligo_length - 1]
        position_in_oligo = max(0, min(position_in_oligo, oligo_length - 1))

        is_3prime = position_in_oligo >= (oligo_length - 3)

        snp_info = SNPInfo(
            rsid=rsid,
            position_in_oligo=position_in_oligo,
            maf=maf,
            alleles=alleles,
            is_3prime=is_3prime,
        )
        snp_list.append(snp_info)

        # Xác định critical SNP
        if is_probe:
            # Với probe: mọi SNP đều critical
            critical_snps.append(snp_info)
        elif is_3prime or maf >= maf_threshold:
            critical_snps.append(snp_info)

    # Xếp loại
    has_snp = len(snp_list) > 0
    if not has_snp:
        recommendation = "PASS"
    elif critical_snps:
        recommendation = "FAIL"
    else:
        recommendation = "WARNING"

    return SNPResult(
        oligo_id=oligo_id,
        has_snp=has_snp,
        snp_list=snp_list,
        critical_snps=critical_snps,
        recommendation=recommendation,
    )


class SNPChecker:
    """Kiểm tra SNP trong vùng primer/probe từ dbSNP.

    Args:
        config: AppConfig đã load.
        maf_threshold: Ngưỡng MAF để flag SNP (default: 0.01).

    Example::

        checker = SNPChecker(config)
        pairs_flagged = checker.check_all(primer_pairs, seq_record)
    """

    def __init__(
        self,
        config: AppConfig,
        maf_threshold: float = _DEFAULT_MAF_THRESHOLD,
    ) -> None:
        self.config = config
        self.maf_threshold = maf_threshold

    def check_all(
        self,
        primer_pairs: list[PrimerPair],
        seq_record: object,
    ) -> list[PrimerPair]:
        """Kiểm tra SNP cho tất cả primer pairs.

        Gán flags dạng string vào ``pair.snp_flags`` (ví dụ:
        ``"LEFT:WARNING:rs123@pos5"``).  Pair bị lỗi khi query sẽ được log
        warning và bỏ qua (snp_flags giữ nguyên).

        Args:
            primer_pairs: Danh sách PrimerPair cần kiểm tra.
            seq_record: BioPython SeqRecord (dùng ``.id`` làm accession).
                        Nếu ``None`` hoặc không có ``.id``, bỏ qua SNP check.

        Returns:
            Danh sách PrimerPair đã được gán ``snp_flags`` (nếu thành công).
        """
        accession: str = ""
        if seq_record is not None and hasattr(seq_record, "id"):
            accession = str(seq_record.id).split(".")[0]  # bỏ version suffix

        if not accession:
            logger.warning("SNPChecker: không có accession, bỏ qua SNP check")
            return primer_pairs

        for pair in primer_pairs:
            try:
                self._check_pair_snps(pair, accession)
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "SNPChecker: lỗi khi kiểm tra pair '%s' — %s",
                    pair.pair_id,
                    exc,
                )
        return primer_pairs

    def _check_pair_snps(self, pair: PrimerPair, accession: str) -> None:
        """Kiểm tra SNP cho một PrimerPair và gán flags.

        Args:
            pair: PrimerPair cần kiểm tra.
            accession: GenBank/RefSeq accession.
        """
        left = pair.left_primer
        right = pair.right_primer

        left_result = self.check_oligo(
            oligo_sequence=left.sequence,
            genomic_start=left.start + 1,  # chuyển 0-based → 1-based
            genomic_end=left.start + left.length,
            accession=accession,
            is_probe=False,
        )
        right_result = self.check_oligo(
            oligo_sequence=right.sequence,
            genomic_start=right.start + 1,
            genomic_end=right.start + right.length,
            accession=accession,
            is_probe=False,
        )

        for result, label in [(left_result, "LEFT"), (right_result, "RIGHT")]:
            if result.has_snp:
                for snp in result.snp_list:
                    flag = (
                        f"{label}:{result.recommendation}:{snp.rsid}"
                        f"@pos{snp.position_in_oligo}"
                        f"(MAF={snp.maf:.3f})"
                    )
                    pair.snp_flags.append(flag)

        if pair.probe is not None:
            probe = pair.probe
            probe_result = self.check_oligo(
                oligo_sequence=probe.sequence,
                genomic_start=probe.start + 1,
                genomic_end=probe.start + probe.length,
                accession=accession,
                is_probe=True,
            )
            if probe_result.has_snp:
                for snp in probe_result.snp_list:
                    flag = (
                        f"PROBE:{probe_result.recommendation}:{snp.rsid}"
                        f"@pos{snp.position_in_oligo}"
                        f"(MAF={snp.maf:.3f})"
                    )
                    pair.snp_flags.append(flag)

    def check_oligo(
        self,
        oligo_sequence: str,
        genomic_start: int,
        genomic_end: int,
        accession: str = "",
        is_probe: bool = False,
    ) -> SNPResult:
        """Kiểm tra SNP cho một oligo qua NCBI dbSNP API.

        Query NCBI eSearch để tìm SNP IDs trong vùng [genomic_start, genomic_end],
        sau đó lấy chi tiết qua esummary và phân loại kết quả.

        Args:
            oligo_sequence: Chuỗi oligo (5'→3').
            genomic_start: Vị trí genomic bắt đầu (1-based, inclusive).
            genomic_end: Vị trí genomic kết thúc (1-based, inclusive).
            accession: GenBank/RefSeq accession để search dbSNP.
                       Nếu rỗng, trả về SNPResult PASS không query.
            is_probe: True nếu là probe (stricter rules).

        Returns:
            :class:`SNPResult` đã được phân loại.

        Raises:
            requests.RequestException: Nếu NCBI API request thất bại.
        """
        oligo_id = f"{oligo_sequence[:8]}...({genomic_start}-{genomic_end})"

        if not accession:
            logger.debug("SNPChecker: không có accession, trả về PASS cho %s", oligo_id)
            return SNPResult(oligo_id=oligo_id, recommendation="PASS")

        email = self.config.ncbi.email
        api_key = self.config.ncbi.api_key

        snp_ids = _search_snps_by_accession_range(
            accession=accession,
            start=genomic_start,
            end=genomic_end,
            email=email,
            api_key=api_key,
        )

        if not snp_ids:
            return SNPResult(oligo_id=oligo_id, has_snp=False, recommendation="PASS")

        summaries = _fetch_snp_summaries(snp_ids, email=email, api_key=api_key)

        oligo_length = len(oligo_sequence)
        return _build_snp_result(
            oligo_id=oligo_id,
            oligo_length=oligo_length,
            snp_summaries=summaries,
            genomic_start=genomic_start,
            oligo_genomic_start=genomic_start,
            maf_threshold=self.maf_threshold,
            is_probe=is_probe,
        )
