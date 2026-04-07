"""Module 8 — SNPChecker: kiểm tra SNP trong vùng primer/probe.

Nguồn dữ liệu:
    A) NCBI Variation Services API (online)
       Endpoint: https://api.ncbi.nlm.nih.gov/variation/v0/
    B) dbSNP VCF file local (khuyến nghị cho batch)
       Dùng pysam hoặc cyvcf2 để parse

Logic:
    1. Lấy genomic coordinates của primer/probe
    2. Query dbSNP cho region đó
    3. Flag SNPs trong primer với MAF > threshold (default: 0.01)
    4. Cảnh báo đặc biệt SNP ở 3' end primer
    5. Cảnh báo tất cả SNP trong probe

Tiêu chí đánh giá:
    - "PASS"    : không có SNP
    - "WARNING" : có SNP nhưng MAF < threshold hoặc ở giữa primer
    - "FAIL"    : SNP ở 3' end hoặc MAF > threshold trong probe
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import requests

from tta_primer_design.config import AppConfig
from tta_primer_design.models import PrimerPair

logger = logging.getLogger("tta_primer_design.modules.snp_checker")

_DEFAULT_MAF_THRESHOLD = 0.01

# NCBI eUtils base URL
_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI Variation Services API
_VARIATION_URL = "https://api.ncbi.nlm.nih.gov/variation/v0"

# 3' prime region length (bases from 3' end considered critical)
_THREE_PRIME_LEN = 3


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
        self._session = requests.Session()
        self._current_accession: str | None = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def check_all(
        self,
        primer_pairs: list[PrimerPair],
        seq_record: object,
    ) -> list[PrimerPair]:
        """Kiểm tra SNP cho tất cả primer pairs.

        Gán flags vào ``pair.snp_flags``.

        Args:
            primer_pairs: Danh sách PrimerPair.
            seq_record: BioPython SeqRecord (có chứa accession/coords).

        Returns:
            Danh sách PrimerPair đã được gán snp_flags.
        """
        self._current_accession = getattr(seq_record, "id", None)

        for pair in primer_pairs:
            flags: list[str] = []

            # Check left primer
            left = pair.left_primer
            left_end = left.start + len(left.sequence)
            left_result = self.check_oligo(left.sequence, left.start, left_end, is_probe=False)
            if left_result.recommendation != "PASS":
                flags.append(f"LEFT:{left_result.recommendation}")

            # Check right primer
            right = pair.right_primer
            right_end = right.start + len(right.sequence)
            right_result = self.check_oligo(right.sequence, right.start, right_end, is_probe=False)
            if right_result.recommendation != "PASS":
                flags.append(f"RIGHT:{right_result.recommendation}")

            # Check probe if present
            if pair.probe is not None:
                probe = pair.probe
                probe_end = probe.start + len(probe.sequence)
                probe_result = self.check_oligo(
                    probe.sequence, probe.start, probe_end, is_probe=True
                )
                if probe_result.recommendation != "PASS":
                    flags.append(f"PROBE:{probe_result.recommendation}")

            pair.snp_flags = flags

        return primer_pairs

    def check_oligo(
        self,
        oligo_sequence: str,
        genomic_start: int,
        genomic_end: int,
        is_probe: bool = False,
    ) -> SNPResult:
        """Kiểm tra SNP cho một oligo qua NCBI dbSNP.

        Sử dụng NCBI eUtils esearch để tìm SNP trong vùng genomic,
        sau đó efetch để lấy thông tin MAF và allele.

        Args:
            oligo_sequence: Chuỗi oligo.
            genomic_start: Vị trí genomic bắt đầu (0-based).
            genomic_end: Vị trí genomic kết thúc.
            is_probe: True nếu là probe (stricter rules).

        Returns:
            SNPResult.
        """
        oligo_id = f"{genomic_start}_{genomic_end}"
        oligo_len = len(oligo_sequence)
        result = SNPResult(oligo_id=oligo_id)

        rsids = self._search_snps_in_region(genomic_start, genomic_end)
        if not rsids:
            return result

        snp_infos = self._fetch_snp_details(rsids, genomic_start, oligo_len, is_probe)
        if not snp_infos:
            return result

        result.has_snp = True
        result.snp_list = snp_infos
        result.critical_snps = [s for s in snp_infos if s.is_3prime or is_probe]
        result.recommendation = self._determine_recommendation(snp_infos, is_probe)

        return result

    # ------------------------------------------------------------------
    # Private helpers — NCBI API
    # ------------------------------------------------------------------

    def _search_snps_in_region(self, genomic_start: int, genomic_end: int) -> list[str]:
        """Tìm rsIDs trong vùng genomic qua NCBI eUtils esearch.

        Args:
            genomic_start: Vị trí bắt đầu (0-based).
            genomic_end: Vị trí kết thúc.

        Returns:
            Danh sách rsID (strings, không có prefix "rs").
        """
        accession = self._current_accession
        if not accession:
            logger.debug("No accession available; skipping SNP search")
            return []

        # NCBI SNP esearch: search by accession + position range
        # CHRPOS uses 1-based coordinates
        term = f"{accession}[accn] AND {genomic_start + 1}:{genomic_end}[CHRPOS]"

        params: dict[str, str | int] = {
            "db": "snp",
            "term": term,
            "retmode": "json",
            "retmax": 100,
            "usehistory": "n",
        }
        if self.config.ncbi.email:
            params["email"] = self.config.ncbi.email
        if self.config.ncbi.api_key:
            params["api_key"] = self.config.ncbi.api_key

        try:
            resp = self._session.get(
                f"{_EUTILS_URL}/esearch.fcgi",
                params=params,
                timeout=self.config.ncbi.timeout,
            )
            resp.raise_for_status()
            data = resp.json()
            ids: list[str] = data.get("esearchresult", {}).get("idlist", [])
            return ids
        except Exception as exc:
            logger.warning(
                "SNP esearch failed for region %d-%d: %s", genomic_start, genomic_end, exc
            )
            return []

    def _fetch_snp_details(
        self,
        rsids: list[str],
        region_start: int,
        oligo_len: int,
        is_probe: bool,
    ) -> list[SNPInfo]:
        """Lấy thông tin chi tiết SNP từ NCBI Variation Services API.

        Args:
            rsids: Danh sách rsIDs (numeric strings).
            region_start: Vị trí bắt đầu của oligo (genomic, 0-based).
            oligo_len: Độ dài oligo (để tính vị trí SNP trong oligo).
            is_probe: True nếu là probe.

        Returns:
            Danh sách SNPInfo.
        """
        snp_list: list[SNPInfo] = []

        for rsid_str in rsids[:50]:  # limit to avoid rate-limiting
            try:
                resp = self._session.get(
                    f"{_VARIATION_URL}/beta/refsnp/{rsid_str}",
                    timeout=self.config.ncbi.timeout,
                )
                if resp.status_code == 404:
                    continue
                resp.raise_for_status()
                data = resp.json()
                snp_info = self._parse_refsnp_response(
                    data, rsid_str, region_start, oligo_len, is_probe
                )
                if snp_info:
                    snp_list.append(snp_info)
            except Exception as exc:
                logger.debug("Failed to fetch rs%s: %s", rsid_str, exc)

        return snp_list

    def _parse_refsnp_response(
        self,
        data: dict,
        rsid_str: str,
        region_start: int,
        oligo_len: int,
        is_probe: bool,
    ) -> SNPInfo | None:
        """Parse RefSNP JSON response thành SNPInfo.

        Args:
            data: JSON dict từ NCBI Variation Services.
            rsid_str: rsID numeric string.
            region_start: Vị trí bắt đầu của oligo (genomic, 0-based).
            oligo_len: Độ dài oligo.
            is_probe: True nếu là probe.

        Returns:
            SNPInfo hoặc None nếu parse thất bại.
        """
        try:
            primary_snapshot = data.get("primary_snapshot_data", {})

            # Extract MAF from allele annotations
            maf = 0.0
            alleles_str = ""
            for ann_item in primary_snapshot.get("allele_annotations", []):
                for freq in ann_item.get("frequency", []):
                    study_maf = freq.get("minor_allele_freq")
                    if study_maf is not None:
                        maf = max(maf, float(study_maf))

            # Extract alleles from placements
            allele_list: list[str] = []
            placements = primary_snapshot.get("placements_with_allele", [])
            if placements:
                for al in placements[0].get("alleles", []):
                    spdi = al.get("allele", {}).get("spdi", {})
                    alt = spdi.get("inserted_sequence", "")
                    if alt:
                        allele_list.append(alt)
            alleles_str = "/".join(allele_list)

            # Estimate position in oligo: use the middle of the SNP region
            snp_genomic_pos = _extract_genomic_pos(data, placements)
            pos_in_oligo = snp_genomic_pos - region_start if snp_genomic_pos >= region_start else 0
            pos_in_oligo = max(0, min(pos_in_oligo, oligo_len - 1))

            # 3' end check: last _THREE_PRIME_LEN bases of oligo
            is_3prime = not is_probe and (oligo_len - pos_in_oligo <= _THREE_PRIME_LEN)

            return SNPInfo(
                rsid=f"rs{rsid_str}",
                position_in_oligo=pos_in_oligo,
                maf=maf,
                alleles=alleles_str,
                is_3prime=is_3prime,
            )
        except Exception as exc:
            logger.debug("Failed to parse RefSNP rs%s: %s", rsid_str, exc)
            return None

    # ------------------------------------------------------------------
    # Private helpers — classification
    # ------------------------------------------------------------------

    def _determine_recommendation(self, snp_list: list[SNPInfo], is_probe: bool) -> str:
        """Xác định mức độ cảnh báo SNP.

        Rules:
            - FAIL  : SNP ở 3' end primer, hoặc MAF > threshold trong probe
            - WARNING: SNP có MAF < threshold và không ở 3' end
            - PASS  : không có SNP

        Args:
            snp_list: Danh sách SNPInfo.
            is_probe: True nếu là probe.

        Returns:
            "PASS" | "WARNING" | "FAIL"
        """
        if not snp_list:
            return "PASS"

        for snp in snp_list:
            # Any SNP in a probe is at least WARNING; high MAF = FAIL
            if is_probe:
                if snp.maf >= self.maf_threshold:
                    return "FAIL"
            else:
                # 3' end SNP = FAIL regardless of MAF
                if snp.is_3prime:
                    return "FAIL"
                # High MAF SNP in primer = FAIL
                if snp.maf >= self.maf_threshold:
                    return "FAIL"

        return "WARNING"


def _extract_genomic_pos(data: dict, placements: list[dict]) -> int:
    """Trích xuất vị trí genomic (0-based) từ RefSNP JSON.

    Args:
        data: Full RefSNP JSON dict.
        placements: placements_with_allele list.

    Returns:
        Genomic position (0-based) hoặc 0 nếu không tìm thấy.
    """
    try:
        if placements:
            alleles = placements[0].get("alleles", [])
            if alleles:
                spdi = alleles[0].get("allele", {}).get("spdi", {})
                pos = spdi.get("position")
                if pos is not None:
                    return int(pos)
    except (KeyError, TypeError, ValueError):
        pass
    return 0
