"""Data models (dataclasses) cho TTA Primer Design pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# Input models
# ---------------------------------------------------------------------------


@dataclass
class DesignTarget:
    """Một target cần thiết kế primer/probe.

    Attributes:
        target_id: ID duy nhất cho target này.
        input_type: Loại input — "accession" | "sequence" | "gene_name".
        accession: NCBI accession (ví dụ: NM_001101).
        sequence: Chuỗi nucleotide trực tiếp (nếu không có accession).
        gene_name: Tên gene (dùng để tra cứu accession).
        organism: Tên loài NCBI (ví dụ: "Homo sapiens").
        design_mode: Chế độ — "pcr" | "qpcr" | "taqman" | "sybr".
        region_include: Khoảng (start, end) bắt buộc có primer.
        region_exclude: Danh sách khoảng cần tránh.
        exon_junction: True nếu yêu cầu primer span exon-exon junction.
        custom_params: Ghi đè tham số Primer3 cho target này.
    """

    target_id: str
    input_type: str = "accession"
    accession: str | None = None
    sequence: str | None = None
    gene_name: str | None = None
    organism: str = "Homo sapiens"
    design_mode: str = "qpcr"
    region_include: tuple[int, int] | None = None
    region_exclude: list[tuple[int, int]] | None = None
    exon_junction: bool = False
    custom_params: dict | None = None


# ---------------------------------------------------------------------------
# Sequence models
# ---------------------------------------------------------------------------


@dataclass
class ProcessedSequence:
    """Sequence đã qua tiền xử lý, sẵn sàng cho Primer3.

    Attributes:
        sequence: Chuỗi nucleotide (đã mask các vùng lặp/low-complexity).
        excluded_regions: Danh sách vùng loại trừ (SEQUENCE_EXCLUDED_REGION).
        target_regions: Danh sách vùng amplicon bắt buộc (SEQUENCE_TARGET).
        included_region: Vùng bao gồm (SEQUENCE_INCLUDED_REGION).
        exon_junctions: Vị trí các exon-exon junction (0-based).
        gc_content: Tỷ lệ GC (0.0–1.0).
        complexity_score: Điểm độ phức tạp (0.0–1.0, cao = phức tạp hơn).
    """

    sequence: str
    excluded_regions: list[tuple[int, int]] = field(default_factory=list)
    target_regions: list[tuple[int, int]] = field(default_factory=list)
    included_region: tuple[int, int] | None = None
    exon_junctions: list[int] = field(default_factory=list)
    gc_content: float = 0.0
    complexity_score: float = 0.0


# ---------------------------------------------------------------------------
# Oligo / Primer models
# ---------------------------------------------------------------------------


@dataclass
class Oligo:
    """Một oligo đơn (primer hoặc probe).

    Attributes:
        sequence: Chuỗi oligo (5'→3').
        start: Vị trí bắt đầu trên template (0-based).
        length: Độ dài oligo.
        tm: Nhiệt độ nóng chảy (°C).
        gc_percent: Tỷ lệ GC (%).
        self_any_th: ΔG self-dimer (thermodynamic, °C equivalent).
        self_end_th: ΔG 3'-end self-dimer.
        hairpin_th: ΔG hairpin.
        end_stability: Độ bền 3'-end (kcal/mol).
        penalty: Điểm phạt Primer3.
    """

    sequence: str
    start: int = 0
    length: int = 0
    tm: float = 0.0
    gc_percent: float = 0.0
    self_any_th: float = 0.0
    self_end_th: float = 0.0
    hairpin_th: float = 0.0
    end_stability: float = 0.0
    penalty: float = 0.0

    def __post_init__(self) -> None:
        if self.length == 0:
            self.length = len(self.sequence)


@dataclass
class PrimerPair:
    """Một cặp primer (left + right) kèm probe tuỳ chọn.

    Attributes:
        pair_id: ID duy nhất cho cặp primer.
        left_primer: Oligo primer bên trái (forward).
        right_primer: Oligo primer bên phải (reverse).
        probe: Internal oligo (TaqMan probe), None nếu không có.
        amplicon_size: Kích thước amplicon (bp).
        amplicon_sequence: Chuỗi amplicon.
        pair_penalty: Điểm phạt tổng hợp của Primer3.
        specificity_result: Kết quả BLAST specificity (None nếu chưa chạy).
        snp_flags: Danh sách cảnh báo SNP.
        score: Điểm tổng hợp sau khi lọc/xếp hạng (0–100).
    """

    pair_id: str
    left_primer: Oligo
    right_primer: Oligo
    probe: Oligo | None = None
    amplicon_size: int = 0
    amplicon_sequence: str = ""
    pair_penalty: float = 0.0
    specificity_result: object | None = None
    snp_flags: list[str] = field(default_factory=list)
    score: float = 0.0


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------


@dataclass
class DesignResult:
    """Kết quả thiết kế primer cho một DesignTarget.

    Attributes:
        target: Target đã xử lý.
        primer_pairs: Danh sách các PrimerPair kết quả (đã xếp hạng).
        status: "success" | "failed" | "no_primers".
        error: Thông báo lỗi (nếu có).
    """

    target: DesignTarget
    primer_pairs: list[PrimerPair] = field(default_factory=list)
    status: str = "success"
    error: str | None = None
