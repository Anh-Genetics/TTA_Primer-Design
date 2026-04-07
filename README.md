# TTA Primer Design

Pipeline thiết kế primer/probe cho PCR, qPCR và TaqMan — sử dụng **Primer3 local** (primer3-py) + **BLAST specificity** + tuỳ chọn **NCBI Primer-BLAST API**.

---

## Mục tiêu

- Thiết kế primer và probe chất lượng cao cho PCR / qPCR / TaqMan
- Kiểm tra tính đặc hiệu (specificity) bằng BLAST
- Phát hiện SNP trong vùng primer/probe (dbSNP)
- Lọc và xếp hạng primer pair theo thuật toán đa tiêu chí
- Xuất báo cáo đa định dạng: CSV, Excel, JSON, FASTA, HTML

---

## Cài đặt

```bash
# Clone repo
git clone https://github.com/Anh-Genetics/TTA_Primer-Design.git
cd TTA_Primer-Design

# Cài đặt package (chế độ phát triển)
pip install -e ".[dev]"
```

### Yêu cầu hệ thống

- Python ≥ 3.11
- pip ≥ 23

---

## Chạy CLI

```bash
# Xem hướng dẫn
tta-primer-design --help
tta-primer-design run --help

# Chạy pipeline với file JSON
tta-primer-design run \
    --config config/pipeline_config.yaml \
    --input data/input/targets.json \
    --output results/run_001/

# Chạy pipeline với file CSV
tta-primer-design run \
    --input data/input/targets.csv \
    --output results/run_002/

# Chạy với log level DEBUG
tta-primer-design run \
    --input data/input/targets.csv \
    --output results/run_003/ \
    --log-level DEBUG
```

### Định dạng input được hỗ trợ

| Định dạng | Mô tả |
|-----------|-------|
| `.json`   | `{"targets": [...]}` — xem `tests/fixtures/sample_input.json` |
| `.csv`    | Cột `id`, `accession`/`sequence`, và các cột tuỳ chọn |
| `.fasta`  | Mỗi FASTA record là một target |
| `.txt`    | Mỗi dòng là một NCBI accession |

---

## Đánh giá cặp mồi sẵn có (`evaluate`)

Lệnh `evaluate` cho phép bạn **kiểm tra chất lượng một cặp primer/probe đã có sẵn**
mà không cần thiết kế lại. Nó thực hiện:

1. **Đánh giá nhiệt động học** — Tm, GC%, ΔG hairpin, ΔG self-dimer, ΔG hetero-dimer,
   GC clamp, repeat runs
2. **Kiểm tra tính đặc hiệu BLAST** (tuỳ chọn, cần internet) — BLAST từng primer
   lên NCBI, phát hiện off-target amplicon tiềm năng, tính specificity score
3. **Điểm tổng hợp** — kết hợp nhiệt động học + specificity

```bash
tta-primer-design evaluate --help
```

### Ví dụ nhanh

```bash
# ① Chỉ đánh giá nhiệt động học (không cần internet, nhanh ~1s)
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --skip-blast

# ② Đánh giá đầy đủ kèm BLAST specificity (~1–2 phút)
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --organism "Homo sapiens" \
    --database nt

# ③ TaqMan probe — có cả probe, lưu kết quả ra JSON
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --probe TGCAGCCACACTTTCTACAATGAGC \
    --name  ACTB_pair1 \
    --output results/eval_ACTB.json

# ④ BLAST bằng database refseq_rna, loài chuột
tta-primer-design evaluate \
    --left  ATGGAGAAAATCTGGCACCAC \
    --right GGGGTGTTGAAGGTCTCAAA \
    --organism "Mus musculus" \
    --database refseq_rna
```

### Giải thích output

```
============================================================
  TTA Primer Design — Evaluate Existing Primer Pair
============================================================
  Pair: ACTB_pair1

─── THERMODYNAMICS ──────────────────────────────────────

  Left  primer: GCAAGGAATGGTTTCAGAAATCCA
    Tm             : 59.8 °C
    GC%            : 41.7%
    Hairpin ΔG     : -1.22 kcal/mol  ✅
    Self-dimer ΔG  : -3.41 kcal/mol  ✅
    3'-end ΔG      : -1.85 kcal/mol
    GC Clamp       : ✅
    Repeat runs    : ✅

  Right primer: CAGGACTCCATGTCGTCCA
    Tm             : 58.5 °C
    GC%            : 57.9%
    ...

  Hetero-dimer ΔG  : -5.23 kcal/mol  ✅

─── BLAST SPECIFICITY ───────────────────────────────────
  Database : nt | Organism : Homo sapiens
  ⏳ Đang chạy BLAST (30–120 giây)…

  Left primer  :  12 hits  | perfect: 1  | 3'-mismatch: 4
  Right primer :   8 hits  | perfect: 1  | 3'-mismatch: 3

  Off-target amplicons : 0  ✅
  Specificity score    : 100.0 / 100

─── SUMMARY ─────────────────────────────────────────────

  Overall Score : 78.3 / 100  →  ✅  PASS
```

### Tùy chọn lệnh

| Tùy chọn | Ngắn | Mô tả | Mặc định |
|----------|------|-------|---------|
| `--left`  | `-l` | Chuỗi primer trái (bắt buộc) | — |
| `--right` | `-r` | Chuỗi primer phải (bắt buộc) | — |
| `--probe` | `-p` | Chuỗi TaqMan probe (tuỳ chọn) | — |
| `--name`  | `-n` | Tên cặp mồi | `pair_001` |
| `--organism` | — | Loài lọc BLAST | `Homo sapiens` |
| `--database` | — | NCBI database: `nt`, `refseq_rna`, `refseq_genomic` | `nt` |
| `--skip-blast` | — | Bỏ qua BLAST (chỉ nhiệt động học) | `False` |
| `--output` | `-o` | Lưu kết quả ra file JSON | — |
| `--config` | `-c` | File YAML config | — |
| `--log-level` | — | Mức log: DEBUG/INFO/WARNING/ERROR | `WARNING` |

### Ngưỡng tham chiếu

| Thông số | Mức tốt | Cảnh báo ⚠️ |
|----------|---------|------------|
| Tm primer | 58–62 °C | < 55 °C hoặc > 65 °C |
| Tm probe (TaqMan) | Tm_primer + 8–10 °C | — |
| GC% | 40–60% | < 30% hoặc > 70% |
| Hairpin ΔG | > −9 kcal/mol | ≤ −9 kcal/mol |
| Self-dimer ΔG | > −9 kcal/mol | ≤ −9 kcal/mol |
| Hetero-dimer ΔG | > −9 kcal/mol | ≤ −9 kcal/mol |
| Off-target amplicons | 0 | ≥ 1 |
| Specificity score | ≥ 80 | < 80 |
| Overall score | ≥ 70 (PASS) | 50–70 (MARGINAL) hoặc < 50 (FAIL) |

---

## Cấu trúc dự án

```
TTA_Primer-Design/
├── src/
│   └── tta_primer_design/          # Package chính (src-layout)
│       ├── __init__.py
│       ├── cli.py                  # CLI (Click)
│       ├── main.py                 # Entry point pipeline
│       ├── config.py               # Load/validate YAML config
│       ├── logging_setup.py        # Logging setup
│       ├── models.py               # Dataclasses: DesignTarget, PrimerPair, ...
│       └── modules/
│           ├── input_parser.py     # Module 1: Parse input
│           ├── sequence_fetcher.py # Module 2: Lấy sequence NCBI
│           ├── sequence_preprocessor.py  # Module 3: Tiền xử lý
│           ├── primer3_runner.py   # Module 4: Chạy Primer3
│           ├── ncbi_primer_blast.py # Module 5: NCBI Primer-BLAST API
│           ├── blast_specificity.py # Module 6: BLAST specificity
│           ├── probe_designer.py   # Module 7: Thiết kế TaqMan probe
│           ├── snp_checker.py      # Module 8: Kiểm tra SNP
│           ├── thermodynamics.py   # Module 9: Nhiệt động học
│           ├── filter_ranker.py    # Module 10: Lọc & xếp hạng
│           └── report_generator.py # Module 11: Báo cáo
├── tests/
│   ├── fixtures/                   # File mẫu cho testing
│   ├── test_config.py
│   ├── test_cli.py
│   ├── test_input_parser.py
│   └── test_models.py
├── config/
│   ├── pipeline_config.yaml        # Cấu hình tổng thể
│   ├── primer3_qpcr_params.yaml    # Tham số qPCR
│   ├── primer3_pcr_params.yaml     # Tham số PCR
│   └── primer3_probe_params.yaml   # Tham số probe
├── data/                           # Input/output data (gitignored)
├── logs/                           # Log files (gitignored)
├── .github/workflows/ci.yml        # GitHub Actions CI
├── pyproject.toml
└── requirements.txt
```

---

## Mô tả các module

| Module | Trạng thái | Mô tả |
|--------|-----------|-------|
| `input_parser.py`        | ✅ Implemented | Parse JSON/CSV/FASTA/TXT → List[DesignTarget] |
| `sequence_fetcher.py`    | ✅ Implemented | Fetch sequence từ NCBI Entrez + cache |
| `sequence_preprocessor.py` | ✅ Implemented | Validate, tính GC, mask low-complexity |
| `primer3_runner.py`      | ✅ Implemented | Chạy primer3-py để thiết kế primer |
| `thermodynamics.py`      | ✅ Implemented | Tính Tm, ΔG, GC clamp, repeat check |
| `probe_designer.py`      | ✅ Implemented | Thiết kế TaqMan probe |
| `filter_ranker.py`       | ✅ Implemented | Lọc và xếp hạng primer pairs (scoring 100 pts) |
| `report_generator.py`    | ✅ Implemented | Xuất báo cáo CSV/Excel/JSON/FASTA/HTML |
| `blast_specificity.py`   | ✅ Implemented | NCBI BLAST API + off-target amplicon detection |
| `ncbi_primer_blast.py`   | 🔲 Sprint 3   | NCBI Primer-BLAST API (submit/poll/parse) |
| `snp_checker.py`         | 🔲 Sprint 3   | Kiểm tra SNP từ dbSNP |

---

## Chạy tests

```bash
# Chạy tất cả tests
pytest

# Chạy với coverage
pytest --cov=tta_primer_design --cov-report=html

# Chạy một test file cụ thể
pytest tests/test_config.py -v
```

---

## Lint & format

```bash
# Kiểm tra code style
ruff check src/ tests/

# Auto-fix (nếu được)
ruff check --fix src/ tests/

# Format code
black src/ tests/
```

---

## Cấu hình pipeline

Xem file `config/pipeline_config.yaml` để tuỳ chỉnh:

```yaml
pipeline:
  mode: "qpcr"              # pcr | qpcr | taqman | sybr | multiplex
  use_local_primer3: true   # true = primer3-py; false = NCBI API

ncbi:
  email: "your@email.com"   # BẮT BUỘC cho NCBI E-utilities
  api_key: ""               # Tùy chọn — tăng rate limit

blast:
  database: "refseq_rna"
  organism: "Homo sapiens"

filters:
  min_specificity_score: 80
  max_off_targets: 0
```

---

## Roadmap

- **Sprint 1** (Tuần 1–2): `sequence_fetcher.py`, `utils/` (NCBI API, rate limiter, validators)
- **Sprint 2** (Tuần 3–4): `sequence_preprocessor.py`, `primer3_runner.py`, `thermodynamics.py`
- **Sprint 3** (Tuần 5–6): `ncbi_primer_blast.py`, `blast_specificity.py`, `probe_designer.py`, `snp_checker.py`
- **Sprint 4** (Tuần 7–8): `filter_ranker.py`, `report_generator.py`, integration tests

---

## License

MIT