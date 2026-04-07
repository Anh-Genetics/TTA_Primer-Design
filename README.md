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
| `sequence_fetcher.py`    | 🔲 Sprint 1   | Fetch sequence từ NCBI Entrez + cache |
| `sequence_preprocessor.py` | 🔲 Sprint 2 | Mask repeats, validate, tính GC |
| `primer3_runner.py`      | 🔲 Sprint 2   | Chạy primer3-py để thiết kế primer |
| `ncbi_primer_blast.py`   | 🔲 Sprint 3   | NCBI Primer-BLAST API (submit/poll/parse) |
| `blast_specificity.py`   | 🔲 Sprint 3   | BLAST specificity check |
| `probe_designer.py`      | 🔲 Sprint 3   | Thiết kế TaqMan probe |
| `snp_checker.py`         | 🔲 Sprint 3   | Kiểm tra SNP từ dbSNP |
| `thermodynamics.py`      | 🔲 Sprint 2   | Tính Tm, ΔG, GC clamp |
| `filter_ranker.py`       | 🔲 Sprint 4   | Lọc và xếp hạng primer pairs |
| `report_generator.py`    | 🔲 Sprint 4   | Xuất báo cáo CSV/Excel/JSON/HTML |

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