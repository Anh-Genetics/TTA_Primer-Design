# TTA Primer Design

> Pipeline thiết kế và đánh giá primer/probe cho PCR, qPCR và TaqMan — kết hợp **Primer3 local** (`primer3-py`), **đánh giá nhiệt động học** (SantaLucia), **BLAST specificity** (NCBI Web API) và **trực quan hóa kết quả BLAST**.

---

## Mục tiêu

- Thiết kế primer và probe chất lượng cao cho **PCR / qPCR / TaqMan / SyBr Green / Multiplex PCR**
- Kiểm tra tính đặc hiệu (specificity) bằng **NCBI BLAST Web API**
- Phát hiện **off-target amplicon** tiềm năng và tính điểm specificity
- Đánh giá **nhiệt động học** đầy đủ: Tm, GC%, ΔG hairpin, self-dimer, hetero-dimer, GC clamp, repeat runs
- Thiết kế **TaqMan probe** tự động với kiểm tra quy tắc thiết kế (strand, Tm, GC%, 5′-G)
- **Lọc và xếp hạng** primer pair theo thuật toán đa tiêu chí (0–100 điểm)
- **Trực quan hóa BLAST** dạng ASCII hoặc file PNG (matplotlib)
- Xuất **báo cáo đa định dạng**: CSV, Excel (.xlsx), JSON, FASTA, HTML

---

## Cài đặt

```bash
# Clone repo
git clone https://github.com/Anh-Genetics/TTA_Primer-Design.git
cd TTA_Primer-Design

# Cài đặt package (chế độ phát triển)
pip install -e ".[dev]"

# Chỉ cài core (không có dev tools)
pip install -e .

# Cài thêm matplotlib để xuất plot PNG
pip install -e ".[viz]"
```

Hoặc cài từ `requirements.txt`:

```bash
pip install -r requirements.txt
```

### Yêu cầu hệ thống

| Yêu cầu | Phiên bản tối thiểu |
|---------|---------------------|
| Python | ≥ 3.11 |
| pip | ≥ 23 |
| Kết nối internet | Cần thiết cho NCBI BLAST/Entrez (tuỳ chọn với `--skip-blast`) |
| matplotlib | Tuỳ chọn — chỉ cần khi dùng `--plot-output` để xuất PNG |

### Dependencies chính

| Package | Mục đích |
|---------|---------|
| `primer3-py` ≥ 2.0 | Thiết kế primer local |
| `biopython` ≥ 1.81 | NCBI Entrez, BLAST, phân tích sequence |
| `pandas` ≥ 2.0 | Xử lý dữ liệu dạng bảng |
| `openpyxl` ≥ 3.1 | Xuất báo cáo Excel (.xlsx) |
| `jinja2` ≥ 3.1 | Render báo cáo HTML |
| `click` ≥ 8.1 | CLI framework |
| `pyyaml` ≥ 6.0 | Đọc file cấu hình YAML |
| `matplotlib` ≥ 3.7 | Xuất plot BLAST dạng PNG (tuỳ chọn) |

---

## Cấu trúc dự án

```
TTA_Primer-Design/
├── src/
│   └── tta_primer_design/              # Package chính (src-layout)
│       ├── __init__.py
│       ├── cli.py                      # CLI entry point (Click)
│       ├── main.py                     # Điều phối pipeline run
│       ├── config.py                   # Load/validate YAML config
│       ├── logging_setup.py            # Cấu hình logging
│       ├── models.py                   # Dataclasses: DesignTarget, PrimerPair, …
│       └── modules/
│           ├── input_parser.py         # Module 1: Parse input JSON/CSV/FASTA/TXT
│           ├── sequence_fetcher.py     # Module 2: Lấy sequence từ NCBI
│           ├── sequence_preprocessor.py # Module 3: Tiền xử lý sequence
│           ├── primer3_runner.py       # Module 4: Chạy Primer3
│           ├── thermodynamics.py       # Module 5: Tính nhiệt động học
│           ├── probe_designer.py       # Module 6: Thiết kế TaqMan probe
│           ├── blast_specificity.py    # Module 7: BLAST specificity check
│           ├── blast_visualizer.py     # Module 8: Trực quan hóa BLAST (ASCII + PNG)
│           ├── filter_ranker.py        # Module 9: Lọc và xếp hạng primer pairs
│           ├── report_generator.py     # Module 10: Xuất báo cáo đa định dạng
│           ├── ncbi_primer_blast.py    # Module 11: NCBI Primer-BLAST API (dự kiến)
│           └── snp_checker.py          # Module 12: Kiểm tra SNP dbSNP (dự kiến)
├── tests/
│   ├── fixtures/                       # File mẫu cho testing
│   ├── test_config.py
│   ├── test_cli.py
│   ├── test_input_parser.py
│   └── test_models.py
├── config/
│   ├── pipeline_config.yaml            # Cấu hình tổng thể pipeline
│   ├── primer3_qpcr_params.yaml        # Tham số Primer3 cho qPCR
│   ├── primer3_pcr_params.yaml         # Tham số Primer3 cho PCR
│   └── primer3_probe_params.yaml       # Tham số Primer3 cho probe
├── data/                               # Input/output data (gitignored)
├── logs/                               # Log files (gitignored)
├── .github/workflows/ci.yml            # GitHub Actions CI
├── pyproject.toml
└── requirements.txt
```

---

## Chạy CLI — Lệnh `run`

Lệnh `run` thực thi toàn bộ pipeline từ đầu đến cuối:  
**Parse input → Fetch sequence NCBI → Tiền xử lý → Primer3 → Thiết kế probe → Lọc/Xếp hạng → Xuất báo cáo**

```bash
# Xem trợ giúp
tta-primer-design --help
tta-primer-design run --help

# Chạy pipeline với file JSON, có file config
tta-primer-design run \
    --config config/pipeline_config.yaml \
    --input data/input/targets.json \
    --output results/run_001/

# Chạy với file CSV (không cần config, dùng giá trị mặc định)
tta-primer-design run \
    --input data/input/targets.csv \
    --output results/run_002/

# Chạy với log level DEBUG để xem chi tiết
tta-primer-design run \
    --input data/input/targets.txt \
    --output results/run_003/ \
    --log-level DEBUG
```

### Tùy chọn lệnh `run`

| Tùy chọn | Ngắn | Bắt buộc | Mô tả | Mặc định |
|----------|------|----------|-------|---------|
| `--config` | `-c` | Không | File cấu hình YAML | — |
| `--input` | `-i` | **Có** | File input (JSON/CSV/FASTA/.txt) | — |
| `--output` | `-o` | **Có** | Thư mục lưu kết quả | — |
| `--log-level` | — | Không | Mức log: `DEBUG`/`INFO`/`WARNING`/`ERROR` | `INFO` |

### Định dạng input được hỗ trợ

| Định dạng | Mô tả | Ví dụ |
|-----------|-------|-------|
| `.json` | Object `{"targets": [...]}`, mỗi target là một object | `tests/fixtures/sample_input.json` |
| `.csv` | Cột `id`, `accession` hoặc `sequence`, và các cột tuỳ chọn | `id,accession,organism,design_mode` |
| `.fasta` | Mỗi FASTA record là một target với sequence trực tiếp | `>ACTB\nATGGAGAAA...` |
| `.txt` | Mỗi dòng là một NCBI accession; dòng `#` và dòng trắng bị bỏ qua | `NM_001101` |

**Ví dụ file `.json`:**

```json
{
  "targets": [
    {
      "id": "ACTB",
      "accession": "NM_001101",
      "design_mode": "qpcr",
      "organism": "Homo sapiens"
    }
  ]
}
```

**Ví dụ file `.csv`:**

```
id,accession,organism,design_mode,exon_junction
ACTB,NM_001101,Homo sapiens,qpcr,true
GAPDH,NM_002046,Homo sapiens,qpcr,false
```

**Ví dụ file `.txt`:**

```
# Danh sách accession cần thiết kế primer
NM_001101
NM_002046
NM_000546
```

### File output của lệnh `run`

Các file được tạo trong thư mục `--output` tùy theo `output.formats` trong config:

| File | Điều kiện | Nội dung |
|------|-----------|---------|
| `results_summary.csv` | `formats` chứa `csv` | Bảng tóm tắt tất cả primer pairs |
| `results.json` | `formats` chứa `json` | Kết quả đầy đủ dạng JSON |
| `results_detailed.xlsx` | `formats` chứa `xlsx` | Báo cáo Excel nhiều sheet |
| `primers.fasta` | `formats` chứa `fasta` | Trình tự primer/probe dạng FASTA |
| `report.html` | `formats` chứa `html` | Báo cáo HTML có thể xem trực tiếp |
| `pipeline.log` | Luôn tạo | File log của lần chạy |

---

## Lệnh `evaluate`

Lệnh `evaluate` cho phép **kiểm tra chất lượng một cặp primer/probe đã có sẵn** mà không cần chạy lại toàn bộ pipeline. Phù hợp để xác minh các primer thiết kế thủ công hoặc lấy từ tài liệu.

```bash
tta-primer-design evaluate --help
```

### Tùy chọn lệnh `evaluate`

| Tùy chọn | Ngắn | Bắt buộc | Mô tả | Mặc định |
|----------|------|----------|-------|---------|
| `--left` | `-l` | **Có** | Chuỗi primer trái (5′→3′) | — |
| `--right` | `-r` | **Có** | Chuỗi primer phải (5′→3′) | — |
| `--probe` | `-p` | Không | Chuỗi TaqMan probe (5′→3′) | — |
| `--name` | `-n` | Không | Tên định danh cặp mồi | `pair_001` |
| `--organism` | — | Không | Loài để lọc kết quả BLAST | `Homo sapiens` |
| `--database` | — | Không | NCBI database: `nt` / `refseq_rna` / `refseq_genomic` | `nt` |
| `--skip-blast` | — | Không | Bỏ qua BLAST, chỉ đánh giá nhiệt động học | `False` |
| `--output` | `-o` | Không | Lưu kết quả ra file JSON | — |
| `--config` | `-c` | Không | File YAML config | — |
| `--log-level` | — | Không | Mức log: `DEBUG`/`INFO`/`WARNING`/`ERROR` | `WARNING` |
| `--plot` | — | Không | In biểu đồ BLAST ASCII ra console | `False` |
| `--plot-output` | — | Không | Lưu biểu đồ BLAST ra file PNG (cần `matplotlib`) | — |

### Ví dụ sử dụng `evaluate`

```bash
# ① Chỉ đánh giá nhiệt động học — không cần internet, hoàn thành ~1 giây
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --skip-blast

# ② Đánh giá đầy đủ kèm BLAST specificity (~1–2 phút, cần internet)
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --organism "Homo sapiens" \
    --database nt

# ③ Có TaqMan probe, lưu kết quả ra file JSON
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --probe TGCAGCCACACTTTCTACAATGAGC \
    --name  ACTB_pair1 \
    --output results/eval_ACTB.json

# ④ BLAST với database refseq_rna, loài chuột
tta-primer-design evaluate \
    --left  ATGGAGAAAATCTGGCACCAC \
    --right GGGGTGTTGAAGGTCTCAAA \
    --organism "Mus musculus" \
    --database refseq_rna

# ⑤ Hiển thị biểu đồ BLAST ASCII ra console (NEW)
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --name  ACTB_pair1 \
    --plot

# ⑥ Lưu biểu đồ BLAST ra file PNG (cần matplotlib) (NEW)
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --name  ACTB_pair1 \
    --plot-output results/blast_ACTB.png
```

### Giải thích output lệnh `evaluate`

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
    Hairpin ΔG     : -0.98 kcal/mol  ✅
    Self-dimer ΔG  : -2.15 kcal/mol  ✅
    3'-end ΔG      : -2.10 kcal/mol
    GC Clamp       : ✅
    Repeat runs    : ✅

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

---

## Trực quan hóa BLAST (`--plot` và `--plot-output`)

Module `blast_visualizer.py` cung cấp hai hình thức trực quan hóa kết quả BLAST:

### `--plot` — Biểu đồ ASCII ra console

Khi truyền flag `--plot`, sau khi BLAST hoàn tất, một bản đồ ASCII sẽ được in ra console, cho thấy số lượng hit của mỗi primer và phân phối số mismatch:

```
═══ BLAST VISUALIZATION ════════════════════════════════════
  Pair: ACTB_pair1

  Left primer hits:  12  [████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░]
  Right primer hits:  8  [████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]

  Mismatch distribution:
  ┌─────────────┬────────────────┬─────────────────┐
  │ Mismatches  │ Left primer    │ Right primer    │
  ├─────────────┼────────────────┼─────────────────┤
  │ Perfect (0) │ 1              │ 1               │
  │ 1 mismatch  │ 4              │ 3               │
  │ 2 mismatches│ 5              │ 2               │
  │ 3+ mismatches│ 2             │ 2               │
  └─────────────┴────────────────┴─────────────────┘
```

### `--plot-output PATH` — Lưu file PNG

Khi truyền `--plot-output <path>`, biểu đồ matplotlib sẽ được lưu ra file PNG. **Yêu cầu `matplotlib` được cài đặt** (xem phần Cài đặt).

```bash
tta-primer-design evaluate \
    --left  GCAAGGAATGGTTTCAGAAATCCA \
    --right CAGGACTCCATGTCGTCCA \
    --name  ACTB_pair1 \
    --plot-output results/blast_ACTB.png
```

> **Lưu ý:** `--plot` và `--plot-output` có thể dùng đồng thời. `--plot` không yêu cầu `matplotlib` (chỉ dùng ký tự Unicode thuần).

---

## Cấu hình pipeline

Pipeline được cấu hình qua file YAML, truyền qua `--config`. Nếu không truyền config, các giá trị mặc định được dùng.

**Ví dụ đầy đủ `config/pipeline_config.yaml`:**

```yaml
pipeline:
  name: "Primer_Probe_Design_Pipeline"
  version: "0.1.0"
  mode: "qpcr"          # pcr | qpcr | taqman | sybr | multiplex
  use_local_primer3: true
  top_n_pairs: 5        # Số lượng primer pair tốt nhất giữ lại

ncbi:
  email: "your@email.com"   # BẮT BUỘC cho NCBI E-utilities
  api_key: ""               # Tùy chọn — tăng rate limit lên 10 req/s
  retries: 3
  timeout: 120
  rate_limit: 3             # Số request/giây (không có API key)

blast:
  database: "refseq_rna"    # nt | refseq_rna | refseq_genomic
  organism: "Homo sapiens"
  max_targets: 500
  word_size: 7
  evalue: 1000.0

output:
  formats: ["csv", "json"]  # csv | json | xlsx | fasta | html
  include_plots: false

filters:
  min_specificity_score: 80
  max_off_targets: 0
  require_exon_junction: false
  avoid_snp_in_primer: true
  avoid_snp_in_probe: true
```

### Mô tả các trường cấu hình

| Trường | Kiểu | Mô tả |
|--------|------|-------|
| `pipeline.mode` | string | Chế độ thiết kế: `pcr`, `qpcr`, `taqman`, `sybr`, `multiplex` |
| `pipeline.top_n_pairs` | int | Giữ lại N cặp mồi tốt nhất cho mỗi target |
| `pipeline.use_local_primer3` | bool | `true` = dùng primer3-py local |
| `ncbi.email` | string | **Bắt buộc** — địa chỉ email đăng ký với NCBI |
| `ncbi.api_key` | string | API key NCBI (tùy chọn, tăng giới hạn rate) |
| `blast.database` | string | Database NCBI để BLAST |
| `blast.evalue` | float | E-value threshold cho BLAST |
| `output.formats` | list | Danh sách định dạng xuất báo cáo |
| `filters.min_specificity_score` | int | Điểm specificity tối thiểu để giữ primer pair |
| `filters.max_off_targets` | int | Số off-target amplicon tối đa cho phép |

### Hành vi theo `pipeline.mode`

| Mode | Thiết kế probe | Mô tả |
|------|---------------|-------|
| `pcr` | Không | PCR thông thường |
| `sybr` | Không | qPCR với SYBR Green (không cần probe) |
| `multiplex` | Không | Multiplex PCR |
| `qpcr` | **Có** (TaqMan) | qPCR định lượng với TaqMan probe |
| `taqman` | **Có** (TaqMan) | Tương đương `qpcr` |

---

## Mô tả các module

| Module | Trạng thái | Mô tả |
|--------|-----------|-------|
| `input_parser.py` | ✅ Hoàn thành | Parse JSON/CSV/FASTA/TXT → `List[DesignTarget]` |
| `sequence_fetcher.py` | ✅ Hoàn thành | Fetch FASTA + GenBank từ NCBI Entrez, retry tự động, in-memory cache |
| `sequence_preprocessor.py` | ✅ Hoàn thành | Validate, uppercase, tính GC%, đánh giá độ phức tạp, mask low-complexity |
| `primer3_runner.py` | ✅ Hoàn thành | Wrapper `primer3-py`, parse kết quả thành `List[PrimerPair]` |
| `thermodynamics.py` | ✅ Hoàn thành | Tính Tm (SantaLucia), ΔG hairpin, self-dimer, hetero-dimer, GC clamp, repeat check |
| `probe_designer.py` | ✅ Hoàn thành | Thiết kế TaqMan probe: Tm = Tm_primer_avg + 9 °C, chọn strand, kiểm tra quy tắc |
| `blast_specificity.py` | ✅ Hoàn thành | NCBI BLAST Web API (biopython), phát hiện off-target amplicon, tính specificity score |
| `blast_visualizer.py` | ✅ Hoàn thành | Biểu đồ ASCII + xuất PNG matplotlib cho kết quả BLAST |
| `filter_ranker.py` | ✅ Hoàn thành | Hard filter (off-target, SNP FAIL, Tm range), soft flag, tính điểm, xếp hạng |
| `report_generator.py` | ✅ Hoàn thành | Xuất báo cáo CSV / Excel / JSON / FASTA / HTML |
| `ncbi_primer_blast.py` | 🔲 Dự kiến | NCBI Primer-BLAST API (submit / poll / parse) |
| `snp_checker.py` | 🔲 Dự kiến | Kiểm tra SNP từ NCBI dbSNP |

---

## Điểm scoring và ngưỡng tham chiếu

### Công thức tính điểm tổng hợp (0–100)

| Thành phần | Điểm tối đa | Cơ sở tính |
|------------|------------|------------|
| Specificity | 30 | Dựa trên kết quả BLAST (số off-target hit, off-target amplicon) |
| Nhiệt động học | 25 | Tm gần 60 °C, ΔG hairpin không quá thấp |
| SNP-free | 20 | Dựa trên `snp_flags` (không có SNP trong vùng primer/probe) |
| Primer3 penalty (nghịch đảo) | 15 | Điểm penalty Primer3 càng thấp càng tốt |
| Kích thước amplicon | 10 | Tối ưu 70–200 bp, tốt nhất ~125 bp |

### Ngưỡng đánh giá kết quả

| Thông số | Mức tốt ✅ | Cảnh báo ⚠️ | Loại ❌ |
|----------|-----------|------------|--------|
| Tm primer | 58–62 °C | 55–58 °C hoặc 62–65 °C | < 55 °C hoặc > 65 °C |
| Tm probe (TaqMan) | Tm_primer + 8–10 °C | — | — |
| GC% primer | 40–60% | 30–40% hoặc 60–70% | < 30% hoặc > 70% |
| Hairpin ΔG | > −9 kcal/mol | −9 đến −12 kcal/mol | ≤ −12 kcal/mol |
| Self-dimer ΔG | > −9 kcal/mol | −9 đến −12 kcal/mol | ≤ −12 kcal/mol |
| Hetero-dimer ΔG | > −9 kcal/mol | −9 đến −12 kcal/mol | ≤ −12 kcal/mol |
| Off-target amplicons | 0 | — | ≥ 1 |
| Specificity score | ≥ 80 | 60–80 | < 60 |
| Overall score | ≥ 70 → **PASS** | 50–70 → **MARGINAL** | < 50 → **FAIL** |

---

## Chạy tests

```bash
# Chạy tất cả tests
pytest

# Chạy với coverage report
pytest --cov=tta_primer_design --cov-report=html

# Chạy một test file cụ thể
pytest tests/test_config.py -v

# Chạy một test function cụ thể
pytest tests/test_input_parser.py::test_parse_json -v
```

---

## Lint & format

```bash
# Kiểm tra code style với ruff
ruff check src/ tests/

# Auto-fix các lỗi có thể tự sửa
ruff check --fix src/ tests/

# Format code với black
black src/ tests/

# Kiểm tra type hints với mypy
mypy src/
```

---

## Roadmap / Trạng thái tính năng

| Tính năng | Trạng thái | Ghi chú |
|-----------|-----------|---------|
| Parse input JSON/CSV/FASTA/TXT | ✅ Hoàn thành | — |
| Fetch sequence từ NCBI Entrez | ✅ Hoàn thành | Retry + in-memory cache |
| Tiền xử lý sequence | ✅ Hoàn thành | Validate, GC%, mask low-complexity |
| Thiết kế primer với Primer3 | ✅ Hoàn thành | Dùng `primer3-py` local |
| Đánh giá nhiệt động học | ✅ Hoàn thành | Tm, ΔG hairpin/dimer, GC clamp |
| Thiết kế TaqMan probe | ✅ Hoàn thành | Tm target = primer avg + 9 °C |
| BLAST specificity check | ✅ Hoàn thành | NCBI Web API qua biopython |
| Trực quan hóa BLAST ASCII | ✅ Hoàn thành | `--plot` flag |
| Trực quan hóa BLAST PNG | ✅ Hoàn thành | `--plot-output` + matplotlib |
| Lọc & xếp hạng primer pairs | ✅ Hoàn thành | Hard filter + soft flag + scoring |
| Báo cáo CSV/JSON/XLSX/FASTA/HTML | ✅ Hoàn thành | — |
| NCBI Primer-BLAST API | 🔲 Dự kiến | `ncbi_primer_blast.py` |
| Kiểm tra SNP từ dbSNP | 🔲 Dự kiến | `snp_checker.py` |

---

## License

MIT — xem file [LICENSE](LICENSE) để biết thêm chi tiết.
