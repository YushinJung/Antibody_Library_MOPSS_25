# Antibody Library MOPSS Analysis

A lightweight command-line tool for analyzing antibody library sequencing data.
Given a sample sheet and raw FASTQ files, it trims primer sequences and outputs:

- Nucleotide length distribution per sample
- Amino acid diversity table per sample
- Run summary across all samples

---

## System Requirements

**Tested on:**
- Python 3.9, 3.10, 3.11 on Linux (Ubuntu 20.04/22.04) and macOS (12+)
- No GPU required; standard CPU only

**Required non-standard hardware:**
- None. Any standard desktop or laptop computer is sufficient.

**Typical install time:**
- ~1–2 minutes (pip install on a normal desktop with internet access)

**Expected run time for demo:**
- ~10–30 seconds for the provided demo file (`data/VH_CDR1.fastq.gz`, ~50k reads) on a standard desktop computer

---

## Package Requirements

- Python 3.8+
- [`regex`](https://pypi.org/project/regex/)
- [`biopython`](https://biopython.org/)

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/YushinJung/Antibody_Library_MOPSS_25.git
cd Antibody_Library_MOPSS_25

# 2. Install dependencies
pip install regex biopython
```

---

## Quick Use Guide

### 1. Prepare your data

Place your raw FASTQ files (`.fastq.gz`) in the `data/` folder.
File names must match the `Sample Name` column in the sample sheet.

```
data/
  VH_CDR1.fastq.gz
  VH_CDR2.fastq.gz
  ...
```

### 2. Prepare the sample sheet

Create a CSV file with three columns:

| Sample Name | H/L | CDR |
|-------------|-----|-----|
| VH_CDR1     | H   | 1   |
| VH_CDR2     | H   | 2   |

- **Sample Name**: must match the FASTQ filename (without `.fastq.gz`)
- **H/L**: `H` for heavy chain, `L` for light chain
- **CDR**: `1`, `2`, or `3`

A demo sheet is provided: [`ABMOPSS_sample_sheet.csv`](ABMOPSS_sample_sheet.csv)

### 3. Run the analysis

```bash
python main.py --sheet ABMOPSS_sample_sheet.csv
```

### 4. Check the output

Results are written to `output/<Sample Name>/`:

```
output/
  summary.csv               # run summary across all samples
  VH_CDR1/
    VH_CDR1_nt_length.csv   # read count per nucleotide length
    VH_CDR1_aa.csv          # amino acid sequences and read counts
```

**`_nt_length.csv`**

| length | read_count |
|--------|------------|
| 70     | 120        |
| 72     | 48302      |
| 74     | 88         |

**`_aa.csv`**

| aa_sequence | read_count | productive |
|-------------|------------|------------|
| ARSGYYDS... | 12045      | True       |
| ARSGYYDS*.. | 32         | False      |

- `productive`: `True` if no stop codon (`*`) is present in the translated sequence.

**`summary.csv`**

| sample_name | chain | cdr | total_reads | primer_matched_reads | primer_match_pct | expected_length | exact_length_reads | exact_length_pct | unique_aa | total_aa_reads |
|-------------|-------|-----|-------------|----------------------|------------------|-----------------|--------------------|------------------|-----------|----------------|
| VH_CDR1     | H     | 1   | 50000       | 46200                | 92.4             | 72              | 44800              | 97.0             | 3210      | 44800          |

- `primer_match_pct`: percentage of reads where both primers were found
- `exact_length_pct`: percentage of primer-matched reads at the expected CDR length
- `unique_aa`: number of unique amino acid sequences
- `total_aa_reads`: reads used for AA analysis (exact-length only)

---

## Supported Chain / CDR Combinations

| Chain | CDR |
|-------|-----|
| H     | 1   |
| H     | 2   |
| H     | 3   |
| L     | 2   |
| L     | 3   |

---

## How It Works

1. Raw FASTQ reads are counted by unique sequence.
2. Primer sequences flanking each CDR are located using fuzzy matching (up to ~1 mismatch per 8 bp).
3. The CDR region between primers is extracted and trimmed by per-CDR offsets.
4. Nucleotide length distribution is written to `_nt_length.csv`.
5. Sequences matching the expected CDR length are translated (BioPython standard codon table) and written to `_aa.csv`.
6. A run-level `summary.csv` is written to `output/` consolidating key statistics across all samples.
