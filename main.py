"""
Antibody Library MOPSS Analysis
Usage: python main.py --sheet <sample_sheet.csv>
"""

import argparse
import csv
import gzip
import os
from collections import Counter

import regex
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Primer configuration: keyed by (chain, cdr_number)
# chain: "H" or "L"
# cdr_number: 1, 2, or 3 (int)
# ---------------------------------------------------------------------------
PRIMER_CONFIG = {
    ("H", 1): {
        "forward": "GAAGACATAGTTGTCCAGTGCGTCTCC",
        "reverse": "GAAGACAT",
        "expected_length": 72,
        "offset_front": 2,
        "offset_end": 1,
    },
    ("H", 2): {
        "forward": "CGTCTCC",
        "reverse": "GAAGACAT",
        "expected_length": 71,
        "offset_front": 0,
        "offset_end": 2,
    },
    ("H", 3): {
        "forward": "CGTCTCC",
        "reverse": "CGTCTCACGGTGAACGATGAGAAGACAT",
        "expected_length": 79,
        "offset_front": 2,
        "offset_end": 2,
    },
    ("L", 2): {
        "forward": "GAAGACATAGTTGTCCAGTGCGTCTCC",
        "reverse": "GAAGACAT",
        "expected_length": 68,
        "offset_front": 1,
        "offset_end": 1,
    },
    ("L", 3): {
        "forward": "CGTCTCA",
        "reverse": "CGTCTCACGGTGAACGATGAGAAGACAT",
        "expected_length": 71,
        "offset_front": 1,
        "offset_end": 1,
    },
}


def read_sample_sheet(path):
    """Parse sample sheet CSV. Returns list of dicts with keys: name, chain, cdr."""
    samples = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["Sample Name"].strip()
            chain = row["H/L"].strip()
            cdr = row["CDR"].strip()
            if not name or not chain or not cdr:
                continue
            samples.append({"name": name, "chain": chain, "cdr": int(cdr)})
    return samples


def count_fastq_reads(fastq_gz_path):
    """Read a gzipped FASTQ file and return Counter of {sequence: count}."""
    counts = Counter()
    with gzip.open(fastq_gz_path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # '+' line
            f.readline()  # quality line
            if seq:
                counts[seq] += 1
    return counts


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def trim_sequences(seq_counts, primer_for, primer_rev_rc):
    """
    Extract mid-sequences between forward and reverse-complement primers using fuzzy matching.
    Offsets are NOT applied here; they are applied only during AA translation.
    Returns Counter of {mid_seq: count}.
    """
    n_err_for = len(primer_for) // 8 + 1
    n_err_rev = len(primer_rev_rc) // 8 + 1
    pattern = regex.compile(
        f"({primer_for}){{s<={n_err_for}}}(?P<mid_seq>.*)({primer_rev_rc}){{s<={n_err_rev}}}"
    )

    trimmed = Counter()
    for seq, count in seq_counts.items():
        m = pattern.search(seq)
        if m is None:
            continue
        mid = m.group("mid_seq")
        if mid:
            trimmed[mid] += count
    return trimmed


def write_nt_length(trimmed, out_path):
    """Write nucleotide length distribution to CSV."""
    length_counts = Counter()
    for seq, count in trimmed.items():
        length_counts[len(seq)] += count

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["length", "read_count"])
        for length in sorted(length_counts):
            writer.writerow([length, length_counts[length]])


def write_aa_diversity(trimmed, out_path, expected_length, offset_front, offset_end):
    """
    Filter to exact-length sequences, translate to amino acids, write diversity CSV.
    Columns: aa_sequence, read_count, productive
    Returns dict of stats: unique_aa, total_aa_reads, productive_reads, productive_pct
    """
    aa_counts = Counter()
    for seq, count in trimmed.items():
        if len(seq) != expected_length:
            continue
        sub = seq[offset_front:-offset_end] if offset_end else seq[offset_front:]
        aa = str(Seq(sub).translate())
        aa_counts[aa] += count

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aa_sequence", "read_count", "productive"])
        for aa, count in sorted(aa_counts.items(), key=lambda x: -x[1]):
            productive = "*" not in aa
            writer.writerow([aa, count, productive])

    return {
        "unique_aa": len(aa_counts),
        "total_aa_reads": sum(aa_counts.values()),
    }


def main():
    parser = argparse.ArgumentParser(description="Antibody library sequence analysis")
    parser.add_argument("--sheet", required=True, help="Path to sample sheet CSV")
    args = parser.parse_args()

    sheet_dir = os.path.dirname(os.path.abspath(args.sheet))
    data_dir = os.path.join(sheet_dir, "data")
    output_dir = os.path.join(sheet_dir, "output")

    samples = read_sample_sheet(args.sheet)
    print(f"Found {len(samples)} sample(s) in sheet.")

    summary_rows = []

    for sample in samples:
        name = sample["name"]
        chain = sample["chain"]
        cdr = sample["cdr"]
        print(f"\n[{name}] chain={chain}, CDR={cdr}")

        # Locate fastq file
        fastq_path = os.path.join(data_dir, f"{name}.fastq.gz")
        if not os.path.exists(fastq_path):
            print(f"  WARNING: {fastq_path} not found — skipping.")
            continue

        # Look up primer config
        config_key = (chain, cdr)
        if config_key not in PRIMER_CONFIG:
            print(f"  WARNING: No primer config for chain={chain}, CDR={cdr} — skipping.")
            continue
        cfg = PRIMER_CONFIG[config_key]

        primer_for = cfg["forward"]
        primer_rev_rc = reverse_complement(cfg["reverse"])
        expected_length = cfg["expected_length"]
        offset_front = cfg["offset_front"]
        offset_end = cfg["offset_end"]

        # Create output directory
        sample_out_dir = os.path.join(output_dir, name)
        os.makedirs(sample_out_dir, exist_ok=True)

        # Step 1: Count reads
        print(f"  Reading {fastq_path} ...")
        seq_counts = count_fastq_reads(fastq_path)
        total_input_reads = sum(seq_counts.values())
        print(f"  Total unique sequences: {len(seq_counts):,}  |  Total reads: {total_input_reads:,}")

        # Step 2: Trim
        print(f"  Trimming primers ...")
        trimmed = trim_sequences(seq_counts, primer_for, primer_rev_rc)
        total_trimmed_reads = sum(trimmed.values())
        primer_match_pct = round(100 * total_trimmed_reads / total_input_reads, 2) if total_input_reads else 0.0
        print(f"  Trimmed: {len(trimmed):,} unique seqs  |  {total_trimmed_reads:,} reads  "
              f"({primer_match_pct}% matched primers)")

        # Step 3: NT length diversity
        nt_out = os.path.join(sample_out_dir, f"{name}_nt_length.csv")
        write_nt_length(trimmed, nt_out)
        exact_length_reads = sum(c for s, c in trimmed.items() if len(s) == expected_length)
        exact_length_pct = round(100 * exact_length_reads / total_trimmed_reads, 2) if total_trimmed_reads else 0.0
        print(f"  NT length distribution -> {nt_out}")

        # Step 4: AA diversity
        aa_out = os.path.join(sample_out_dir, f"{name}_aa.csv")
        aa_stats = write_aa_diversity(trimmed, aa_out, expected_length, offset_front, offset_end)
        print(f"  AA diversity           -> {aa_out}")

        summary_rows.append({
            "sample_name": name,
            "chain": chain,
            "cdr": cdr,
            "total_reads": total_input_reads,
            "primer_matched_reads": total_trimmed_reads,
            "primer_match_pct": primer_match_pct,
            "expected_length": expected_length,
            "exact_length_reads": exact_length_reads,
            "exact_length_pct": exact_length_pct,
            "unique_aa": aa_stats["unique_aa"],
            "total_aa_reads": aa_stats["total_aa_reads"],
        })

    # Write summary CSV
    if summary_rows:
        os.makedirs(output_dir, exist_ok=True)
        summary_path = os.path.join(output_dir, "summary.csv")
        fieldnames = list(summary_rows[0].keys())
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summary_rows)
        print(f"\nSummary -> {summary_path}")

    print("Done.")


if __name__ == "__main__":
    main()
