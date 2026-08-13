#!/usr/bin/env python3
"""
Generate test data for inSVert BND processing tests.

Reference sequences use deterministic repeating patterns:
  chr1: "ACGT" * 2500  (10,000 bp)  → pos p (1-indexed): "ACGT"[(p-1) % 4]
  chr2: "TGCA" * 2500  (10,000 bp)  → pos p (1-indexed): "TGCA"[(p-1) % 4]
"""
import os
import sys

DATA_DIR = os.path.dirname(os.path.abspath(__file__))


def generate_test_reference():
    """Create a small test reference FASTA with chr1 and chr2."""
    fasta_path = os.path.join(DATA_DIR, "test_ref.fa")

    chr1_seq = "ACGT" * 2500  # 10000bp
    chr2_seq = "TGCA" * 2500  # 10000bp

    with open(fasta_path, 'w') as f:
        f.write(">chr1\n")
        for i in range(0, len(chr1_seq), 60):
            f.write(chr1_seq[i:i+60] + "\n")
        f.write(">chr2\n")
        for i in range(0, len(chr2_seq), 60):
            f.write(chr2_seq[i:i+60] + "\n")

    print(f"Generated: {fasta_path}")
    return fasta_path


if __name__ == "__main__":
    fasta_path = generate_test_reference()
    try:
        import pysam
        pysam.faidx(fasta_path)
        print(f"Indexed:   {fasta_path}.fai")
    except ImportError:
        print("WARNING: pysam not available, .fai index not created")
