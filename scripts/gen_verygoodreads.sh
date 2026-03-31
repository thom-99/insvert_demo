#!/bin/bash 


# === quick check ===
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reference.fasta> <output_reads.fastq>"
    exit 1
fi

# === inputs ===
REFERENCE_FA=$(realpath "$1")
OUTPUT=$(realpath "$2")

badread simulate --reference $REFERENCE_FA --quantity 50x --error_model random \
    --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
    --identity 30,3 --length 40000,20000 --start_adapter_seq "" --end_adapter_seq "" \
    | gzip > $OUTPUT
