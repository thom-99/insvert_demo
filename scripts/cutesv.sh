#!/bin/bash

# === usage check === 
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input.bam input.fasta output.vcf"
    exit 1
fi

# === input arguments ===
BAM=$(realpath "$1")
FASTA=$(realpath "$2")
OUTPUT_VCF="$3"

mkdir -p cutesv_workingdir

conda run -n demo cuteSV --threads 8 "$BAM" "$FASTA" "$OUTPUT_VCF" cutesv_workingdir

rm -rf cutesv_workingdir
