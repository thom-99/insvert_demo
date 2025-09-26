#!/bin/bash

# === usage check === 
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.bam output.vcf"
    exit 1
fi

# === input arguments ===
BAM=$(realpath "$1")
OUTPUT_VCF="$2"

conda run -n demo sniffles -i "$BAM" -v "$OUTPUT_VCF" --allow-overwrite
























