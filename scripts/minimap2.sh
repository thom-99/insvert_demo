#!/bin/bash

# === quick check ===
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference.fasta> <reads_1.fastx> <oudir>"
    exit 1
fi

# === inputs ===
REFERENCE_FA=$(realpath "$1")
READS1=$(realpath "$2")
mkdir -p "$3"
OUTDIR=$(realpath "$3")



#running minmap2 for long reads, forming a sorted & index bam file

conda run -n demo samtools faidx "$REFERENCE_FA"
conda run -n demo minimap2 -d "${OUTDIR}/reference.mmi" "$REFERENCE_FA"
conda run -n demo bash -c "
minimap2 -x map-pb -t 8 -a '$REFERENCE_FA' '$READS1' \
    | samtools view -bS -@8 - \
    | samtools sort -@8 -o '${OUTDIR}/sorted.bam' -
"
conda run -n demo samtools index -@8 "${OUTDIR}/sorted.bam"





