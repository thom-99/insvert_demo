#!/bin/bash 
set -euo pipefail


# === quick check ===
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference.fasta> <configfile.yaml> <oudir>"
    exit 1
fi

#activate env

#variables setup
ref="$1"
config="$2"
mkdir -p "$3"
outdir="$3"
ploidy=$(python -c "import yaml; cfg=yaml.safe_load(open('$config')); print(cfg['genome']['ploidy'])")

# 1: SIMULATE
conda run -n demo inSVert simulate $config $ref -o "$outdir/simulated.vcf"

# 2. INSERT
conda run -n demo inSVert insert $ref "$outdir/simulated.vcf" --ploidy $ploidy -o "$outdir/simulated.fa"

# 3. GENERATE READS
badread simulate --reference "$outdir/simulated.fa" --quantity 40x --error_model random \
    --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
    --identity 30,3 --length 40000,20000 --start_adapter_seq "" --end_adapter_seq "" \
    2>"$outdir/badread_stderr.log" \
    | gzip > "$outdir/simulated_reads.fastq.gz"

# 4. MAP
conda run -n demo samtools faidx $ref
conda run -n demo minimap2 -d "$outdir/reference.mmi" $ref
conda run -n demo minimap2 -x map-pb -t 8 -a $ref "$outdir/simulated_reads.fastq.gz" \
    | samtools view -bS -@8 - \
    | samtools sort -@8 -o "$outdir/sorted.bam" -

conda run -n demo samtools index -@8 "$outdir/sorted.bam"

# 5 VARIANT CALL
conda run -n demo sniffles -i "$outdir/sorted.bam" -v "$outdir/sniffles.vcf" --allow-overwrite




