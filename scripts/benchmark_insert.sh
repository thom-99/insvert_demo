#!/bin/bash
set -euo pipefail

# ==============================================================================
# benchmark_insert.sh - Performance Benchmark Script for inSVert insert
# 
# Runs 'inSVert insert' using 'conda run -n demo', tracking wall-clock time,
# peak memory (RSS), and CPU usage specifically for the insert module.
# ==============================================================================

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <reference.fasta> <variants.vcf> [ploidy|config.yaml] [output_dir]"
    echo ""
    echo "Arguments:"
    echo "  <reference.fasta>      Path to reference FASTA file"
    echo "  <variants.vcf>         Path to VCF file containing structural variants"
    echo "  [ploidy|config.yaml]   Ploidy integer (e.g. 2) OR path to config.yaml (default: 2)"
    echo "  [output_dir]           Optional output directory (default: benchmark_insert_results)"
    exit 1
fi

REF="$1"
VCF="$2"
PLOIDY_ARG="${3:-2}"
OUTDIR="${4:-benchmark_insert_results}"

mkdir -p "$OUTDIR"

# Validate inputs
if [ ! -f "$REF" ]; then
    echo "Error: Reference FASTA file '$REF' not found!" >&2
    exit 1
fi

if [ ! -f "$VCF" ]; then
    echo "Error: VCF file '$VCF' not found!" >&2
    exit 1
fi

# Determine ploidy (from number or config.yaml)
if [ -f "$PLOIDY_ARG" ]; then
    PLOIDY=$(python3 -c "import yaml; cfg=yaml.safe_load(open('$PLOIDY_ARG')); print(cfg.get('genome', {}).get('ploidy', 2))")
    PLOIDY_SOURCE="Extracted from $PLOIDY_ARG"
elif [[ "$PLOIDY_ARG" =~ ^[0-9]+$ ]]; then
    PLOIDY="$PLOIDY_ARG"
    PLOIDY_SOURCE="CLI argument"
else
    echo "Error: Ploidy argument '$PLOIDY_ARG' is neither a valid integer nor a readable config file!" >&2
    exit 1
fi

INS_FASTA="$OUTDIR/simulated.fasta"
LOGFILE="$OUTDIR/insert_performance.log"

echo "================================================================================"
echo "                inSVert insert PERFORMANCE BENCHMARK RUNNER                    "
echo "================================================================================"
echo "Reference   : $REF"
echo "VCF File    : $VCF"
echo "Ploidy Used : $PLOIDY ($PLOIDY_SOURCE)"
echo "Output Dir  : $OUTDIR"
echo "================================================================================"

# Function to parse and format stats from GNU time output log
parse_time_log() {
    local logfile="$1"
    local wall_time=$(grep "Elapsed (wall clock) time" "$logfile" | awk -F': ' '{print $2}')
    local user_time=$(grep "User time (seconds)" "$logfile" | awk -F': ' '{print $2}')
    local sys_time=$(grep "System time (seconds)" "$logfile" | awk -F': ' '{print $2}')
    local cpu_pct=$(grep "Percent of CPU" "$logfile" | awk -F': ' '{print $2}')
    local max_rss_kb=$(grep "Maximum resident set size" "$logfile" | awk -F': ' '{print $2}')
    
    local max_rss_mb=$(awk -v kb="$max_rss_kb" 'BEGIN { printf "%.2f", kb / 1024 }')
    
    echo "  Elapsed Wall Time   : $wall_time"
    echo "  Peak Memory (RSS)   : ${max_rss_mb} MB (${max_rss_kb} KB)"
    echo "  User CPU Time       : ${user_time} s"
    echo "  System CPU Time     : ${sys_time} s"
    echo "  Average CPU Usage   : $cpu_pct"
}

# RUN INSERT MODULE
echo ""
echo "[+] Running inSVert insert (ploidy=$PLOIDY)..."
if ! /usr/bin/time -v conda run -n demo inSVert insert "$REF" "$VCF" --ploidy "$PLOIDY" -o "$INS_FASTA" 2>"$LOGFILE"; then
    echo "[-] Error: 'inSVert insert' failed! Check log below:" >&2
    cat "$LOGFILE" >&2
    exit 1
fi

# PRINT STATS REPORT
echo ""
echo "================================================================================"
echo "                   inSVert insert BENCHMARK RESULTS REPORT                      "
echo "================================================================================"
parse_time_log "$LOGFILE"
echo ""
echo "================================================================================"
echo "Benchmark completed successfully!"
echo "Outputs stored in: $OUTDIR"
echo "  - Inserted FASTA : $INS_FASTA"
echo "  - Timing Log     : $LOGFILE"
echo "================================================================================"
