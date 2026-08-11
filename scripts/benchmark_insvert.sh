#!/bin/bash
set -euo pipefail

# ==============================================================================
# benchmark_insvert.sh - Performance Benchmark Script for inSVert
# 
# Runs 'inSVert simulate' and 'inSVert insert' using 'conda run -n demo',
# tracking time, peak memory (RSS), and CPU usage for each module.
# Defaults --ploidy for 'inSVert insert' to the value assigned in config.yaml.
# ==============================================================================

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <config.yaml> <reference.fasta> [ploidy] [output_dir]"
    echo ""
    echo "Arguments:"
    echo "  <config.yaml>      Path to inSVert YAML config file"
    echo "  <reference.fasta>  Path to reference FASTA file"
    echo "  [ploidy]           Optional ploidy override (defaults to genome.ploidy in config.yaml)"
    echo "  [output_dir]       Optional output directory (default: benchmark_results)"
    exit 1
fi

CONFIG="$1"
REF="$2"
USER_PLOIDY="${3:-}"
OUTDIR="${4:-benchmark_results}"

mkdir -p "$OUTDIR"

# Validate inputs
if [ ! -f "$CONFIG" ]; then
    echo "Error: Config file '$CONFIG' not found!" >&2
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "Error: Reference file '$REF' not found!" >&2
    exit 1
fi

# Extract default ploidy from config.yaml if not explicitly supplied
if [ -n "$USER_PLOIDY" ]; then
    PLOIDY="$USER_PLOIDY"
    PLOIDY_SOURCE="CLI argument"
else
    PLOIDY=$(python3 -c "import yaml; cfg=yaml.safe_load(open('$CONFIG')); print(cfg['genome']['ploidy'])")
    PLOIDY_SOURCE="config.yaml"
fi

SIM_VCF="$OUTDIR/simulated.vcf"
INS_FASTA="$OUTDIR/simulated.fasta"

SIM_LOG="$OUTDIR/simulate_time.log"
INS_LOG="$OUTDIR/insert_time.log"

echo "================================================================================"
echo "                   inSVert PERFORMANCE BENCHMARK RUNNER                         "
echo "================================================================================"
echo "Config File : $CONFIG"
echo "Reference   : $REF"
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
    echo "  System CPU Time     : ${sys_cpu} s"
    echo "  Average CPU Usage   : $cpu_pct"
}

# 1. RUN SIMULATE
echo ""
echo "[+] Step 1: Running inSVert simulate..."
/usr/bin/time -v conda run -n demo inSVert simulate "$CONFIG" "$REF" --seed 1234 -o "$SIM_VCF" 2>"$SIM_LOG"

# 2. RUN INSERT
echo ""
echo "[+] Step 2: Running inSVert insert (ploidy=$PLOIDY)..."
/usr/bin/time -v conda run -n demo inSVert insert "$REF" "$SIM_VCF" --ploidy "$PLOIDY" -o "$INS_FASTA" 2>"$INS_LOG"

# 3. PRINT STATS REPORT
echo ""
echo "================================================================================"
echo "                       inSVert BENCHMARK RESULTS REPORT                        "
echo "================================================================================"
echo ""
echo "--------------------------------------------------------------------------------"
echo " 1. MODULE: inSVert simulate"
echo "--------------------------------------------------------------------------------"
parse_time_log "$SIM_LOG"

echo ""
echo "--------------------------------------------------------------------------------"
echo " 2. MODULE: inSVert insert"
echo "--------------------------------------------------------------------------------"
parse_time_log "$INS_LOG"

echo ""
echo "================================================================================"
echo "Benchmark completed successfully!"
echo "Outputs stored in: $OUTDIR"
echo "  - Simulated VCF  : $SIM_VCF"
echo "  - Inserted FASTA : $INS_FASTA"
echo "================================================================================"
