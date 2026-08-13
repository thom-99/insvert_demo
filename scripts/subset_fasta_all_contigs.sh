#!/usr/bin/env bash
# Create a small FASTA subset containing sequence from every input contig.
# Usage: subset_fasta_all_contigs.sh INPUT.fa OUTPUT.fa [MAX_BYTES]

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 INPUT.fa OUTPUT.fa [MAX_BYTES]" >&2
    exit 2
fi

input=$1
output=$2
max_bytes=${3:-50000000}

if [[ ! -f "$input" ]]; then
    echo "Error: input FASTA not found: $input" >&2
    exit 1
fi
if [[ -e "$output" ]]; then
    echo "Error: output already exists: $output" >&2
    exit 1
fi
if [[ ! "$max_bytes" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: MAX_BYTES must be a positive integer." >&2
    exit 1
fi

output_dir=$(dirname "$output")
if [[ ! -d "$output_dir" ]]; then
    echo "Error: output directory does not exist: $output_dir" >&2
    exit 1
fi

lengths=$(mktemp)
selection=$(mktemp)
temp_output=$(mktemp "$output_dir/.subset_fasta.XXXXXX")
trap 'rm -f -- "$lengths" "$selection" "$temp_output"' EXIT

# Record sequence lengths while retaining the FASTA identifier (the first word
# after '>'). Empty or duplicate identifiers are rejected because they would
# create an invalid subset FASTA.
awk '
    /^>/ {
        header = substr($0, 2)
        split(header, fields, /[[:space:]]+/)
        name = fields[1]
        if (name == "") {
            print "Error: FASTA header without a contig identifier." > "/dev/stderr"
            failed = 1
            next
        }
        if (seen[name]++) {
            print "Error: duplicate FASTA contig identifier: " name > "/dev/stderr"
            failed = 1
            next
        }
        names[++count] = name
        lengths[name] = 0
        current = name
        next
    }
    {
        if (current == "") {
            print "Error: sequence found before the first FASTA header." > "/dev/stderr"
            failed = 1
            next
        }
        gsub(/[[:space:]]/, "")
        lengths[current] += length($0)
    }
    END {
        if (count == 0) {
            print "Error: input contains no FASTA records." > "/dev/stderr"
            failed = 1
        }
        for (i = 1; i <= count; i++) {
            if (lengths[names[i]] == 0) {
                print "Error: empty FASTA contig: " names[i] > "/dev/stderr"
                failed = 1
            }
        }
        if (failed) exit 1
        for (i = 1; i <= count; i++) print names[i] "\t" lengths[names[i]]
    }
' "$input" > "$lengths"

contig_count=$(wc -l < "$lengths")
header_bytes=$(awk -F '\t' '{ total += length($1) + 2 } END { print total }' "$lengths")
sequence_budget=$((max_bytes - header_bytes))

if (( sequence_budget <= contig_count )); then
    echo "Error: MAX_BYTES is too small to include sequence from all $contig_count contigs." >&2
    exit 1
fi

# Every emitted base takes one byte; FASTA wrapping adds one newline per 60
# bases. Reserve one additional newline per contig to keep the final file
# safely within MAX_BYTES.
per_contig_bases=$(( (sequence_budget - contig_count) * 60 / (contig_count * 61) ))
if (( per_contig_bases < 1 )); then
    echo "Error: MAX_BYTES is too small to include sequence from all $contig_count contigs." >&2
    exit 1
fi

awk -F '\t' -v quota="$per_contig_bases" '{
    bases = ($2 < quota) ? $2 : quota
    print $1 "\t" bases
}' "$lengths" > "$selection"

# A second streaming pass writes the first allocated bases from every contig.
# It normalizes output to 60-base lines, which makes the size calculation above
# exact for the emitted sequence.
awk -v selection="$selection" '
    function flush_full_lines() {
        while (length(buffer) >= 60) {
            print substr(buffer, 1, 60)
            buffer = substr(buffer, 61)
        }
    }
    function finish_contig() {
        flush_full_lines()
        if (length(buffer) > 0) print buffer
        buffer = ""
    }
    BEGIN {
        while ((getline line < selection) > 0) {
            split(line, fields, "\t")
            take[fields[1]] = fields[2]
        }
        close(selection)
    }
    /^>/ {
        if (selected) finish_contig()
        header = substr($0, 2)
        split(header, fields, /[[:space:]]+/)
        name = fields[1]
        remaining = take[name]
        selected = remaining > 0
        if (selected) print ">" name
        next
    }
    selected && remaining > 0 {
        sequence = $0
        gsub(/[[:space:]]/, "", sequence)
        piece = substr(sequence, 1, remaining)
        buffer = buffer piece
        remaining -= length(piece)
        flush_full_lines()
    }
    END {
        if (selected) finish_contig()
    }
' "$input" > "$temp_output"

actual_bytes=$(wc -c < "$temp_output")
if (( actual_bytes > max_bytes )); then
    echo "Error: generated subset is ${actual_bytes} bytes, exceeding MAX_BYTES=${max_bytes}." >&2
    exit 1
fi

mv -- "$temp_output" "$output"
echo "Wrote $output (${actual_bytes} bytes; ${contig_count} contigs; up to ${per_contig_bases} bases per contig)."
