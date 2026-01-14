#!/bin/bash

# Check if VCF file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <vcf_file>"
    echo "Example: $0 variants.vcf"
    exit 1
fi

VCF_FILE=$1

# Check if file exists
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: File $VCF_FILE not found!"
    exit 1
fi

echo "Counting variant types in $VCF_FILE"
echo "===================================="

# Count each variant type (excluding header lines)
INS_COUNT=$(grep -v "^#" "$VCF_FILE" | grep -c "INS")
DEL_COUNT=$(grep -v "^#" "$VCF_FILE" | grep -c "DEL")
DUP_COUNT=$(grep -v "^#" "$VCF_FILE" | grep -c "DUP")
INV_COUNT=$(grep -v "^#" "$VCF_FILE" | grep -c "INV")

# Display results
echo "INS: $INS_COUNT"
echo "DEL: $DEL_COUNT"
echo "DUP: $DUP_COUNT"
echo "INV: $INV_COUNT"
