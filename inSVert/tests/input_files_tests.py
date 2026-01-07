import os 
import pysam

# --------------------------------------------------------------------------------
# MODULE TO CHECK VALIDITY OF INPUT FILES OF inSVert
# --------------------------------------------------------------------------------




### FUNCTIONS FOR HANDLING THE INPUT VCF ###

''' removes a file and its associated index (.tbi)'''
def cleanup_indexed_vcf(filepath):
    if filepath and os.path.exists(filepath):
        os.remove(filepath)
    if filepath and os.path.exists(filepath + ".tbi"):
        os.remove(filepath + ".tbi")
import os
import pysam



# 1. Helper: Sorts VCF using Python RAM (No external tools needed)
def sort_vcf_pure_python(input_path, output_path):
    """
    Reads all variants, sorts them in memory, and writes a compressed VCF.
    This replaces 'bcftools sort' and avoids missing library errors.
    """
    # Open the input
    infile = pysam.VariantFile(input_path)
    header = infile.header
    
    # Create a map of Chromosome -> Order (e.g., chr1=0, chr2=1)
    # This ensures we sort exactly as the VCF header defines the chromosomes
    contig_map = {contig: i for i, contig in enumerate(header.contigs)}
    
    # Read all variants into a list
    # (Safe for simulation files with < 500,000 variants)
    try:
        records = list(infile)
    finally:
        infile.close()
    
    # Sort: First by Chromosome ID, then by Position
    records.sort(key=lambda r: (contig_map.get(r.chrom, 9999), r.pos))
    
    # Write to output with 'wz' mode (Force BGZF compression)
    outfile = pysam.VariantFile(output_path, mode='wz', header=header)
    for rec in records:
        outfile.write(rec)
    outfile.close()

# 2. Main Function
def prepare_vcf(vcf_path):
    # Case 1: Already compressed and indexed?
    if vcf_path.endswith(".gz") and os.path.exists(vcf_path + ".tbi"):
        print("Using existing indexed VCF.")
        return vcf_path
    
    # Case 2: We need to sort it. 
    # We define the target name (e.g. data.vcf -> data.vcf.gz)
    sorted_vcf_path = vcf_path + ".gz"
    
    print(f"Sorting and indexing variants to {sorted_vcf_path}...")
    
    try:
        # STEP A: Sort using Python (Cannot fail due to missing libs)
        sort_vcf_pure_python(vcf_path, sorted_vcf_path)
        
        # STEP B: Index using Pysam
        pysam.tabix_index(sorted_vcf_path, preset="vcf", force=True)
        
        return sorted_vcf_path 
    
    except Exception as e:
        # Cleanup if we failed halfway
        if os.path.exists(sorted_vcf_path):
            os.remove(sorted_vcf_path)
        raise RuntimeError(f"Failed to prepare VCF: {e}")