import os 
import pysam


### INPUT VCF VALIDATION ###

def validate_vcf(vcf_path, fasta_path, target_ploidy=None):
    """
    Validates the input VCF for structural integrity, compatibility with the reference FASTA,
    and genotype ploidy.
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"Input VCF file not found at: {vcf_path}")
    
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Reference FASTA file not found at: {fasta_path}")

    try:
        # Load both files
        with pysam.VariantFile(vcf_path) as vcf, pysam.FastaFile(fasta_path) as fasta:
            vcf_contigs = vcf.header.contigs
            fasta_chroms = set(fasta.references)

            # Check 1: Do the VCF chromosomes exist in the FASTA?
            for vcf_chrom in vcf_contigs:
                if vcf_chrom not in fasta_chroms:
                    raise ValueError(f"Chromosome '{vcf_chrom}' found in VCF but missing from reference FASTA. "
                                     f"Are you sure the VCF was called against this exact FASTA?")

            # Check 2: Origin FASTA matching (Compare contig lengths)
            for vcf_chrom, contig_record in vcf_contigs.items():
                if contig_record.length is not None:
                    fasta_length = fasta.get_reference_length(vcf_chrom)
                    if contig_record.length != fasta_length:
                        raise ValueError(
                            f"Length mismatch for '{vcf_chrom}'. "
                            f"VCF expects {contig_record.length}bp, but FASTA has {fasta_length}bp. "
                            f"Are you sure the VCF was called against this exact FASTA?"
                        )

            if target_ploidy is not None:
                for rec in vcf:

                    # Check 3: presence of a sample 
                    if not rec.samples:
                        raise ValueError(f"VCF record at {rec.chrom}:{rec.pos} has no sample genotype.")

                    # Check 4: presence of a genotype
                    gt = rec.samples[0].get('GT')
                    if gt is None or len(gt) == 0:
                        raise ValueError(f"VCF record at {rec.chrom}:{rec.pos} has no GT value")

                    # Check 5: ploidy matching with target_ploidy
                    if len(gt) != target_ploidy:
                        raise ValueError(f"Ploidy mismatch at {rec.chrom}:{rec.pos}, GT={gt} has {len(gt)} alleles but --ploidy {target_ploidy} was specified")
                        

    except Exception as e:
        raise ValueError(f"VCF Validation Failed: {e}")




### INPUT FASTA VALIDATION ###

def validate_fasta(fasta_path):
    """Ensures the reference FASTA exists and can be parsed."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Reference FASTA file not found at: {fasta_path}")
    try:
        with pysam.FastaFile(fasta_path) as fasta:
            if not fasta.references:
                raise ValueError(f"FASTA file '{fasta_path}' contains no sequences.")
    except Exception as e:
        raise ValueError(f"Failed to parse FASTA file. Details: {e}")
    return True


### INPUT BED VALIDATION ###

def validate_bed(bed_path, fasta_path=None):
    """Validates BED format (3 tab-separated columns) and alerts if chroms mismatch."""
    if not bed_path:
        return True
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED exclusion file not found at: {bed_path}")

    fasta_chroms = set()
    if fasta_path and os.path.exists(fasta_path):
        with pysam.FastaFile(fasta_path) as fasta:
            fasta_chroms = set(fasta.references)

    with open(bed_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue
            columns = line.strip().split('\t')
            if len(columns) < 3:
                raise ValueError(f"Malformed BED line {line_num}: Must have 3 tab-separated columns.")
            
            chrom, start, end = columns[0], columns[1], columns[2]
            if not start.isdigit() or not end.isdigit():
                raise ValueError(f"Malformed BED line {line_num}: Coordinates must be integers.")
            if int(start) >= int(end):
                raise ValueError(f"Malformed BED line {line_num}: Start coordinate cannot be >= End.")
            if fasta_chroms and chrom not in fasta_chroms:
                print(f"[Warning] BED line {line_num}: Chromosome '{chrom}' not in reference FASTA.")
    return True