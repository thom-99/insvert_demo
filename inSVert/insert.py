import sys
import utils2
import pysam

'''
insvert module that performs the insertion of a VCF file into a fasta reference
'''

if len(sys.argv) != 5:
    print("Usage: python insert.py <config.yaml_path> <fasta_path> <vcf_path> <output_fasta>")
    sys.exit(1)

config_path = sys.argv[1]
ref_fasta = sys.argv[2]
vcf_file = sys.argv[3]
output_fasta = sys.argv[4]

# reading configfile
print('reading config file')
target_gc = utils2.get_GC(config_path)

# loading genome
print('loading the refetrence genome...')

genome = utils2.parse_fasta(ref_fasta)

print(genome.keys())
print()

# parsing vcf & applying SVs
print('parsing VCF...')
vcf = pysam.VariantFile(vcf_file)

for SV in vcf:

    chrom = SV.chrom
    pos = SV.pos - 1
    svtype = SV.info.get("SVTYPE")

    chrom_seq = genome[chrom]['sequence']
    chrom_offset = genome[chrom]['offset']

    # --- DELETION ---
    if svtype == 'DEL':
        svlen = SV.info.get("SVLEN")
        
        # applying deletions and updating offset
        genome[chrom]['offset'] = utils2.apply_deletion(chrom_seq, svlen, pos, chrom_offset)

    # --- INSERTION ---
    elif svtype == 'INS':
        svlen =  SV.info.get("SVLEN")
        # generate sequence using the config GC content
        ins_seq = utils2.generate_seq(svlen, target_gc) 

        genome[chrom]['offset'] = utils2.apply_insertion(chrom_seq, ins_seq, pos, chrom_offset)

    # --- INVERSION ---
    elif svtype == 'INV':
        svlen = SV.info.get("SVLEN")
        end = pos + svlen
        genome[chrom]['offset'] = utils2.apply_inversion(chrom_seq, pos, end, chrom_offset)

    # --- DUPLICATION ---
    elif svtype == 'DUP':
        svlen= SV.info.get("SVLEN")

        # EXTRACT COPY NUMBER
        # PySAM stores samples in a dict-like object. 
        # We grab the first sample (index 0) and look up the 'CN' key.
        sample_name = list(SV.samples.keys())[0]
        copy_number = SV.samples[sample_name]['CN']        

        genome[chrom]['offset'] = utils2.apply_duplication(chrom_seq, pos, svlen, copy_number, chrom_offset)


vcf.close()


# writing output 
print(f'writing edited genome to {output_fasta}')

with open(output_fasta, 'w') as fasta:
    for chrom, data in genome.items():
        
        fasta.write(f'>{chrom}\n')

        # decode bytearray as string
        sequence = data['sequence'].decode('ascii')

        for i in range(0, len(sequence), 60):
            fasta.write(sequence[i:i+60] + '\n')

print('Done.')