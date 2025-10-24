import sys
import utils2
import pysam

'''
insvert module that performs the insertion of a VCF file into a fasta reference
'''

if len(sys.argv) != 4:
    print("Usage: python insert.py <fasta_path> <vcf_path> <output_fasta>")
    sys.exit(1)

ref_fasta = sys.argv[1]
vcf_file = sys.argv[2]
output_fasta = sys.argv[3]


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
    id = SV.info.get("SVTYPE")

    chrom_seq = genome[chrom]['sequence']
    chrom_offset = genome[chrom]['offset']

    print(f"Before: {chrom} length = {len(chrom_seq)}")  # DEBUG

    if id == 'DEL':
        svlen = SV.info.get("SVLEN")
        
        # applying deletiona and updating offset
        genome[chrom]['offset'] = utils2.apply_deletion(chrom_seq, svlen, pos, chrom_offset)

    print(f"After: {chrom} length = {len(chrom_seq)}")  # DEBUG


vcf.close()


with open(output_fasta, 'w') as fasta:
    for chrom, data in genome.items():
        
        fasta.write(f'>{chrom}\n')

        # decode bytearray as string
        sequence = data['sequence'].decode('ascii')

        for i in range(0, len(sequence), 60):
            fasta.write(sequence[i:i+60] + '\n')