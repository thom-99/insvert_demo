import VariantObjects
import utils
import datetime
import sys

# -------------------------------------------
# implement the code below with args from sys
# implement a large number of tests (consider building a separate module for just tests)
# implement the header building as a FUNCTION in UTILS.PY (INCLUDE CHROMSOME/CONTIG LENGTHS) (take as input the chrooms, lengths from read.fai)
if len(sys.argv) != 4:
    print("Usage: python script.py <vcf_path> <fasta_path> <output_file>")
    sys.exit(1)

vcf_path = sys.argv[1]
fasta_path = sys.argv[2]
output_file = sys.argv[3]




realdict = utils.parse_vcf(vcf_path)
fakedict = utils.simdict(realdict)

chroms, lengths = utils.read_fai(fasta_path)

with open(output_file, 'w') as vcf:
 
    header = utils.buildheader(chroms, lengths, fasta_path)
    vcf.write(header)

    count = 1

    #start to build SVs programmatically
    for svtype in fakedict:
        if svtype == 'INS':
            for l in fakedict['INS']['lengths']:
                chrom, length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, length)

                id = f'inSVert.{svtype}.{count}'
                count += 1

                INS = VariantObjects.Insertion(chrom, pos, l, id)
                vcf.write(INS.format() + '\n')
    
        if svtype == 'DEL':
            for l in fakedict['DEL']['lengths']:    
                chrom, length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, length)


                id = f'inSVert.{svtype}.{count}'
                count += 1 

                DEL = VariantObjects.Deletion(chrom, pos, l, id)
                vcf.write(DEL.format() + '\n')

        if svtype == 'INV':
            for l in fakedict['INV']['lengths']:    
                chrom, length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, length)


                id = f'inSVert.{svtype}.{count}'
                count += 1 

                INV = VariantObjects.Inversion(chrom, pos, l, id)
                vcf.write(INV.format() + '\n')

        if svtype == 'DUP':
            for l in fakedict['DUP']['lengths']:    
                chrom, length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, length)


                id = f'inSVert.{svtype}.{count}'
                count += 1 

                DUP = VariantObjects.Duplication(chrom, pos, l, id)
                vcf.write(DUP.format() + '\n')


