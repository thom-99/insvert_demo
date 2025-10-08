import VariantObjects
import utils
import datetime

# -------------------------------------------
# implement the code below with args from sys
# implement a large number of tests (consider building a separate module for just tests)
# implement the header building as a FUNCTION in UTILS.PY (INCLUDE CHROMSOME/CONTIG LENGTHS) (take as input the chrooms, lengths from read.fai)

realdict = utils.parse_vcf('data/sniffles.vcf')
fakedict = utils.simdict(realdict)

chroms, lengths = utils.read_fai('data/cerevisiae.fa')



with open('data/inSVert.vcf', 'w') as vcf:
    # Write header
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write('##source=inSVert')
    vcf.write(f"##fileDate={datetime.date.today().strftime('%Y%m%d')}\n")
    vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n')
    vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">\n')

    vcf.write('##ALT=<ID=DEL,Description="Deletion">\n')
    vcf.write('##ALT=<ID=INS,Description="Insertion">\n')
    vcf.write('##ALT=<ID=DUP,Description="Duplication">\n')
    vcf.write('##ALT=<ID=INV,Description="Inversion">\n')

    vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">\n')

    
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    
    count = 1

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

                INV = VariantObjects.Deletion(chrom, pos, l, id)
                vcf.write(INV.format() + '\n')

        if svtype == 'DUP':
            for l in fakedict['DUP']['lengths']:    
                chrom, length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, length)


                id = f'inSVert.{svtype}.{count}'
                count += 1 

                DUP = VariantObjects.Deletion(chrom, pos, l, id)
                vcf.write(DUP.format() + '\n')