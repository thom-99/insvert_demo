import VariantObjects
import utils


# implement the code below with args from sys

realdict = utils.parse_vcf('data/sniffles.vcf')
fakedict = utils.simdict(realdict)

chroms, lengths = utils.read_fai('data/cerevisiae.fa')


lines = []

for svtype in fakedict:
    if svtype == 'INS':
        for l in fakedict['INS']['lengths']:
            chrom, length = utils.select_chr(chroms, lengths)
            pos = utils.select_pos(chrom, length)

            count=1
            id = f'inSVert.{svtype}.{count}'

            INS = VariantObjects.Insertion(chrom, pos, l, id)
            lines.append(INS.format())
    
    if svtype == 'DEL':
        for l in fakedict['DEL']['lengths']:    
            chrom, length = utils.select_chr(chroms, lengths)
            pos = utils.select_pos(chrom, length)

            count=1
            id = f'inSVert.{svtype}.{count}'

            DEL = VariantObjects.Deletion(chrom, pos, l, id)
            lines.append(DEL.format())

print(lines)