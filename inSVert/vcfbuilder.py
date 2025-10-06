import VariantObjects
import utils

# implement the code below with args from sys

realdict = utils.parse_vcf('data/sniffles.vcf')
fakedict = utils.simdict(realdict)

chroms, lengths = utils.read_fai('data/cerevisiae.fa')


for svtype in fakedict:
    if svtype == 'INS':
        for l in fakedict['INS']['lengths']:
