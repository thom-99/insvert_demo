import vcfparser

'''
from the output dictionary of vcfparser.py it gets the data and
1. fits those to a power-law distribution
2. samples randomly N SVLEN for each SVtype
3. builds a reflection of the input dictionary with the sampled SVLEN
'''

# I can get my dict with vcfparser.parse_vcf('vcfpath')
