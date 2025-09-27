import vcfparser
from scipy import stats
import numpy as np

'''
testing out different distributions and comparing them
'''

vcf_data = vcfparser.parse_vcf('data/sniffles.vcf')

distributions = ['gamma','lognorm','expon','weibull_min','pareto']
fits = {}

for distr_name in distributions:
    distr = getattr(stats, distr_name)
    params = distr.fit(vcf_data['INS']['lengths'], floc=49)
    ll = np.sum(distr.logpdf(vcf_data['INS']['lengths'], *params))
    fits[distr_name] = (params, ll)

print(sorted(fits.items(), key= lambda x: x[1][1], reverse=True))
