import vcfparser
from scipy import stats

'''
from the output dictionary of vcfparser.py it gets the data and
1. fits those to a power-law distribution
2. samples randomly N SVLEN for each SVtype
3. builds a reflection of the input dictionary with the sampled SVLEN
'''

# I can get my dict with vcfparser.parse_vcf('vcfpath')
# lengths are stored as a list in the dict dict['INS']['lenghts']

def fit(data:list):
    distr = getattr(stats,'lognorm')
    params = distr.fit(data)
    return distr, params

    

def sample(distr, params:tuple, n:int):
    # sample n values from a fitted distribution
    samples = distr.rvs(*params, size=n)
    samples = samples.astype(int)
    return samples[:n].tolist()


vcfdata = vcfparser.parse_vcf('data/sniffles.vcf')
insdata = vcfdata['INS']['lengths']

distr, params = fit(insdata)
samples = sample(distr, params, 10000)

print(max(samples))
print(min(samples))
print(sum(samples)/len(samples))