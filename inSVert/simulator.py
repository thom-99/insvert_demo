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


# FIX, DOES NOT PRODODUCE THE OTHER SVTYPES EXCEPT INS

def simdict(realdict:dict):
    # builds a simulated dictionary pulling lengths / copy numbers from fitted distributions
    fakedict = {}

    #looping through the first dict and applying fit & sample to the values
    for svtype, fields in realdict.items():
        fakedict[svtype] = {}

        for field, values in fields.items():
            n = len(values)
            if n > 5:
                # if there are sufficient entries, fit
                distr, params = fit(values)
                fakevalues = sample(distr, params, n)
                fakedict[svtype][field] = fakevalues
            else:
                # keep the original values but in reverse order
                fakevalues = reversed(values)
                fakedict[svtype][field] = fakevalues
        
    return fakedict



        

#checks

vcfdata = vcfparser.parse_vcf('data/sniffles.vcf')
print(vcfdata)
#insdata = vcfdata['INS']['lengths']

#distr, params = fit(insdata)
#samples = sample(distr, params, 10000)

#print(max(samples))
#print(min(samples))
#print(sum(samples)/len(samples))

fakedict = simdict(vcfdata)
print(fakedict)
