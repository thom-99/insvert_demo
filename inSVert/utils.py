'''
MISCELLANEOUS COLLECTION OF FUNCTIONS OF INSVERT
'''
import pysam
import os 
import sys
import numpy as np
import datetime
from scipy import stats


'''
PARSES A VCF FILE WITH PYSAM AND INSERTS INFORMATION INTO A DICTIONARY
'''
def parse_vcf(vcf_path:str):
    #parses a vcf and inserts information into a dictonary
    vcf = pysam.VariantFile(vcf_path)
    
    sv_data = {
        'INS': {
            'lengths': [],
        },
        'DEL': {
            'lengths': [],
        },
        'DUP': {
            'lengths': [],
            'copy_numbers': [],  
        },
        'INV': {
            'lengths': [],
        }
    }

    for SV in vcf:
        svtype = SV.info.get("SVTYPE")
        if svtype in sv_data:
            # get SV length
            sv_length = SV.info.get("SVLEN")
            if sv_length is not None:
                # SVLEN can be negative for deletions, take absolute value
                sv_length = abs(sv_length)
                sv_data[svtype]['lengths'].append(sv_length)
            
            # for duplications, get also copy number
            if svtype == 'DUP':
                if "CN" in SV.info:
                    copy_number = SV.info.get("CN")
                    if copy_number is not None:
                        sv_data[svtype]['copy_numbers'].append(copy_number)

    vcf.close()
    return sv_data


'''
PRINTS OPTIONAL SUMMARY OF THE ABOVE
'''
def print_summary(sv_data:dict):
    # optional summary statistics for parsed SV
    print("Structural Variant Summary:")
    print("-" * 30)
    for sv_type in sv_data:
        count = len(sv_data[sv_type]['lengths'])
        print(f"{sv_type}: {count} variants")
        if count > 0:
            lengths = sv_data[sv_type]['lengths']
            print(f"  Length range: {min(lengths)} - {max(lengths)} bp")
            print(f"  Mean length: {int(sum(lengths) / len(lengths))} bp")
            
            if sv_type == 'DUP' and sv_data[sv_type]['copy_numbers']:
                cn_values = sv_data[sv_type]['copy_numbers']
                print(f"  Copy numbers: {min(cn_values)} - {max(cn_values)}")
        print()
    print("-" * 30)


'''
-------------------------------------------------------------------------------------------------
'''

'''
FITS AN ARRAY OF NUMBERS TO A LOGNORMAL DISTRIBUTION
returns the distribution and the resulting estimated parameters
'''
def fit(data:list):
    distr = getattr(stats,'lognorm')
    params = distr.fit(data)
    return distr, params 

def sample(distr, params:tuple, n:int):
    # sample n values from a fitted distribution
    samples = distr.rvs(*params, size=n)
    samples = samples.astype(int)
    return samples[:n].tolist()

'''
PRODUCES A SIMULATED DICTIONARY WITH LENGTHS / COPY-NUMBERS SAMPLED FROM THE LOGNORMAL
the input is the dictionary of real values extracted from the VCF,
it fits the values for each field of each svtype and samples from it.
'''
def simdict(realdict:dict):
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
                fakevalues = values[::-1]
                fakedict[svtype][field] = fakevalues
        
    return fakedict

'''
-----------------------------------------------------------------------------------------------
'''

'''
IF IT EXISTS, READS A FAI FILE FROM THE PATH OF THE FASTA; ELSE IT CREATES IT WITH PYSAM
returns a tuple of all the ([chromosomes],[lengths])
'''
def read_fai(fasta_path):
 
    #check suffix 
    if fasta_path.endswith('.fasta') or fasta_path.endswith('.fa'):
        fasta_file = fasta_path
        fai_file = fasta_path + '.fai'
    else:
        raise ValueError(f'input file must be a fasta file ending with .fa or .fasta : {fasta_path}')

    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f'fasta file not found : {fasta_file}')
    
    if not os.path.exists(fai_file):
        print(f'index file not found, creating index for {fasta_file}')
        try:
            pysam.faidx(fasta_file)
            print(f'index file created successfully : {fai_file}', file=sys.sterr)
        except Exception as e:
            raise RuntimeError(f'failed to create fasta index file: {e}')

    # read the .fai file and process info

    chroms = []
    lengths = []

    with open(fai_file) as fai:
        for line in fai:
            fields = line.strip().split('\t')
            chrom = fields[0]
            length = int(fields[1])
            chroms.append(chrom)
            lengths.append(length)

    print(f"Loaded {len(chroms)} chromosomes from {fai_file}", file=sys.stderr)

    return chroms, lengths

'''
PICKS A CHROMOSOME RANDOMLY TAKING INTO ACCOUNT THE LENGTH: WEIGHTED RANDOM SELECTION
outputs the chromosome and the respective length
'''
def select_chr(chroms, lengths):

    weights = np.array(lengths)
    probabilities = weights / weights.sum()

    chrom = np.random.choice(chroms, p=probabilities)
    chrom_idx = chroms.index(chrom)

    return chrom, lengths[chrom_idx]    

'''
SELECTS A RANDOM POSITION ALONG THE SPECIFIED CHROMSOME
optionally it can avoid choosing positions too close to the chr end
'''
def select_pos(chrom, length, buffer=1000):

    if buffer > length/10:
        raise ValueError(f'edge buffer :{buffer} is too big, try decreasing it')
    
    position = np.random.randint(1 + buffer, length - buffer + 1)

    return position



'''
BUILDS A COMPLETE HEADER FOR A VCF FILE
given the chroms, lengths (easily accessible through read_fai it also builds contigs lines)
'''
def buildheader(chroms, lengths, reference_path=None):

    header_lines = []

    header_lines.append("##fileformat=VCFv4.2")
    header_lines.append("##source=inSVert")
    header_lines.append(f"##fileDate={datetime.date.today().strftime('%Y%m%d')}")

    # optional reference path
    if reference_path:
        header_lines.append(f"##fileDate={datetime.date.today().strftime('%Y%m%d')}")
    
    # contig lines 
    for chrom, length in zip(chroms, lengths):
        header_lines.append(f"##contig=<ID={chrom},length={length}>")    

    # INFO fields
    header_lines.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header_lines.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">')
    header_lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">')

    # ALT fields
    header_lines.append('##ALT=<ID=DEL,Description="Deletion">')
    header_lines.append('##ALT=<ID=INS,Description="Insertion">')
    header_lines.append('##ALT=<ID=DUP,Description="Duplication">')
    header_lines.append('##ALT=<ID=INV,Description="Inversion">')

    # FORMAT fields
    header_lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header_lines.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">')

    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")

    return '\n'.join(header_lines) + '\n'


'''
FUNCTION TO CHECK WETHER A SV OVERLAPS WITH ANOTHER PRE-EXISTING ONE
the pre-existing SVs are stored in a dictionary {chrom : [(pos, end)]}
'''
def overlaps(chrom, start, end, sv_positions:dict):
    for existing_start, existing_end in sv_positions[chrom]:
        # the intervals do not overlap if the new interval ends before an existing one starts or if starts before the existing one ends
        if not (end < existing_start or start > existing_end):
            return True # it overlaps
    return False