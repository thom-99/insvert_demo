'''
MISCELLANEOUS COLLECTION OF FUNCTIONS OF INSVERT
'''
import pysam
import os 
import sys
import yaml
import numpy as np
import datetime
from scipy import stats


import yaml
import numpy as np

def calculate_lognormal_params(mean_bp, sigma_bp):
    """
    Helper: Converts arithmetic mean and standard deviation (in bp) 
    to lognormal underlying parameters (mu and sigma).
    """
    variance = sigma_bp ** 2
    sigma = np.sqrt(np.log(variance / (mean_bp ** 2) + 1))
    mu = np.log(mean_bp) - (sigma ** 2 / 2)
    return mu, sigma

def parse_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    
    sv_data = {}
    
    # Iterate through the variants defined in YAML
    for sv_type, settings in config['variants'].items():
        count = settings['count']
        dist_type = settings['distribution']
        params = settings['parameters']
        
        lengths = []
        
        # GENERATE LENGTHS
        if dist_type == "lognormal":
            
            user_mean = params['mean_length']
            user_sigma = params.get('sigma', user_mean * 0.5) # Default sigma to half of mean if missing
            
            mu, sigma = calculate_lognormal_params(user_mean, user_sigma)
            
            lengths = []
            while len(lengths) < count:
                
                # sample 1 obs
                sample = np.random.lognormal(mu, sigma)
                val = int(sample)
    
                # check if between min and max length
                if params['min_length'] <= val <= params['max_length']:
                    lengths.append(val)

            
        sv_data[sv_type] = {'lengths': lengths}
        

        # GENERATE COPY NUMBERS (specific to DUP)
        if sv_type == 'DUP' and 'copy_number' in settings:
            cn_settings = settings['copy_number']
            cn_min = cn_settings['min']
            cn_max = cn_settings['max']
            weights = cn_settings.get('weights')

            # Create array of possible values (e.g. [2, 3, 4, 5])
            possible_cns = np.arange(cn_min, cn_max + 1)
            
            # check and normalize weights to sum 1
            probs = None
            if weights:
                if len(weights) != len(possible_cns):
                    raise ValueError(
                        f"Config Error: DUP weights count ({len(weights)}) "
                        f"does not match copy number range {cn_min}-{cn_max} ({len(possible_cns)} values)."
                    )
                
                weights = np.array(weights)
                probs = weights / weights.sum()
            
            # Sample Copy Numbers based on weights
            # choice() picks 'count' items from 'possible_cns' using 'probs'
            copy_numbers = np.random.choice(possible_cns, size=count, p=probs).tolist()
            
            sv_data[sv_type]['copy_numbers'] = copy_numbers
    
    return sv_data

print(parse_config('inSVert/config.yaml'))

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
            print(f'index file created successfully : {fai_file}', file=sys.stderr)
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
        header_lines.append(f"##reference={reference_path}")
    
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