'''
MISCELLANEOUS COLLECTION OF FUNCTIONS OF INSVERT
'''
import pysam
import os 
import sys
import yaml
import numpy as np
import datetime
import bisect
import random
from scipy import stats


import yaml
import numpy as np


def calculate_pareto_alpha(median_bp, min_bp):
    """
    compute the alpha (shape) parameter for a Pareto distribution
    given a desired median and a minimum length
    reference formula: Median = min * 2^(1/alpha)
    """
    if median_bp <= min_bp:
        raise ValueError("median must be greater than min_length for a valid Pareto distribution")
    
    alpha = np.log(2) / np.log(median_bp/min_bp)
    return alpha
    

def parse_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    
    #extracting ploidy and heterozygosity values
    genome_settings = config.get('genome')
    ploidy = genome_settings.get('ploidy',2)
    heterozygosity = genome_settings.get('heterozygosity',0.5) 

    sv_data = {}
    
    # Iterate through the variants defined in YAML
    for sv_type, settings in config['variants'].items():
        count = settings['count']
        dist_type = settings['distribution']
        params = settings['parameters']
        
        lengths = []
        
        # GENERATE LENGTHS
        if dist_type == "pareto":
            median_length = params['median_length']
            min_length = params['min_length']
            alpha = calculate_pareto_alpha(median_length, min_length)

            lengths = []
            while len(lengths) < count:
                # np.random.pareto(alpha) generates a type II pareto, so I adjust
                sample = (np.random.pareto(alpha) + 1) * min_length
                val = int(sample)
                if params['min_length'] <= val <= params['max_length']:
                    lengths.append(val)
        
        elif dist_type == "normal":
            mu = params["median_length"]
            sigma = params.get('sigma')
            if sigma is None:
                sigma = mu * 0.1 #default to 10% of the median

            lengths = []
            while len(lengths) < count:
                sample = np.random.normal(mu, sigma)
                val = int(sample)
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
    
    return {'variants' : sv_data,
            'ploidy' : ploidy,
            'heterozygosity' : heterozygosity
            }

#print(parse_config('inSVert/config.yaml'))

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

    header_lines.append('##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">')
    header_lines.append('##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated with breakend">')
    header_lines.append('##INFO=<ID=TRA_ROLE,Number=1,Type=String,Description="Role in translocation: SOURCE or SINK">')

    # ALT fields
    header_lines.append('##ALT=<ID=DEL,Description="Deletion">')
    header_lines.append('##ALT=<ID=INS,Description="Insertion">')
    header_lines.append('##ALT=<ID=DUP,Description="Duplication">')
    header_lines.append('##ALT=<ID=INV,Description="Inversion">')
    header_lines.append('##ALT=<ID=BND,Description="Breakend">')

    # FORMAT fields
    header_lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header_lines.append('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">')

    header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE")

    return '\n'.join(header_lines) + '\n'


'''
FUNCTION TO CHECK WETHER A SV OVERLAPS WITH ANOTHER PRE-EXISTING ONE O(NxL) -- suboptimal solution
the pre-existing SVs are stored in a dictionary {chrom : [(pos, end)]}
'''
def overlaps_suboptimal(chrom, start, end, sv_positions:dict):
    for existing_start, existing_end in sv_positions[chrom]:
        # the intervals do not overlap if the new interval ends before an existing one starts or if starts before the existing one ends
        if not (end < existing_start or start > existing_end):
            return True # it overlaps
    return False



"""
Checks for overlaps using Binary Search. O(log N).
Requires sv_positions[chrom] to be sorted by start position.
Now haplotype-aware: only checks the haplotypes where the variant is placed.
"""
def overlaps(chrom, start, end, genotype_str, sv_positions: dict):

    alleles = genotype_str.split('/')

    for hap_idx, allele in enumerate(alleles):
        if allele == "1": #check only haplotypes that carry the variant

            intervals = sv_positions[chrom][hap_idx]
        
            if not intervals: # base case, empty list == no overlaps 
                return False

            # 1) Binary Search to find the "Insertion Point" 
            # We ask: "If I were to insert this new variant into the list while keeping it sorted, at which index would it go?"
            idx = bisect.bisect_right(intervals, (start, end))

            # 4. Check the Left Neighbor == the variant immediately BEFORE the new one
            if idx > 0:
                prev_start, prev_end = intervals[idx - 1]
                # If the previous variant's END extends past our START, there's overlap.
                if prev_end > start:
                    return True

            # 5. Check the Right Neighbor == the variant immediately AFTER the new one
            if idx < len(intervals):
                next_start, next_end = intervals[idx]
                # If the next variant's START begins before our END, we overlap.
                if next_start < end:
                    return True

    # otherwise there's no overlap
    return False






'''
Generates a genotype string (e.g., '0/1/0') based on ploidy and heterozygosity.
Ensures that at least one allele is '1'
'''
def generate_genotype(ploidy:int, heterozygosity:float) -> str:

    #randomly assign alleles based on heterozygosity prob
    alleles = [1 if random.random() < heterozygosity else 0 for copy in range(ploidy)]

    #if all alleles are 0, we force at least one to be a 1
    if sum(alleles) == 0:
        random_idx = random.randint(0,ploidy-1)
        alleles[random_idx] = 1
    
    #format the output as a VCF genotype string (ex. "0/1")
    return "/".join(map(str,alleles))

