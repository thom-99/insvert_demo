'''
MISCELLANEOUS COLLECTION OF FUNCTIONS OF INSVERT
'''
import pysam
import os 
import sys
import yaml
import datetime
import bisect
import math
import random
from . import VariantObjects
from .utils_ins import reverse_complement



def calculate_pareto_alpha(median_bp, min_bp):
    """
    compute the alpha (shape) parameter for a Pareto distribution
    given a desired median and a minimum length
    reference formula: Median = min * 2^(1/alpha)
    """
    if median_bp <= min_bp:
        raise ValueError("median must be greater than min_length for a valid Pareto distribution")
    
    alpha = math.log(2) / math.log(median_bp/min_bp)
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
                # paretovariate(a) returns a type I pareto where min_val = 
                sample = random.paretovariate(alpha) * min_length
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
                sample = random.gauss(mu, sigma)
                val = int(sample)
                if params['min_length'] <= val <= params['max_length']:
                    lengths.append(val)

           

            
        sv_data[sv_type] = {'lengths': lengths}

        # PARSE REVERSE_RATIO (specific to TRA_CUT / TRA_COPY)
        if sv_type in ('TRA_CUT', 'TRA_COPY'):
            rr = settings.get('reverse_ratio')
            if rr is None and isinstance(params, dict):
                rr = params.get('reverse_ratio')
            if rr is not None:
                if not (0.0 <= rr <= 1.0):
                    raise ValueError(
                        f"Config Error: '{sv_type}.reverse_ratio' must be between 0.0 and 1.0, got {rr}"
                    )
                sv_data[sv_type]['reverse_ratio'] = rr
        

        # GENERATE COPY NUMBERS (specific to DUP)
        if sv_type == 'DUP' and 'copy_number' in settings:
            cn_settings = settings['copy_number']
            cn_min = cn_settings['min']
            cn_max = cn_settings['max']
            weights = cn_settings.get('weights')

            # Create array of possible values (e.g. [2, 3, 4, 5])
            possible_cns = list(range(cn_min, cn_max + 1))
            
            # check and normalize weights to sum 1
            probs = None
            if weights:
                if len(weights) != len(possible_cns):
                    raise ValueError(
                        f"Config Error: DUP weights count ({len(weights)}) "
                        f"does not match copy number range {cn_min}-{cn_max} ({len(possible_cns)} values)."
                    )
                
            
            # Sample Copy Numbers based on weights
            copy_numbers = random.choices(possible_cns, weights=weights, k=count)
            
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
def read_fai(fasta_path:str):
 
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


def parse_bed(bed_path:str):
    """
    Parses a BED file into a dictionary: {chrom: [(start, end), ...]}
    Intervals are stored as sorted tuples for binary search.
    """
    excluded_ranges = {}
    if not bed_path or not os.path.exists(bed_path):
        return None

    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')): continue
            parts = line.strip().split('\t')
            # to add: print a warning to the screen that the line is not properly formatted
            if len(parts) != 3: continue
            
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom not in excluded_ranges:
                excluded_ranges[chrom] = []
            excluded_ranges[chrom].append((start, end))

    # Sort intervals by start position for bisect to work
    for chrom in excluded_ranges:
        excluded_ranges[chrom].sort()
    
    return excluded_ranges


def overlaps_excluded_region(chrom, start, end, excluded_regions:dict):
    """
    Checks if a proposed SV interval overlaps with any region in the BED file.
    Uses O(log N) binary search.
    The structure of this is like a base for the overlaps() function, which is haplotype aware
    and is designed to keep track of variants specifically, in the future it makes sense to make one that can be used for both purposes
    """
    if not excluded_regions or chrom not in excluded_regions:
        return False
    
    if excluded_regions is None:
        return False
    
    intervals = excluded_regions[chrom]
    # finding the insertion point in the sorted excluded ranges of the chromosome of interest
    idx = bisect.bisect_right(intervals, (start, end))

    # check interval immediately before
    if idx > 0: #if it is not the first
        prev_start, prev_end = intervals[idx - 1]
        if prev_end > start:
            return True
        
    # check interval immediately after
    if idx < len(intervals): #if it is not the last
        next_start, next_end = intervals[idx]
        if next_start < end:
            return True 
        
    return False








'''
PICKS A CHROMOSOME RANDOMLY TAKING INTO ACCOUNT THE LENGTH: WEIGHTED RANDOM SELECTION
outputs the chromosome and the respective length
'''
def select_chr(chroms:list, lengths:list):

    # random.choices always returns a list, so I set the len = 1 and take the first element
    chrom = random.choices(chroms, weights=lengths, k=1)[0]
    chrom_idx = chroms.index(chrom)

    return chrom, lengths[chrom_idx]    

'''
SELECTS A RANDOM POSITION ALONG THE SPECIFIED CHROMSOME
optionally it can avoid choosing positions too close to the chr end
'''
def select_pos(length, buffer=1000):

    if buffer > length/10:
        raise ValueError(f'edge buffer :{buffer} is too big, try decreasing it')
    
    # random.randint is inclusive of both ends
    position = random.randint(1 + buffer, length - buffer)

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




"""
Checks for overlaps using Binary Search. O(log N).
Requires sv_positions[chrom] to be sorted by start position.
Now haplotype-aware: only checks the haplotypes where the variant is placed.
"""
def overlaps(chrom, start, end, genotype_str, sv_positions: dict):

    alleles = genotype_str.split('|')

    for hap_idx, allele in enumerate(alleles):
        if allele == "1": #check only haplotypes that carry the variant

            intervals = sv_positions[chrom][hap_idx]
        
            if not intervals: # base case, empty list == no overlaps 
                continue

            # 1) Binary Search to find the "Insertion Point" 
            # We ask: "If I were to insert this new variant into the list while keeping it sorted, at which index would it go?"
            idx = bisect.bisect_right(intervals, (start, end))

            # 4. Check the Left Neighbor == the variant immediately BEFORE the new one
            if idx > 0:
                prev_start, prev_end = intervals[idx - 1]
                # If the previous variant's END extends past our START, there's overlap.
                if prev_end >= start:
                    return True

            # 5. Check the Right Neighbor == the variant immediately AFTER the new one
            if idx < len(intervals):
                next_start, next_end = intervals[idx]
                # If the next variant's START begins before our END, we overlap.
                if next_start <= end:
                    return True

    # otherwise there's no overlap
    return False


"""
Generates a genotype string (e.g., '0|1|0') with PHASED variants based on ploidy and heterozygosity.
Ensures that at least one allele is '1'
"""
def generate_genotype(ploidy:int, heterozygosity:float) -> str:

    #randomly assign alleles based on heterozygosity prob
    alleles = [1 if random.random() < heterozygosity else 0 for copy in range(ploidy)]

    #if all alleles are 0, we force at least one to be a 1
    if sum(alleles) == 0:
        random_idx = random.randint(0,ploidy-1)
        alleles[random_idx] = 1
    
    #format the output as a VCF genotype string (ex. "0/1")
    return "|".join(map(str,alleles))


"""
Parses genotype and inserts SV coordinates into the tracker.
"""
def track_sv(sv_positions, chrom, start, end, genotype_str):
    alleles = genotype_str.split('|') # Or '|' depending on generate_genotype output
    for hap_idx, allele in enumerate(alleles):
        if allele == "1":
            import bisect
            bisect.insort(sv_positions[chrom][hap_idx], (start, end))



### BNDs HANDLING - HELPER FUNCTIONS FOR THE MAIN LOOP ###


def find_valid_tra_coords(chroms, lengths, sv_length, gt, sv_positions, excluded_regions):
    """
    Performs a single attempt to generate and validate TRA coordinates.
    Returns (chrom_src, pos_src, chrom_dst, pos_dst) if valid, else None.
    """
    chrom_src, len_src = select_chr(chroms, lengths)
    pos_src = select_pos(len_src)
    
    chrom_dst, len_dst = select_chr(chroms, lengths)
    pos_dst = select_pos(len_dst)
    
    # 1. Check boundaries
    if pos_src + sv_length > len_src or pos_dst + 1 > len_dst:
        return None
        
    # 2. Check overlaps
    if (overlaps(chrom_src, pos_src, pos_src + sv_length, gt, sv_positions) or
        overlaps(chrom_dst, pos_dst, pos_dst + 1, gt, sv_positions) or
        overlaps_excluded_region(chrom_src, pos_src, pos_src + sv_length, excluded_regions) or
        overlaps_excluded_region(chrom_dst, pos_dst, pos_dst + 1, excluded_regions)):
        return None
        
    return chrom_src, pos_src, chrom_dst, pos_dst



def fetch_ref_base(chrom, pos, ref_fasta):
    """
    Fetches the true biological reference base at a 1-indexed VCF position.
    """
    # pysam is 0-indexed, so we fetch from (pos - 1) to (pos)
    return ref_fasta.fetch(chrom, pos - 1, pos)


def fetch_ref_span(chrom, start_1idx, end_1idx, ref_fasta):
    """
    Fetches a reference sequence span for 1-indexed inclusive VCF coordinates.
    Used by non-symbolic VCF generation to build explicit REF/ALT strings.
    """
    return ref_fasta.fetch(chrom, start_1idx - 1, end_1idx).upper()


def generate_tra_bnds(svtype, chrom_src, pos_src, chrom_dst, pos_dst, length, event_id, gt, ref_fasta, is_reverse=False):
    """
    Builds the Breakend objects for TRA_CUT or TRA_COPY, using the true reference bases.
    
    Args:
        is_reverse (bool): If True, the translocated segment is pasted in reverse-complement
                           orientation. Swaps paste BND brackets from t[p[ / ]p]t (forward)
                           to t]p] / [p[t (reverse) per VCF 4.2 spec.
    
    Returns a list of Breakend objects.
    """
    # Fetch true bases using the helper
    ref_dst = fetch_ref_base(chrom_dst, pos_dst, ref_fasta)
    ref_src_start = fetch_ref_base(chrom_src, pos_src + 1, ref_fasta)
    ref_src_end = fetch_ref_base(chrom_src, pos_src + length, ref_fasta)
    
    bnds = []
    
    if not is_reverse:
        # FORWARD ORIENTATION (current behavior, unchanged)
        # 1. PASTE START
        bnds.append(VariantObjects.Breakend(
            chrom_dst, pos_dst, f"{event_id}.P1", gt, f"{event_id}.P2", event_id, 
            f"{ref_dst}[{chrom_src}:{pos_src + 1}[", ref_dst))
            
        bnds.append(VariantObjects.Breakend(
            chrom_src, pos_src + 1, f"{event_id}.P2", gt, f"{event_id}.P1", event_id, 
            f"]{chrom_dst}:{pos_dst}]{ref_src_start}", ref_src_start))
            
        # 2. PASTE END
        bnds.append(VariantObjects.Breakend(
            chrom_src, pos_src + length, f"{event_id}.P3", gt, f"{event_id}.P4", event_id, 
            f"{ref_src_end}[{chrom_dst}:{pos_dst + 1}[", ref_src_end))
            
        bnds.append(VariantObjects.Breakend(
            chrom_dst, pos_dst + 1, f"{event_id}.P4", gt, f"{event_id}.P3", event_id, 
            f"]{chrom_src}:{pos_src + length}]{ref_dst}", ref_dst))
    else:
        # REVERSE ORIENTATION
        # The translocated segment is pasted in reverse-complement.
        # Per VCF 4.2 spec: t[p[ / ]p]t → t]p] / [p[t
        # Source coordinates also swap: dst now joins to src_end first (reading right-to-left)
        
        # re-fetch ref for dst_pos+1 (needed for P3/P4 in reverse)
        ref_dst_plus1 = fetch_ref_base(chrom_dst, pos_dst + 1, ref_fasta)

        # 1. PASTE START — dst joins to src_end (reverse extending left)
        bnds.append(VariantObjects.Breakend(
            chrom_dst, pos_dst, f"{event_id}.P1", gt, f"{event_id}.P2", event_id,
            f"{ref_dst}]{chrom_src}:{pos_src + length}]", ref_dst))
        
        bnds.append(VariantObjects.Breakend(
            chrom_src, pos_src + length, f"{event_id}.P2", gt, f"{event_id}.P1", event_id,
            f"{ref_src_end}]{chrom_dst}:{pos_dst}]", ref_src_end))
        
        # 2. PASTE END — src_start joins back to dst (reverse extending right)
        bnds.append(VariantObjects.Breakend(
            chrom_dst, pos_dst + 1, f"{event_id}.P3", gt, f"{event_id}.P4", event_id,
            f"[{chrom_src}:{pos_src + 1}[{ref_dst_plus1}", ref_dst_plus1))
        
        bnds.append(VariantObjects.Breakend(
            chrom_src, pos_src + 1, f"{event_id}.P4", gt, f"{event_id}.P3", event_id,
            f"[{chrom_dst}:{pos_dst + 1}[{ref_src_start}", ref_src_start))
    
    if svtype == "TRA_CUT":
        # 3. HEAL THE SOURCE (orientation-independent — always forward)
        ref_src_heal_start = fetch_ref_base(chrom_src, pos_src, ref_fasta)
        ref_src_heal_end = fetch_ref_base(chrom_src, pos_src + length + 1, ref_fasta)
        
        bnds.append(VariantObjects.Breakend(
            chrom_src, pos_src, f"{event_id}.H1", gt, f"{event_id}.H2", event_id, 
            f"{ref_src_heal_start}[{chrom_src}:{pos_src + length + 1}[", ref_src_heal_start))
            
        bnds.append(VariantObjects.Breakend(
            chrom_src, pos_src + length + 1, f"{event_id}.H2", gt, f"{event_id}.H1", event_id, 
            f"]{chrom_src}:{pos_src}]{ref_src_heal_end}", ref_src_heal_end))
            
    return bnds
