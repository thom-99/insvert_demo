from . import VariantObjects
from . import utils_sim
from collections import defaultdict
import bisect
import random
import pysam


def run(config_path, fasta_path, output_file, seed=None, excluded_bed=None, non_symbolic=False):

    # setting up the seed for reproducibility
    if seed is not None:
        print(f"Setting global random seed to: {seed}")
        random.seed(seed)

    print(f"Parsing config: {config_path}")
    config_info = utils_sim.parse_config(config_path)
    ploidy = config_info['ploidy']
    heterozygosity = config_info['heterozygosity']
    fakedict = config_info['variants']

    print(f"Reading index from: {fasta_path}")
    chroms, lengths = utils_sim.read_fai(fasta_path)

    excluded_regions=None
    if excluded_bed is not None:
        print(f"Reading the bed file with excluded genomic regions")
        excluded_regions = utils_sim.parse_bed(excluded_bed)
        if excluded_regions:
            print(f"Loaded excluded genomic regions from {excluded_bed}")


    # {chrom: {haplotype_index: [(pos, end)]}}
    # the lambda is a necessity because deafultdict requires a function... it's a bit convoluted
    sv_positions = defaultdict(lambda: defaultdict(list)) 


    with open(output_file, 'w') as vcf:

        header = utils_sim.buildheader(chroms, lengths, fasta_path)
        print(f"Writing {output_file} VCF header")
        vcf.write(header)

        count = 1

        #start to build SVs programmatically
        print("Building SVs...")
        ref_fasta = pysam.FastaFile(fasta_path)

        for svtype in fakedict:
            if svtype == 'INS':
                for l in fakedict[svtype]['lengths']:

                    # collect arguments of the SV orbject
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                    count += 1
                    
                    INS = VariantObjects.Insertion(chrom, pos, l, id, gt, ref_base)

                    # if either the SV is placed out of chromsome bounds or overlapping with another SV
                    # then another chromsome and position are chosen for the SV
                    attempts = 0
                    while (INS.get_end() > chrom_length or 
                           utils_sim.overlaps(chrom, pos, INS.get_end(), gt, sv_positions) or 
                           utils_sim.overlaps_excluded_region(chrom, pos, INS.get_end(), excluded_regions)):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count-1} could not be placed after 3 attempts, skipping')
                            break    
                        print(f'{svtype} exceeds the chormsome boundaries or overlaps with another SV, fetching a new position')
                        pos = utils_sim.select_pos(chrom_length)
                        ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                        INS = VariantObjects.Insertion(chrom, pos, l, id, gt, ref_base)
                    
                    # log SV and format a VCF line
                    if attempts <= 3:
                        if non_symbolic:
                            from .utils_ins import generate_seq
                            ins_sequence = generate_seq(l, 0.5)
                            alt_seq = f"{ref_base}{ins_sequence}"
                            INS = VariantObjects.Insertion(chrom, pos, l, id, gt, ref_base, alt_seq=alt_seq)
                        
                        # track SVs
                        utils_sim.track_sv(sv_positions, chrom, pos, INS.get_end(), gt)
                        vcf.write(INS.format() + '\n')
        
            if svtype == 'DEL':
                for l in fakedict[svtype]['lengths']:   

                    # collect arguments of the SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                    count += 1 
                    
                    l = -l # deletions require negative lengths, previously they have been made positive for statistical fitting
                    DEL = VariantObjects.Deletion(chrom, pos, l, id, gt, ref_base)

                    attempts = 0
                    while (DEL.get_end() > chrom_length or 
                           utils_sim.overlaps(chrom, pos, DEL.get_end(), gt, sv_positions) or
                           utils_sim.overlaps_excluded_region(chrom, pos, DEL.get_end(), excluded_regions)):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count-1} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom_length)
                        ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                        DEL = VariantObjects.Deletion(chrom, pos, l, id, gt, ref_base)

                    if attempts <= 3:
                        if non_symbolic:
                            full_ref = utils_sim.fetch_ref_span(chrom, pos, pos + abs(l), ref_fasta)
                            alt_seq = full_ref[0]
                            DEL = VariantObjects.Deletion(chrom, pos, l, id, gt, ref_base, alt_seq=alt_seq, ref_seq=full_ref)

                        # track SV
                        utils_sim.track_sv(sv_positions, chrom, pos, DEL.get_end(), gt)
                        vcf.write(DEL.format() + '\n')

            if svtype == 'INV':
                for l in fakedict[svtype]['lengths']:  

                    # collect arguments of the SV object  
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                    count += 1 

                    INV = VariantObjects.Inversion(chrom, pos, l, id, gt, ref_base)

                    # loop to produce valid pos to allow END to be within chromsome bounds
                    attempts = 0
                    while (INV.get_end() > chrom_length or 
                           utils_sim.overlaps(chrom, pos, INV.get_end(), gt, sv_positions) or
                           utils_sim.overlaps_excluded_region(chrom, pos, INV.get_end(), excluded_regions)):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count-1} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom_length)
                        ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                        INV = VariantObjects.Inversion(chrom, pos, l, id, gt, ref_base)

                    if attempts <= 3:
                        if non_symbolic:
                            anchor = ref_base
                            forward_seq = utils_sim.fetch_ref_span(chrom, pos + 1, pos + l, ref_fasta)
                            ref_seq = anchor + forward_seq
                            alt_seq = anchor + utils_sim.reverse_complement(forward_seq)
                            INV = VariantObjects.Inversion(chrom, pos, l, id, gt, ref_base, alt_seq=alt_seq, ref_seq=ref_seq)

                        # track SV
                        utils_sim.track_sv(sv_positions, chrom, pos, INV.get_end(), gt)
                        vcf.write(INV.format() + '\n')

            if svtype == 'DUP':
                for l,cn in zip(fakedict[svtype]['lengths'], fakedict[svtype]['copy_numbers']):   

                    # collect arguments for SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                    count += 1 

                    DUP = VariantObjects.Duplication(chrom, pos, l, id, gt, ref_base, copy_number=cn)

                    # loop to produce valid pos to allow END to be within chromsome bounds
                    attempts = 0
                    while (DUP.get_end() > chrom_length or 
                           utils_sim.overlaps(chrom, pos, DUP.get_end(), gt, sv_positions) or
                           utils_sim.overlaps_excluded_region(chrom, pos, DUP.get_end(), excluded_regions)):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count-1} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom_length)
                        ref_base = utils_sim.fetch_ref_base(chrom, pos, ref_fasta)
                        DUP = VariantObjects.Duplication(chrom, pos, l, id, gt, ref_base, cn)
                    
                    if attempts <= 3:
                        if non_symbolic:
                            anchor = ref_base
                            unit_seq = utils_sim.fetch_ref_span(chrom, pos + 1, pos + l, ref_fasta)
                            ref_seq = anchor + unit_seq
                            alt_seq = anchor + (unit_seq * cn)
                            DUP = VariantObjects.Duplication(chrom, pos, l, id, gt, ref_base, cn, alt_seq=alt_seq, ref_seq=ref_seq)

                        # track SV
                        utils_sim.track_sv(sv_positions, chrom, pos, DUP.get_end(), gt)
                        vcf.write(DUP.format() + '\n')


            # TRA_COPY & TRA_CUT processing
            # part of the processing is shared, differences are in the number of BND lines and the position logging
            # TRA_COPY & TRA_CUT processing
            if svtype in ["TRA_CUT", "TRA_COPY"]:
                for l in fakedict[svtype]['lengths']:
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    event_id = f'inSVert.{svtype}.{count}'
                    count += 1
                    
                    attempts = 0
                    coords = None
                    
                    # 1. Coordinate Hunting
                    while coords is None:
                        attempts += 1
                        if attempts > 10:
                            print(f'{svtype} n: {count-1} could not be placed after 10 attempts, skipping')
                            break
                        
                        # Roll the dice once for valid coordinates
                        coords = utils_sim.find_valid_tra_coords(
                            chroms, lengths, l, gt, sv_positions, excluded_regions
                        )

                    # skip this iteration if we didn't find coordinates
                    if coords is None:
                        continue
                        
                    chrom_src, pos_src, chrom_dst, pos_dst = coords
                    
                    # 2. Track the positions
                    utils_sim.track_sv(sv_positions, chrom_dst, pos_dst, pos_dst + 1, gt)
                    if svtype == "TRA_CUT":
                        utils_sim.track_sv(sv_positions, chrom_src, pos_src, pos_src + l, gt)
                        
                    # 3. Generate and write BNDs
                    reverse_ratio = fakedict[svtype].get('reverse_ratio', 0.0)
                    is_reverse = random.random() < reverse_ratio
                    bnds = utils_sim.generate_tra_bnds(
                        svtype, chrom_src, pos_src, chrom_dst, pos_dst, 
                        l, event_id, gt, ref_fasta, is_reverse
                    )
                    
                    for bnd in bnds:
                        vcf.write(bnd.format() + '\n')


    ref_fasta.close()
    print(f"VCF simulated. Output written to {output_file}")










