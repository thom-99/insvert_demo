from . import VariantObjects
from . import utils_sim
from collections import defaultdict
import bisect


def run(config_path, fasta_path, output_file):

    print(f"Parsing config: {config_path}")
    config_info = utils_sim.parse_config(config_path)
    ploidy = config_info['ploidy']
    heterozygosity = config_info['heterozygosity']
    fakedict = config_info['variants']

    print(f"Reading index from: {fasta_path}")
    chroms, lengths = utils_sim.read_fai(fasta_path)

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
        for svtype in fakedict:
            if svtype == 'INS':
                for l in fakedict['INS']['lengths']:

                    # collect arguments of the SV orbject
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom, chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    count += 1

                    INS = VariantObjects.Insertion(chrom, pos, l, id, gt)

                    # if either the SV is placed out of chromsome bounds or overlapping with another SV
                    # then another chromsome and position are chosen for the SV
                    attempts = 0
                    while INS.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, INS.get_end(), gt, sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break    
                        print(f'{svtype} exceeds the chormsome boundaries or overlaps with another SV, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        INS = VariantObjects.Insertion(chrom, pos, l, id, gt)
                    
                    # log SV and format a VCF line
                    if attempts <= 3:

                        # Split genotype and insert into the correct haplotype lists
                        alleles = gt.split('/')
                        for hap_idx, allele in enumerate(alleles):
                            if allele == "1":
                                bisect.insort(sv_positions[chrom][hap_idx], (pos, INS.get_end()))
                        vcf.write(INS.format() + '\n')
        
            if svtype == 'DEL':
                for l in fakedict['DEL']['lengths']:   

                    # collect arguments of the SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom, chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    count += 1 
                    
                    l = -l # deletions require negative lengths, previously they have been made positive for statistical fitting
                    DEL = VariantObjects.Deletion(chrom, pos, l, id, gt)

                    attempts = 0
                    while DEL.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, DEL.get_end(), gt, sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        DEL = VariantObjects.Deletion(chrom, pos, l, id, gt)

                    if attempts <= 3:
                        alleles = gt.split('/')
                        for hap_idx, allele in enumerate(alleles):
                            if allele == "1":
                                bisect.insort(sv_positions[chrom][hap_idx], (pos, DEL.get_end()))
                        vcf.write(DEL.format() + '\n')

            if svtype == 'INV':
                for l in fakedict['INV']['lengths']:  

                    # collect arguments of the SV object  
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom, chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    count += 1 

                    INV = VariantObjects.Inversion(chrom, pos, l, id, gt)

                    # loop to produce valid pos to allow END to be within chromsome bounds
                    attempts = 0
                    while INV.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, INV.get_end(), gt, sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        INV = VariantObjects.Inversion(chrom, pos, l, id, gt)

                    if attempts <= 3:
                        alleles = gt.split('/')
                        for hap_idx, allele in enumerate(alleles):
                            if allele == "1":
                                bisect.insort(sv_positions[chrom][hap_idx], (pos, INV.get_end()))
                        vcf.write(INV.format() + '\n')

            if svtype == 'DUP':
                for l,cn in zip(fakedict['DUP']['lengths'], fakedict['DUP']['copy_numbers']):   

                    # collect arguments for SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom, chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    count += 1 

                    DUP = VariantObjects.Duplication(chrom, pos, l, id, gt, copy_number=cn)

                    # loop to produce valid pos to allow END to be within chromsome bounds
                    attempts = 0
                    while DUP.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, DUP.get_end(), gt, sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        DUP = VariantObjects.Duplication(chrom, pos, l, id, gt, cn)
                    
                    if attempts <= 3:
                        alleles = gt.split('/')
                        for hap_idx, allele in enumerate(alleles):
                            if allele == "1":
                                bisect.insort(sv_positions[chrom][hap_idx], (pos, DUP.get_end()))
                        vcf.write(DUP.format() + '\n')


    print(f"VCF simulated. Output written to {output_file}")










