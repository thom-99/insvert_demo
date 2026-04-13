from . import VariantObjects
from . import utils_sim
from collections import defaultdict
import bisect
import random


def run(config_path, fasta_path, output_file, seed=None):

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
                for l in fakedict[svtype]['lengths']:

                    # collect arguments of the SV orbject
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
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
                        pos = utils_sim.select_pos(chrom_length)
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
                for l in fakedict[svtype]['lengths']:   

                    # collect arguments of the SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
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
                        pos = utils_sim.select_pos(chrom_length)
                        DEL = VariantObjects.Deletion(chrom, pos, l, id, gt)

                    if attempts <= 3:
                        alleles = gt.split('/')
                        for hap_idx, allele in enumerate(alleles):
                            if allele == "1":
                                bisect.insort(sv_positions[chrom][hap_idx], (pos, DEL.get_end()))
                        vcf.write(DEL.format() + '\n')

            if svtype == 'INV':
                for l in fakedict[svtype]['lengths']:  

                    # collect arguments of the SV object  
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
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
                        pos = utils_sim.select_pos(chrom_length)
                        INV = VariantObjects.Inversion(chrom, pos, l, id, gt)

                    if attempts <= 3:
                        alleles = gt.split('/')
                        for hap_idx, allele in enumerate(alleles):
                            if allele == "1":
                                bisect.insort(sv_positions[chrom][hap_idx], (pos, INV.get_end()))
                        vcf.write(INV.format() + '\n')

            if svtype == 'DUP':
                for l,cn in zip(fakedict[svtype]['lengths'], fakedict[svtype]['copy_numbers']):   

                    # collect arguments for SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom_length)
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
                        pos = utils_sim.select_pos(chrom_length)
                        DUP = VariantObjects.Duplication(chrom, pos, l, id, gt, cn)
                    
                    if attempts <= 3:
                        alleles = gt.split('/')
                        for hap_idx, allele in enumerate(alleles):
                            if allele == "1":
                                bisect.insort(sv_positions[chrom][hap_idx], (pos, DUP.get_end()))
                        vcf.write(DUP.format() + '\n')


            # TRA_COPY & TRA_CUT processing
            # part of the processing is shared, differences are in the number of BND lines and the position logging
            if svtype == "TRA_CUT" or svtype=="TRA_COPY":
                for l in fakedict[svtype]['lengths']:

                    # collecting the arguments to create the objects 
                    chrom_src, len_src = utils_sim.select_chr(chroms, lengths)
                    pos_src = utils_sim.select_pos(len_src)

                    chrom_dst, len_dst = utils_sim.select_chr(chroms, lengths)
                    pos_dst = utils_sim.select_pos(len_dst)

                    event_id = f'inSVert.{svtype}.{count}'
                    gt = utils_sim.generate_genotype(ploidy, heterozygosity)
                    count += 1

                    # overlap checks 
                    attempts = 0 
                    while (pos_src + l > len_src or 
                           pos_dst + 1 > len_dst or 
                           utils_sim.overlaps(chrom_src, pos_src, pos_src + l, gt, sv_positions) or
                           utils_sim.overlaps(chrom_dst, pos_dst, pos_dst + 1, gt, sv_positions)):
                        
                        attempts += 1
                        if attempts > 10:
                            print(f'{svtype} n: {count} could not be placed after 10 attempts, skipping')
                            break
                        print(f'{svtype} exceeds boundaries or overlaps, fetching new positions')
                        # re-selecting both source and destination
                        chrom_src, len_src = utils_sim.select_chr(chroms, lengths)
                        pos_src = utils_sim.select_pos(len_src)
                        chrom_dst, len_dst = utils_sim.select_chr(chroms, lengths)
                        pos_dst = utils_sim.select_pos(len_dst)

                    if attempts <= 10:
                        # instantiate paste BND IDs
                        id_p1, id_p2 = f"{event_id}.P1", f"{event_id}.P2"
                        id_p3, id_p4 = f"{event_id}.P3", f"{event_id}.P4"

                        # adjacency 1 : PASTE START
                        # joins destination falnk to the beginning of the moved segment
                        bnd_p1 = VariantObjects.Breakend(chrom_dst, pos_dst, id_p1, gt, id_p2, event_id, f"N[{chrom_src}:{pos_src + 1}[")
                        bnd_p2 = VariantObjects.Breakend(chrom_src, pos_src + 1, id_p2, gt, id_p1, event_id, f"]{chrom_dst}:{pos_dst}]N")
                        # adjacency 2 : PASTE END
                        # joins the end of the moved segment to the destination right flank
                        bnd_p3 = VariantObjects.Breakend(chrom_src, pos_src + l, id_p3, gt, id_p4, event_id, f"N[{chrom_dst}:{pos_dst + 1}[")
                        bnd_p4 = VariantObjects.Breakend(chrom_dst, pos_dst + 1, id_p4, gt, id_p3, event_id, f"]{chrom_src}:{pos_src+l}]N")
                        
                        if svtype=="TRA_CUT":

                            # adjacency 3 : HEAL THE SOURCE
                            # joins flank before segment to flank after segment, healing the rupture
                            id_h1, id_h2 = f"{event_id}.H1", f"{event_id}.H2"
                            bnd_h1 = VariantObjects.Breakend(chrom_src, pos_src, id_h1, gt, id_h2, event_id, f"N[{chrom_src}:{pos_src + l + 1}[")
                            bnd_h2 = VariantObjects.Breakend(chrom_src, pos_src + l + 1, id_h2, gt, id_h1, event_id, f"]{chrom_src}:{pos_src}]N")                        


                            # update overlap tracker with both positions
                            alleles = gt.split('/')
                            for hap_idx, allele in enumerate(alleles):
                                if allele == "1":
                                    # Track the source as a deleted interval
                                    bisect.insort(sv_positions[chrom_src][hap_idx], (pos_src, pos_src + l))
                                    # Track the destination as a point insertion
                                    bisect.insort(sv_positions[chrom_dst][hap_idx], (pos_dst, pos_dst + 1))
                            
                            # write 6 BNDs lines (4 paste + 2 heal)
                            vcf.write(bnd_p1.format() + '\n')
                            vcf.write(bnd_p2.format() + '\n')
                            vcf.write(bnd_p3.format() + '\n')
                            vcf.write(bnd_p4.format() + '\n')
                            vcf.write(bnd_h1.format() + '\n')
                            vcf.write(bnd_h2.format() + '\n')
                        
                        
                        elif svtype=="TRA_COPY":
                            
                            alleles = gt.split('/')
                            for hap_idx, allele in enumerate(alleles):
                                if allele == "1":
                                    # track only the point insertion
                                    bisect.insort(sv_positions[chrom_dst][hap_idx], (pos_dst, pos_dst + 1))

                            # write 4 BND lines (paste)
                            vcf.write(bnd_p1.format() + '\n')
                            vcf.write(bnd_p2.format() + '\n')
                            vcf.write(bnd_p3.format() + '\n')
                            vcf.write(bnd_p4.format() + '\n')




    print(f"VCF simulated. Output written to {output_file}")










