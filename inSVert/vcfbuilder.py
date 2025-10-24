import VariantObjects
import utils
import sys
from collections import defaultdict

if len(sys.argv) != 4:
    print("Usage: python script.py <vcf_path> <fasta_path> <output_file>")
    sys.exit(1)

vcf_path = sys.argv[1]
fasta_path = sys.argv[2]
output_file = sys.argv[3]

realdict = utils.parse_vcf(vcf_path)
fakedict = utils.simdict(realdict)

chroms, lengths = utils.read_fai(fasta_path)

sv_positions = defaultdict(list) # {chrom : [(pos, end),(pos2, end2)]}


# IT MIGHT BE THAT THE CONTIG IS TOO SHORT, IMPLEMENT A COUNTER : IF IT FAILS 3 TIMES TO PLACE THE SV IN THE CHR, CHOOSE ANOTHER CHR 
# IMPLEMENT THIS IN A FINCTION IN UTILS.PY

with open(output_file, 'w') as vcf:
 
    header = utils.buildheader(chroms, lengths, fasta_path)
    vcf.write(header)

    count = 1

    #start to build SVs programmatically
    for svtype in fakedict:
        if svtype == 'INS':
            for l in fakedict['INS']['lengths']:

                # collect arguments of the SV orbject
                chrom, chrom_length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, chrom_length)
                id = f'inSVert.{svtype}.{count}'
                count += 1

                INS = VariantObjects.Insertion(chrom, pos, l, id)

                # if either the SV is placed out of chromsome bounds or overlapping with another SV
                # then another chromsome and position are chosen for the SV
                attempts = 0
                while INS.get_end() > chrom_length or utils.overlaps(chrom, pos, INS.get_end(), sv_positions):
                    attempts += 1
                    if attempts > 3:
                        print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                        break    
                    print(f'{svtype} exceeds the chormsome boundaries or overlaps with another SV, fetching a new position')
                    pos = utils.select_pos(chrom, chrom_length)
                    INS = VariantObjects.Insertion(chrom, pos, l, id)
                
                # log SV and format a VCF line
                if attempts <= 3:
                    sv_positions[chrom].append((pos, INS.get_end()))
                    vcf.write(INS.format() + '\n')
    
        if svtype == 'DEL':
            for l in fakedict['DEL']['lengths']:   

                # collect arguments of the SV object 
                chrom, chrom_length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, chrom_length)
                id = f'inSVert.{svtype}.{count}'
                count += 1 
                
                l = -l # deletions require negative lengths, previously they have been made positive for statistical fitting
                DEL = VariantObjects.Deletion(chrom, pos, l, id)

                attempts = 0
                while DEL.get_end() > chrom_length or utils.overlaps(chrom, pos, DEL.get_end(), sv_positions):
                    attempts += 1
                    if attempts > 3:
                        print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                        break                     
                    print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                    pos = utils.select_pos(chrom, chrom_length)
                    DEL = VariantObjects.Deletion(chrom, pos, l, id)

                if attempts <= 3:
                    sv_positions[chrom].append((pos, DEL.get_end()))
                    vcf.write(DEL.format() + '\n')

        if svtype == 'INV':
            for l in fakedict['INV']['lengths']:  

                # collect arguments of the SV object  
                chrom, chrom_length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, chrom_length)
                id = f'inSVert.{svtype}.{count}'
                count += 1 

                INV = VariantObjects.Inversion(chrom, pos, l, id)

                # loop to produce valid pos to allow END to be within chromsome bounds
                attempts = 0
                while INV.get_end() > chrom_length or utils.overlaps(chrom, pos, INV.get_end(), sv_positions):
                    attempts += 1
                    if attempts > 3:
                        print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                        break                     
                    print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                    pos = utils.select_pos(chrom, chrom_length)
                    INV = VariantObjects.Inversion(chrom, pos, l, id)

                if attempts <= 3:
                    sv_positions[chrom].append((pos, INV.get_end()))                
                    vcf.write(INV.format() + '\n')

        if svtype == 'DUP':
            for l in fakedict['DUP']['lengths']:   

                # collect arguments for SV object 
                chrom, chrom_length = utils.select_chr(chroms, lengths)
                pos = utils.select_pos(chrom, chrom_length)
                id = f'inSVert.{svtype}.{count}'
                count += 1 

                DUP = VariantObjects.Duplication(chrom, pos, l, id)

                # loop to produce valid pos to allow END to be within chromsome bounds
                attempts = 0
                while DUP.get_end() > chrom_length or utils.overlaps(chrom, pos, DUP.get_end(), sv_positions):
                    attempts += 1
                    if attempts > 3:
                        print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                        break                     
                    print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                    pos = utils.select_pos(chrom, chrom_length)
                    DUP = VariantObjects.Duplication(chrom, pos, l, id)
                
                if attempts <= 3:
                    sv_positions[chrom].append((pos, DUP.get_end()))
                    vcf.write(DUP.format() + '\n')


