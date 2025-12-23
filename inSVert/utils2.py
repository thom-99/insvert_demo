'''
UTILITIES FOR THE INSERT MODULE OF INSVERT
'''
import yaml
import random

'''
PARSES A FASTA FILE AND RETURNS A DICTIONARY WITH THE CHR : SEQUENCES
sequences are stored as bytearrays, also the offset is initialized as 0 for each chr
{'chr' : {'sequence : 'AATT..' , 'offset' : 0}}
'''
def parse_fasta(fasta_path:str):

    genome = {}
    current_chrom = None
    current_seq = []

    with open(fasta_path, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            
            # skip empy lines
            if not line:
                continue
            
            # process header
            if line.startswith('>'):
                # if it exists, save previous chr
                # concatenate a list of strings with nothing ('') 
                if current_chrom is not None:
                    genome[current_chrom] = {'sequence' : bytearray(''.join(current_seq).encode('ascii')),
                                             'offset' : 0}
                    
                # else start new chr
                current_chrom = line[1:].split()[0]
                current_seq = []
            
            else:
                # accumulate sequence
                current_seq.append(line.upper())

        # last chr processing
        if current_chrom is not None:
            genome[current_chrom] = {'sequence' : bytearray(''.join(current_seq).encode('ascii')),
                                             'offset' : 0}
        
    return genome


'''
RETRIEVES GLOBAL GC CONTENT FROM THE CONFIG.YAML 
defaults to 0.41 (Human)
'''

def get_GC(config_path:str) -> float:

    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    gc = config.get('gc_content')

    if gc is None:
        print('GC content non specified, defaulting to 0.41 (Human)')
        gc = 0.41
    
    if gc <= 0 or gc >= 1.0:
        raise ValueError(f"GC content must be between 0.0 and 1.0, got {gc} instead")
    
    return gc

'''
GENERATES A RANDOM DNA SEQUENCE GIVEN A LENGTH AND A GC CONTENT
returns a bytearray
'''
def generate_seq(length:int, gc_content:float) -> bytearray:

    if length < 0:
        raise ValueError("Length must be bigger than zero")
    
    prob_GC = gc_content / 2
    prob_AT = (1-gc_content) / 2

    bases = [b'A', b'C', b'G', b'T'] # we use bytes so we avoid encoding later 
    weights = [prob_AT, prob_GC, prob_GC, prob_AT]

    seq_list = random.choices(bases, weights, k=length)

    return bytearray(b''.join(seq_list))


'''
APPLIES A INSERTION TO A CHROMSOME 
'''
def apply_insertion(chrom_seq : bytearray, ins_seq : str, pos : int, offset : int):

    adjusted_pos = pos + offset 

    # convert insertion seq to a bytearray
    if isinstance(ins_seq, str):
        ins_seq = bytearray(ins_seq.encode('ascii'))
    elif isinstance(ins_seq, (bytearray, bytes)):
        pass
    else:
        raise TypeError(f'Error: INS sequence must be str or bytearray, got {type(ins_seq)}')
    
    # perform insertion in-place
    chrom_seq[adjusted_pos:adjusted_pos] = ins_seq

    offset += len(ins_seq)
    return offset

'''
APPLIES A TANDEM DUPLICATION BY CALLING APPLY_INSERTION
Handles the logic of extracting the sequence and calculating the correct tandem position.
'''
def apply_duplication(chrom_seq: bytearray, pos:int, length:int, copy_number:int, offset):

    if copy_number<= 1: # sanity check: cn=1 means == reference
        return offset 
    
    # extract sequence to be duplicated from genome
    current_start_idx = pos + offset
    current_end_idx = current_start_idx + length

    dup_seq = chrom_seq[current_start_idx:current_end_idx]
    repeats = copy_number - 1  # if CN=3, we have 1 original copy and 2 new copies to insert
    ins_seq = dup_seq * repeats 

    # calculate where to insert them
    # tandem dups hppen just AFTER the original unit
    ins_pos = pos + length

    return apply_insertion(chrom_seq, ins_seq, ins_pos, offset)


'''
APPLIES A DELETION TO A CHROMSOME
the edit is done inplace [O(1) space] thaks to the bytearray datastruct
'''
def apply_deletion(chrom_seq : bytearray, length : int, pos : int, offset : int):

    if length > 0:
        raise ValueError('DEL SV_LENGTHS need to be negative')

    adjusted_pos = pos + offset

    # perform deletion in-place 
    del chrom_seq[adjusted_pos : adjusted_pos + abs(length)]

    offset += length
    return offset


'''
reverse complement a bytearray, helper for apply_inversion
'''
def reverse_complement(sequence):

    complement = {
        ord('A'): ord('T'),
        ord('T'): ord('A'),
        ord('G'): ord('C'),
        ord('C'): ord('G'),
        ord('N'): ord('N'),
        ord('a'): ord('t'),
        ord('t'): ord('a'),
        ord('g'): ord('c'),
        ord('c'): ord('g'),
        ord('n'): ord('n')
    }
    
    # reverse the sequence and complement each base
    rev_comp = bytearray(len(sequence))
    for i, base in enumerate(reversed(sequence)):
        rev_comp[i] = complement.get(base, base)
    
    return rev_comp


'''
APPLIES AN INVERSION TO A CHROMOSOME
the edit is done inplace [O(1) space] by replacing the region with its reverse complement
'''
def apply_inversion(chrom_seq : bytearray, pos : int, end : int, offset : int):
    
    adjusted_pos = pos + offset
    adjusted_end = end + offset
    
    # get the region to invert
    region = chrom_seq[adjusted_pos:adjusted_end]
    
    # replace with reverse complement in place
    chrom_seq[adjusted_pos:adjusted_end] = reverse_complement(region)
    
    return offset


