'''
UTILITIES FOR THE INSERT MODULE OF INSVERT
'''


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
APPLIES A INSERTION TO A CHROMSOME 
'''
def apply_insertion(chrom_seq : bytearray, ins_seq : str, pos : int, offset : int):

    adjusted_pos = pos + offset 

    # convert insertion seq to a bytearray
    if isinstance(ins_seq, str):
        ins_seq = bytearray(ins_seq.encode('ascii'))
    else:
        raise TypeError('Error : the sequence of the INS is not a str')
    
    # perform insertion
    new_chrom_seq = chrom_seq[:adjusted_pos] + ins_seq + chrom_seq[adjusted_pos:]

    # modify the original bytearray
    chrom_seq.clear()
    chrom_seq.extend(new_chrom_seq)

    offset += len(ins_seq)
    return offset

'''
APPLIES A DELETION TO A CHROMSOME
the edit is done inplace [O(1) space] thaks to the bytearray datastruct
'''
def apply_deletion(chrom_seq : bytearray, length : int, pos : int, offset : int):

    if length > 0:
        raise ValueError('DEL SVLENGTHS need to be negative')

    adjusted_pos = pos + offset

    # perform deletion in place 
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


