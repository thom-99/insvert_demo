'''
UTILITIES FOR THE INSERT MODULE OF INSVERT
'''


'''
PARSES A FASTA FILE AND RETURNS A DICTIONARY WITH THE CHR : SEQUENCES
sequences are stored as bytearrays, also the offset is initialized as 0 for each chr
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


genome = parse_fasta('data/cerevisiae.fa')
print(genome)