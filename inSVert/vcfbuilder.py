import os
import sys
import numpy as np
import random
import pysam

'''
takes as input the dictionary resuliting from simulator.py
builds a functional VCF file from that dictionary
'''

def read_fai(fasta_path):
    '''
    general function to read a fai file or create it if not already present
    '''
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
            print(f'index file created successfully : {fai_file}', file=sys.sterr)
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


def select_chr(chroms, lengths):
    '''
    selects chromsome randomly, but takes into account the length : weighted selection
    '''
    weights = np.array(lengths)
    probabilities = weights / weights.sum()

    chrom = np.random.choice(chroms, p=probabilities)
    chrom_idx = chroms.index(chrom)

    return chrom, lengths[chrom_idx]

#chroms, lengths = read_fai('data/cerevisiae.fa')
#randomchr = select_chr(chroms, lengths)
#print(randomchr)

