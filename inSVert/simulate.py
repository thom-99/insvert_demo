from . import VariantObjects
from . import utils_sim
from collections import defaultdict


def run(config_path, fasta_path, output_file):

    print(f"Parsing config: {config_path}")
    fakedict = utils_sim.parse_config(config_path)
    print(f"Reading index from: {fasta_path}")
    chroms, lengths = utils_sim.read_fai(fasta_path)

    sv_positions = defaultdict(list) # {chrom : [(pos, end),(pos2, end2)]}


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
                    count += 1

                    INS = VariantObjects.Insertion(chrom, pos, l, id)

                    # if either the SV is placed out of chromsome bounds or overlapping with another SV
                    # then another chromsome and position are chosen for the SV
                    attempts = 0
                    while INS.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, INS.get_end(), sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break    
                        print(f'{svtype} exceeds the chormsome boundaries or overlaps with another SV, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        INS = VariantObjects.Insertion(chrom, pos, l, id)
                    
                    # log SV and format a VCF line
                    if attempts <= 3:
                        sv_positions[chrom].append((pos, INS.get_end()))
                        vcf.write(INS.format() + '\n')
        
            if svtype == 'DEL':
                for l in fakedict['DEL']['lengths']:   

                    # collect arguments of the SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom, chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    count += 1 
                    
                    l = -l # deletions require negative lengths, previously they have been made positive for statistical fitting
                    DEL = VariantObjects.Deletion(chrom, pos, l, id)

                    attempts = 0
                    while DEL.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, DEL.get_end(), sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        DEL = VariantObjects.Deletion(chrom, pos, l, id)

                    if attempts <= 3:
                        sv_positions[chrom].append((pos, DEL.get_end()))
                        vcf.write(DEL.format() + '\n')

            if svtype == 'INV':
                for l in fakedict['INV']['lengths']:  

                    # collect arguments of the SV object  
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom, chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    count += 1 

                    INV = VariantObjects.Inversion(chrom, pos, l, id)

                    # loop to produce valid pos to allow END to be within chromsome bounds
                    attempts = 0
                    while INV.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, INV.get_end(), sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        INV = VariantObjects.Inversion(chrom, pos, l, id)

                    if attempts <= 3:
                        sv_positions[chrom].append((pos, INV.get_end()))                
                        vcf.write(INV.format() + '\n')

            if svtype == 'DUP':
                for l,cn in zip(fakedict['DUP']['lengths'], fakedict['DUP']['copy_numbers']):   

                    # collect arguments for SV object 
                    chrom, chrom_length = utils_sim.select_chr(chroms, lengths)
                    pos = utils_sim.select_pos(chrom, chrom_length)
                    id = f'inSVert.{svtype}.{count}'
                    count += 1 

                    DUP = VariantObjects.Duplication(chrom, pos, l, id, copy_number=cn)

                    # loop to produce valid pos to allow END to be within chromsome bounds
                    attempts = 0
                    while DUP.get_end() > chrom_length or utils_sim.overlaps(chrom, pos, DUP.get_end(), sv_positions):
                        attempts += 1
                        if attempts > 3:
                            print(f'{svtype} n: {count} could not be placed after 3 attempts, skipping')
                            break                     
                        print(f'{svtype} exceeds the chormsome boundaries, fetching a new position')
                        pos = utils_sim.select_pos(chrom, chrom_length)
                        DUP = VariantObjects.Duplication(chrom, pos, l, id, cn)
                    
                    if attempts <= 3:
                        sv_positions[chrom].append((pos, DUP.get_end()))
                        vcf.write(DUP.format() + '\n')


    print(f"VCF simulated. Output written to {output_file}")











# experimentals 

import pysam
import sys
from . import utils_ins

'''
inSVert module that performs the insertion of a VCF file into a fasta reference
using a robust, low-memory STREAMING approach.
'''

def write_wrapped(file_handle, sequence, line_width=60):
    """Writes sequence to file with FASTA line wrapping."""
    for i in range(0, len(sequence), line_width):
        file_handle.write(sequence[i:i+line_width] + '\n')

def run2(gc_content, ref_fasta, vcf_file, output_fasta):

    print(f'Streaming variants from {vcf_file}...')
    
    # 1. Open Reference and VCF
    ref = pysam.FastaFile(ref_fasta)
    vcf = pysam.VariantFile(vcf_file)
    
    with open(output_fasta, 'w') as out_f:
        
        # 2. Iterate over each chromosome in the reference
        for chrom in ref.references:
            print(f"Processing {chrom}...", end="\r")
            
            # Write chromosome header
            out_f.write(f">{chrom}\n")
            
            # Fetch variants for ONLY this chromosome
            try:
                # We assume VCF is sorted (handled by cli.py validation)
                chrom_variants = list(vcf.fetch(chrom))
            except ValueError:
                chrom_variants = []

            # 3. Stream the sequence
            ref_pos = 0 # Current position in the reference (0-based)
            chrom_len = ref.get_reference_length(chrom)
            
            # Buffer for FASTA line wrapping
            # We accumulate sequence here and flush to file when len > 60
            line_buffer = "" 
            
            for var in chrom_variants:
                start = var.pos - 1 # VCF is 1-based, Python is 0-based
                
                # --- OVERLAP CHECK ---
                if start < ref_pos:
                    # If this variant starts before our current position, it overlaps 
                    # with a previous variant we just processed.
                    print(f"\n[Warning] Skipping overlapping variant at {chrom}:{var.pos}. "
                          f"Current Ref Position: {ref_pos}")
                    continue

                # --- VALIDATE SVLEN ---
                svlen = var.info.get("SVLEN")
                # Handle tuple return (some VCFs return (len,)) or None
                if isinstance(svlen, tuple): svlen = svlen[0]
                
                # Fallback if SVLEN is missing: infer from coordinates (less accurate for INS)
                if svlen is None:
                    svlen = var.stop - var.start
                
                svtype = var.info.get("SVTYPE")
                
                # A. Write Reference Sequence up to the variant start
                if start > ref_pos:
                    # Fetch chunk from reference
                    chunk = ref.fetch(chrom, ref_pos, start)
                    
                    # Instead of writing directly, add to buffer and wrap
                    line_buffer += chunk
                    while len(line_buffer) >= 60:
                        out_f.write(line_buffer[:60] + '\n')
                        line_buffer = line_buffer[60:]
                        
                    ref_pos = start
                
                # B. Handle Variants
                # We calculate the sequence to ADD to the buffer
                alt_seq_fragment = ""
                
                if svtype == 'INS':
                    # Generate random sequence
                    length = abs(svlen)
                    alt_seq_fragment = utils_ins.generate_seq(length, gc_content).decode('ascii')
                    # Ref_pos does NOT move for simple insertions (start == end in VCF usually, 
                    # but logic holds: we insert text, but consume 0 reference bases unless it's a replacement)
                    
                elif svtype == 'DEL':
                    # We write nothing, just skip reference
                    ref_pos += abs(svlen)
                    
                elif svtype == 'DUP':
                    dup_len = abs(svlen)
                    # Fetch the sequence to be duplicated (the "original")
                    dup_src = ref.fetch(chrom, start, start + dup_len)
                    
                    # Determine Copy Number
                    # Safely get first sample or default to 2
                    cn = 2 
                    if var.samples:
                        # Grab first sample name dynamically
                        sample_name = list(var.samples.keys())[0]
                        cn_val = var.samples[sample_name].get('CN')
                        if cn_val is not None:
                            cn = cn_val
                    
                    # Logic: We write the sequence CN times total.
                    # Then we ADVANCE ref_pos past the original region so it isn't written again by the loop.
                    # Memory Safety: Write in a loop, don't build huge string
                    for _ in range(cn):
                        line_buffer += dup_src
                        while len(line_buffer) >= 60:
                            out_f.write(line_buffer[:60] + '\n')
                            line_buffer = line_buffer[60:]
                            
                    ref_pos += dup_len # Skip the original on reference read

                elif svtype == 'INV':
                    inv_len = abs(svlen)
                    end = start + inv_len
                    
                    # Fetch sequence
                    seq_to_invert = ref.fetch(chrom, start, end)
                    
                    # Reverse Complement
                    seq_bytes = bytearray(seq_to_invert.encode('ascii'))
                    inv_bytes = utils_ins.reverse_complement(seq_bytes)
                    alt_seq_fragment = inv_bytes.decode('ascii')
                    
                    # Move reference pointer to skip the original (non-inverted) sequence
                    ref_pos = end

                # If we generated a fragment (INS/INV), buffer it
                if alt_seq_fragment:
                    line_buffer += alt_seq_fragment
                    while len(line_buffer) >= 60:
                        out_f.write(line_buffer[:60] + '\n')
                        line_buffer = line_buffer[60:]

            # 4. Write any remaining reference sequence after the last variant
            if ref_pos < chrom_len:
                chunk = ref.fetch(chrom, ref_pos, chrom_len)
                line_buffer += chunk
                while len(line_buffer) >= 60:
                    out_f.write(line_buffer[:60] + '\n')
                    line_buffer = line_buffer[60:]
            
            # Flush remaining buffer
            if line_buffer:
                out_f.write(line_buffer + '\n')
            
    vcf.close()
    ref.close()
    print(f'\nDone! Output written to {output_fasta}')