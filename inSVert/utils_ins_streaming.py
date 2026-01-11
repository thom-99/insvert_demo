'''
UTILITIES FOR THE INSERT MODULE OF INSVERT
(Streaming Version)
'''
import random

def generate_seq(length:int, gc_content:float) -> str:
    """Generates a random DNA sequence as a string."""
    if length < 0:
        raise ValueError("Length must be bigger than zero")
    
    prob_GC = gc_content / 2
    prob_AT = (1-gc_content) / 2
    
    # Return string directly for streaming
    bases = ['A', 'C', 'G', 'T']
    weights = [prob_AT, prob_GC, prob_GC, prob_AT]
    return ''.join(random.choices(bases, weights, k=length))

def reverse_complement(sequence: str) -> str:
    """Returns reverse complement string."""
    complement = str.maketrans("ATGCNatgcn", "TACGNtacgn")
    return sequence.translate(complement)[::-1]

# --- STREAMING APPLICATION FUNCTIONS ---

def apply_insertion(out_buffer, ins_seq: str):
    """
    Applies an INSERTION by simply writing the sequence to the buffer.
    Does not consume reference position.
    """
    out_buffer.write(ins_seq)

def apply_deletion(ref_pos: int, length: int):
    """
    Applies a DELETION by advancing the reference position.
    Returns the new reference position.
    """
    # Simply skip the next 'length' bases
    return ref_pos + abs(length)

def apply_inversion(ref_file, chrom: str, start: int, length: int, out_buffer):
    """
    Applies an INVERSION: fetches ref sequence, reverse complements it, writes it.
    Returns the new reference position (end of inversion).
    """
    end = start + length
    # 1. Fetch original sequence
    seq = ref_file.fetch(chrom, start, end)
    # 2. Invert
    inv_seq = reverse_complement(seq)
    # 3. Write
    out_buffer.write(inv_seq)
    
    return end

def apply_duplication(ref_file, chrom: str, start: int, length: int, copy_number: int, out_buffer):
    """
    Applies a DUPLICATION: fetches sequence, multiplies it, writes it.
    Uses 'apply_insertion' logic implicitly (writing to buffer).
    Returns new ref position (skipping the original unit to avoid double-writing).
    """
    if copy_number <= 0:
        return start # Should be a deletion if CN=0, but strictly handling DUP here
        
    # 1. Fetch the sequence unit
    dup_unit = ref_file.fetch(chrom, start, start + length)
    
    # 2. Create the full sequence (Original + Copies)
    # We construct the payload similar to an insertion
    full_seq = dup_unit * copy_number
    
    # 3. Delegate to insertion logic
    apply_insertion(out_buffer, full_seq)
    
    # 4. Advance reference pointer past the original unit 
    # (because we just wrote it as part of full_seq)
    return start + length