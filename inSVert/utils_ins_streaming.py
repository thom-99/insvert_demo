'''
UTILITIES FOR THE INSERT MODULE OF INSVERT
(Streaming Version)
'''
import random
import pysam
from collections import defaultdict

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




######################################################
### TRASLOCATION VALIDATION AND HANDLING FUNCTIONS ###
######################################################



def is_valid_tra(event_id, adjacencies):
    """
    Rigorously validates the biological structure of a translocation event.
    Supports both Inter-chromosomal and Intra-chromosomal (transposition) events.
    
    Returns: (type, source_chrom, source_range, sink_chrom, sink_pos) or None
    """
    count = len(adjacencies)
    
    # -------------------------------------------------------------------------
    # CATEGORY 1: COPY & PASTE (4 VCF lines / 2 Adjacencies)
    # -------------------------------------------------------------------------
    if count == 2:
        # Step 1: Group positions by chromosome
        chrom_map = defaultdict(list)
        for r, m in adjacencies:
            chrom_map[r.chrom].append(r.pos)
            chrom_map[m.chrom].append(m.pos)
            
        # Identify Source vs Sink
        if len(chrom_map) == 2:
            # INTER-chromosomal logic (existing)
            source_chrom = next(c for c, p in chrom_map.items() if len(set(p)) > 1)
            sink_chrom = next(c for c in chrom_map if c != source_chrom)
            source_pos = sorted(list(set(chrom_map[source_chrom])))
            sink_pos = min(chrom_map[sink_chrom])
        
        elif len(chrom_map) == 1:
            # INTRA-chromosomal logic
            # Heuristic: The two positions closest to each other (dist 1) are the Sink.
            # The other two positions define the Source segment.
            source_chrom = sink_chrom = list(chrom_map.keys())[0]
            all_pos = sorted(chrom_map[source_chrom])
            
            # Find which pair in the sorted list are adjacent (the sink join point)
            found_sink = False
            for i in range(len(all_pos) - 1):
                if all_pos[i+1] - all_pos[i] == 1:
                    sink_pos = all_pos[i]
                    # The other two positions are the source
                    source_pos = [p for p in all_pos if p != all_pos[i] and p != all_pos[i+1]]
                    found_sink = True
                    break
            if not found_sink: return None
        else:
            return None
        
        # Rigorous Logic: Source must be > 1bp
        if source_pos[1] - source_pos[0] <= 1: 
            return None 
            
        return ("TRA_COPY", source_chrom, (source_pos[0], source_pos[1]), 
                sink_chrom, sink_pos)

    # -------------------------------------------------------------------------
    # CATEGORY 2: CUT & PASTE (6 VCF lines / 3 Adjacencies)
    # -------------------------------------------------------------------------
    elif count == 3:
        # Step 1: Identify the 'HEAL' adjacency
        # we can consider the HEAL adj as the one that encompasses inside another adj, the COPY/SOURCE one

        heal_adj = None
        all_bnds = []
        for r, m in adjacencies:
            all_bnds.extend([(r.chrom, r.pos), (m.chrom, m.pos)])

        for r, m in adjacencies:
            # a heal adj must be interchromosomal
            if r.chrom == m.chrom:
                # find all the other breakpoints present on this same chrom
                other_bnds_on_this_chrom = [pos for chrom, pos in all_bnds if chrom == r.chrom and pos != r.pos and pos != m.pos]
                start, end = min(r.pos, m.pos), max(r.pos, m.pos)
                
                # TOPOLOGY CHECK !!!
                # if this adjacency spans (start to end) any other breakpoints it it the HEAL
                # the 'len==0' conditions covers the Inter-chromosomal case in which only the HEAL adj exists on the source chrom
                if any(start < p < end for p in other_bnds_on_this_chrom) or len(other_bnds_on_this_chrom) == 0:
                    heal_adj = (r, m)
                    break
        
        if not heal_adj: return None # if no HEAL adjacency is found, the TRA cannot be characterized
        
        source_chrom = heal_adj[0].chrom
        source_pos = sorted([heal_adj[0].pos, heal_adj[1].pos])
        
        # Rigorous Logic: Source must be > 1bp
        if source_pos[1] - source_pos[0] - 1 <= 1: 
            return None 

        # Step 2: Find the Sink from 'PASTE' adjacencies
        # the remaining 2 adj are the 'PASTE' junctions, they connect the segment ends to a sinlge insertion point
        paste_adjs = [a for a in adjacencies if a[0].id != heal_adj[0].id and a[0].id != heal_adj[1].id]
        
        
        sink_chrom = sink_pos = None
        for r, m in paste_adjs:
            for bnd in [r, m]:
                # the sink is the side of the junction that is NOT the segment
                # it is identifiable either by being on a different chrom OR by being outside of the range of the heal
                if bnd.chrom != source_chrom or (bnd.pos < source_pos[0] or bnd.pos > source_pos[1]):
                    sink_chrom, sink_pos = bnd.chrom, bnd.pos
                    break
            if sink_chrom: break #stop once the destination is found

        if not sink_chrom: return None
        
        return ("TRA_CUT", source_chrom, (source_pos[0], source_pos[1]), 
                sink_chrom, sink_pos)

    return None





def prefetch_translocations(vcf_path, ref_path):
    """
    ===========================================================================
    Scans the VCF for BND (breakend) records to reconstruct biological 
    translocation events. It pre-fetches translocated sequences and maps 
    them to specific genomic "triggers."

    Args:
        vcf_path (str): Path to the input VCF containing structural variants.
        ref_path (str): Path to the reference FASTA genome.

    Returns:
        dict: A nested dictionary (Action Map) structured as:
            {
                "deletions": { 
                    chrom (str): { pos (int): (length (int), event_id (str)) }
                },
                "insertions": { 
                    chrom (str): { pos (int): (sequence (str), event_id (str)) }
                }
            }
        - 'deletions' trigger a reference skip (the "cut").
        - 'insertions' trigger a sequence paste (the "paste").
    ===========================================================================
    """
    vcf = pysam.VariantFile(vcf_path)
    ref = pysam.FastaFile(ref_path)
    events = defaultdict(list)
    for var in vcf:
        if var.info.get("SVTYPE") == "BND":
            events[var.info.get("EVENT")].append(var)

    # Actions keyed by [chrom][pos] = (type, data, event_id)
    tra_map = {"deletions": defaultdict(dict), "insertions": defaultdict(dict)}

    for event_id, records in events.items():
        # Helper to pair BNDs by MATEID
        adjacencies = []
        seen = set()
        for r in records:
            if r.id not in seen:
                mate = next((m for m in records if m.id == r.info.get("MATEID")[0]), None)
                if mate: 
                    adjacencies.append((r, mate))
                    seen.update([r.id, mate.id])
        
        # check if the traslocation is valid and categorize it
        res = is_valid_tra(event_id, adjacencies)
        if not res: continue
        
        tra_type, src_chr, (s_start, s_end), snk_chr, snk_pos = res
        sequence = ref.fetch(src_chr, s_start, s_end - 1)

        if tra_type == "TRA_CUT":
            tra_map["deletions"][src_chr][s_start] = (len(sequence), event_id)
        
        tra_map["insertions"][snk_chr][snk_pos] = (sequence, event_id)

    vcf.close(); ref.close()
    return tra_map