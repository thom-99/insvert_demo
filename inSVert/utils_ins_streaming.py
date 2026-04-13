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
    PHASE 1: TRANSLOCATION PREFETCH & VALIDATION
    ===========================================================================
    Scans the VCF to reconstruct biological translocation events from BND lines.
    Builds an action map (deletions/insertions) for the streaming insertion loop.
    ===========================================================================
    """
    vcf = pysam.VariantFile(vcf_path)
    ref = pysam.FastaFile(ref_path)
    
    # Step 1: Group all BND records by their EVENT ID
    events = defaultdict(list)
    for var in vcf:
        if var.info.get("SVTYPE") == "BND":
            event_id = var.info.get("EVENT")
            if event_id:
                events[event_id].append(var)

    # Output structure for O(1) lookup during streaming
    tra_map = {
        "deletions": defaultdict(dict),  # {chrom: {event_id: length}}
        "insertions": defaultdict(dict)  # {chrom: {event_id: sequence}}
    }

    # Step 2: Categorize and Extract Sequences
    for event_id, records in events.items():
        # Pair records by MATEID to find novel adjacencies (junctions)
        adjacencies = []
        seen = set()
        for r in records:
            mate_id = r.info.get("MATEID")[0] if r.info.get("MATEID") else None
            if r.id not in seen and mate_id:
                mate = next((m for m in records if m.id == mate_id), None)
                if mate:
                    adjacencies.append((r, mate))
                    seen.update([r.id, mate.id])

        # Step 3: Run Rigorous Validation
        result = is_valid_tra(event_id, adjacencies)
        if not result:
            continue # Skip invalid or unsupported (e.g., single breakends)
            
        tra_type, src_chr, (src_start, src_end), snk_chr, snk_pos = result

        # Step 4: Cache sequence and log actions
        # Fetch 0-based sequence between the source breakends
        sequence = ref.fetch(src_chr, src_start, src_end - 1)
        
        if tra_type == "TRA_CUT":
            # For CUT: We delete from source and insert at sink
            tra_map["deletions"][src_chr][event_id] = len(sequence)
            tra_map["insertions"][snk_chr][event_id] = sequence
        else:
            # For COPY: We only insert at sink
            tra_map["insertions"][snk_chr][event_id] = sequence

    vcf.close()
    ref.close()
    return tra_map


valid_tra = prefetch_translocations('test_tra.vcf','data/cerevisiae_test.fa')
total_deletions = sum(len(events) for events in valid_tra['deletions'].values())
total_insertions = sum(len(events) for events in valid_tra['insertions'].values())
print(f"Total Deletion Events: {total_deletions}")
print(f"Total Insertion Events: {total_insertions}")