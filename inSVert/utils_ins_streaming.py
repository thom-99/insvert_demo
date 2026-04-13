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
    It reconstructs the 'logic' of the BNDs to see if they form a closed loop.
    
    Returns: (type, source_chrom, source_range, sink_chrom, sink_pos) or None
    """
    count = len(adjacencies)
    
    # -------------------------------------------------------------------------
    # CATEGORY 1: COPY & PASTE (4 VCF lines / 2 Adjacencies)
    # -------------------------------------------------------------------------
    if count == 2:
        # Step 1: Map positions to chromosomes to identify 'Source' vs 'Sink'
        chroms = defaultdict(list)
        for r, m in adjacencies:
            chroms[r.chrom].append(r.pos)
            chroms[m.chrom].append(m.pos)
            
        # A valid Copy&Paste must involve exactly 2 chromosomes
        if len(chroms) != 2: 
            return None
        
        # Step 2: Identify the 'Source' (The chrom where the two BNDs define a segment)
        # We look for the chromosome that has more than one unique position
        try:
            source_chrom = next(c for c, p in chroms.items() if len(set(p)) > 1)
        except StopIteration:
            return None # Fails if all BNDs are at the exact same coordinate
            
        sink_chrom = next(c for c in chroms if c != source_chrom)
        
        # Step 3: Calculate length and ensure it's not a 1bp 'ambiguous' segment
        source_pos = sorted(list(set(chroms[source_chrom])))
        source_len = source_pos[1] - source_pos[0]
        
        if source_len <= 1: 
            return None # Filter out 1bp segments to avoid destination/source confusion
        
        # Logically: source_pos[0] is segment start, source_pos[1] is segment end.
        # sink_chrom[0] is the point of insertion.
        return ("TRA_COPY", source_chrom, (source_pos[0], source_pos[1]), 
                sink_chrom, chroms[sink_chrom][0])

    # -------------------------------------------------------------------------
    # CATEGORY 2: CUT & PASTE (6 VCF lines / 3 Adjacencies)
    # -------------------------------------------------------------------------
    elif count == 3:
        # Step 1: Find the 'HEAL' adjacency (The one where both ends stay on the same chrom)
        # This adjacency 'stitches' the hole left by the cut segment.
        heal_adj = [a for a in adjacencies if a[0].chrom == a[1].chrom]
        if len(heal_adj) != 1: 
            return None # A TRA_CUT must have exactly one healing adjacency
        
        source_chrom = heal_adj[0][0].chrom
        # The coordinates of the 'HEAL' adjacency define the boundaries of the cut
        source_pos = sorted([heal_adj[0][0].pos, heal_adj[0][1].pos])
        
        # Calculate true sequence length (bases between the two heal points)
        source_len = source_pos[1] - source_pos[0] - 1
        
        if source_len <= 1: 
            return None # Rigorous filter for ambiguous 1bp segments
        
        # Step 2: Find the 'Sink' using the inter-chromosomal adjacencies (The 'PASTE' lines)
        # We look for an adjacency where the two ends are on DIFFERENT chromosomes.
        try:
            sink_adj = next(a for a in adjacencies if a[0].chrom != a[1].chrom)
        except StopIteration:
            return None
            
        # The sink is whichever chromosome is NOT the source chromosome
        if sink_adj[0].chrom != source_chrom:
            sink_chrom = sink_adj[0].chrom
            sink_pos = sink_adj[0].pos
        else:
            sink_chrom = sink_adj[1].chrom
            sink_pos = sink_adj[1].pos
        
        return ("TRA_CUT", source_chrom, (source_pos[0], source_pos[1]), 
                sink_chrom, sink_pos)

    # Ditch any EVENT that doesn't strictly match 2 or 3 adjacencies
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


valid_tra = prefetch_translocations('data/corrected-cut-paste-tra.vcf','data/cerevisiae_test.fa')
total_deletions = sum(len(events) for events in valid_tra['deletions'].values())
total_insertions = sum(len(events) for events in valid_tra['insertions'].values())
print(f"Total Deletion Events: {total_deletions}")
print(f"Total Insertion Events: {total_insertions}")