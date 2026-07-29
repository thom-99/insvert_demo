'''
UTILITIES FOR THE INSERT MODULE OF INSVERT
(Streaming Version)
'''
import re 
import os
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
    Rigorously validates the biological structure of a translocation event
    using a mate-pairing graph approach.
    
    Supports both Inter-chromosomal and Intra-chromosomal (transposition) events:
    - TRA_COPY: Copy & Paste / Interspersed Duplication (2 adjacencies / 4 BNDs)
    - TRA_CUT: Cut & Paste / Transposition (3 adjacencies / 6 BNDs)

    Args:
        event_id (str): The ID of the translocation event.
        adjacencies (list): A list of tuples [(rec1, mate1), (rec2, mate2), ...]
                            representing paired BND records connected by MATEID.

    Returns:
        tuple or None: 
            (tra_type, source_chrom, (s_start_0idx, s_end_0idx), sink_chrom, sink_pos_1idx, del_trigger_pos_1idx)
            Where:
            - tra_type: "TRA_COPY" or "TRA_CUT"
            - source_chrom: Chromosome of the translocated source segment
            - (s_start_0idx, s_end_0idx): 0-indexed half-open slice [start, end) for ref.fetch()
            - sink_chrom: Chromosome of the insertion target site
            - sink_pos_1idx: 1-indexed insertion position on sink chromosome
            - del_trigger_pos_1idx: 1-indexed VCF position of the HEAL record (H1) that triggers deletion skip (None for TRA_COPY)
            - is_inverted (bool): True if the ALT string specifies reverse-strand joining  
    """
    count = len(adjacencies)
    
    # -------------------------------------------------------------------------
    # CATEGORY 1: COPY & PASTE (4 VCF lines / 2 Adjacencies)
    # -------------------------------------------------------------------------
    if count == 2:
        (r1, m1), (r2, m2) = adjacencies
        
        # CRITICAL TOPOLOGY CHECK:
        # In a 2-adjacency translocation event, each adjacency (r, m) has one breakend
        # on the source side and one breakend on the sink side.
        # Therefore, grouping the records across the two adjacencies yields two candidate sides:
        # Candidate 1: Side A = (r1, r2), Side B = (m1, m2)
        # Candidate 2: Side A = (r1, m2), Side B = (m1, r2)
        candidates = [
            ([r1, r2], [m1, m2]),
            ([r1, m2], [m1, r2])
        ]
        
        valid_pair = None
        for side1, side2 in candidates:
            # Both breakends on a given side must reside on the same chromosome
            if side1[0].chrom == side1[1].chrom and side2[0].chrom == side2[1].chrom:
                d1 = abs(side1[0].pos - side1[1].pos)
                d2 = abs(side2[0].pos - side2[1].pos)
                
                # Distinguish SOURCE vs SINK:
                # - SOURCE side: P2 and P3 bound the copied segment (distance > 1 bp)
                # - SINK side: P1 and P4 bound the insertion site (distance <= 1 bp)
                if (d1 > 1 and d2 <= 1) or (d2 > 1 and d1 <= 1):
                    valid_pair = (side1, side2, d1, d2)
                    break
                    
        if not valid_pair:
            return None
            
        side1, side2, d1, d2 = valid_pair
        if d1 > d2:
            source_bnds, sink_bnds = side1, side2
        else:
            source_bnds, sink_bnds = side2, side1
            
        source_chrom = source_bnds[0].chrom
        s_min = min(source_bnds[0].pos, source_bnds[1].pos)
        s_max = max(source_bnds[0].pos, source_bnds[1].pos)
        
        # Convert to 0-indexed half-open coordinates [start, end) for pysam ref.fetch:
        # P2 is at s_min (1-indexed start base of segment). Index is s_min - 1.
        # P3 is at s_max (1-indexed end base of segment). Half-open end index is s_max.
        s_start_0idx = s_min - 1
        s_end_0idx = s_max
        
        sink_chrom = sink_bnds[0].chrom
        
        # Sort sink_bnds by position to guarantee index 0 is the minimum position
        sink_bnds.sort(key=lambda b: b.pos)
        sink_pos = sink_bnds[0].pos

        # Check ALT bracket orientation for reverse-strand joining and attachment position
        is_inverted = False
        attach_after = True
        if sink_bnds[0].alts:
            rev, attach = parse_bnd_orientation(sink_bnds[0].alts[0])
            if rev is not None:
                is_inverted = rev
                attach_after = attach
        
        return ("TRA_COPY", source_chrom, (s_start_0idx, s_end_0idx), sink_chrom, sink_pos, None, is_inverted, attach_after)

    # -------------------------------------------------------------------------
    # CATEGORY 2: CUT & PASTE (6 VCF lines / 3 Adjacencies)
    # -------------------------------------------------------------------------
    elif count == 3:
        # CRITICAL TOPOLOGY CHECK:
        # The HEAL adjacency (H1, H2) bridges the excised segment on source_chrom.
        # Therefore:
        # 1. H1 and H2 MUST be intrachromosomal (r.chrom == m.chrom).
        # 2. The interval (h_min, h_max) MUST strictly enclose the 2 source breakends
        #    belonging to the PASTE adjacencies, while excluding the sink breakends.
        
        heal_adj = None
        for i, (r, m) in enumerate(adjacencies):
            if r.chrom != m.chrom:
                continue
            
            h_min, h_max = min(r.pos, m.pos), max(r.pos, m.pos)
            paste_adjs = [adj for j, adj in enumerate(adjacencies) if j != i]
            
            # Count breakends in paste_adjs that fall inside vs outside the HEAL interval
            inside_bnds = []
            outside_bnds = []
            for pr, pm in paste_adjs:
                for pb in (pr, pm):
                    if pb.chrom == r.chrom and h_min < pb.pos < h_max:
                        inside_bnds.append(pb)
                    else:
                        outside_bnds.append(pb)
                        
            # In a valid TRA_CUT, exactly 2 paste breakends lie inside the HEAL interval.
            # Additionally, the 2 outside breakends must be adjacent (distance <= 1) 
            # because they represent the single insertion point (sink).
            if len(inside_bnds) == 2 and len(outside_bnds) == 2:
                if outside_bnds[0].chrom == outside_bnds[1].chrom:
                    if abs(outside_bnds[0].pos - outside_bnds[1].pos) <= 1:
                        heal_adj = (r, m)
                        valid_sink_bnds = outside_bnds
                        # Note: For intra-chromosomal transpositions, there are 2 mathematically 
                        # identical interpretations (cut B paste after C == cut C paste before B).
                        # We just take the first one we find that is valid.
                        break
                
        if not heal_adj:
            return None
            
        source_chrom = heal_adj[0].chrom
        h_min = min(heal_adj[0].pos, heal_adj[1].pos)
        h_max = max(heal_adj[0].pos, heal_adj[1].pos)
        
        # Excised source segment in 0-indexed half-open coordinates [start, end):
        # H1 is at h_min (1-indexed base before cut). Index is h_min.
        # H2 is at h_max (1-indexed base after cut). Half-open end index is h_max - 1.
        s_start_0idx = h_min
        s_end_0idx = h_max - 1
        del_trigger_pos = h_min  # 1-indexed VCF position of H1 record triggering deletion skip
        
        # Find sink location from paste adjacencies (breakends outside the HEAL interval or on another chrom)
        paste_adjs = [a for a in adjacencies if a != heal_adj]
        sink_bnds = []
        for pr, pm in paste_adjs:
            for pb in (pr, pm):
                if pb.chrom != source_chrom or pb.pos < h_min or pb.pos > h_max:
                    sink_bnds.append(pb)
                    
        if not sink_bnds:
            return None
            
        sink_chrom = sink_bnds[0].chrom
        
        # Sort sink_bnds by position to guarantee index 0 is the minimum position
        sink_bnds.sort(key=lambda b: b.pos)
        sink_pos = sink_bnds[0].pos

        # Check ALT bracket orientation for reverse-strand joining and attachment position
        is_inverted = False
        attach_after = True
        if sink_bnds[0].alts:
            rev, attach = parse_bnd_orientation(sink_bnds[0].alts[0])
            if rev is not None:
                is_inverted = rev
                attach_after = attach
        
        return ("TRA_CUT", source_chrom, (s_start_0idx, s_end_0idx), sink_chrom, sink_pos, del_trigger_pos, is_inverted, attach_after)

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
    orphan_bnds = 0
    for var in vcf:
        if var.info.get("SVTYPE") == "BND":
            event_tag = var.info.get("EVENT")
            if event_tag is None:
                orphan_bnds += 1
                continue
            events[event_tag].append(var)
    
    if orphan_bnds > 0:
        print(f"WARNING: {orphan_bnds} BND record(s) found without EVENT tags. "
              f"Translocation reconstruction requires EVENT grouping. "
              f"Skipping these BND records.")

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
        
        # Check if the translocation is valid and categorize it
        res = is_valid_tra(event_id, adjacencies)
        if not res: continue
        
        tra_type, src_chr, (s_start, s_end), snk_chr, snk_pos, del_pos, is_inverted, attach_after = res
        sequence = ref.fetch(src_chr, s_start, s_end)
        if is_inverted:
            sequence = reverse_complement(sequence)

        if tra_type == "TRA_CUT" and del_pos is not None:
            tra_map["deletions"][src_chr][del_pos] = (len(sequence), event_id)
        
        tra_map["insertions"][snk_chr][snk_pos] = (sequence, attach_after, event_id)

    vcf.close(); ref.close()
    return tra_map




### FUNCTIONS FOR HANDLING THE INPUT VCF ###

# checks if the vcf file is compressed and indexed, if not it provides to do so and return


''' removes a file and its associated index (.tbi)'''
def cleanup_indexed_vcf(filepath):
    if filepath and os.path.exists(filepath):
        os.remove(filepath)
    if filepath and os.path.exists(filepath + ".tbi"):
        os.remove(filepath + ".tbi")



# 1. Helper: Sorts VCF using Python RAM (No external tools needed)
def sort_vcf(input_path, output_path):
    """
    Reads all variants, sorts them in memory, and writes a compressed VCF.
    """
    
    infile = pysam.VariantFile(input_path)
    header = infile.header
    
    # Create a map of Chromosome -> Order (e.g., chr1=0, chr2=1)
    # This ensures we sort exactly as the VCF header defines the chromosomes
    contig_map = {contig: i for i, contig in enumerate(header.contigs)}
    
    # Read all variants into a list
    try:
        records = list(infile)
    finally:
        infile.close()
    
    # Sort: First by Chromosome ID, then by Position
    records.sort(key=lambda r: (contig_map.get(r.chrom, 9999), r.pos))
    
    # Write to output with 'wz' mode (Force BGZF compression)
    outfile = pysam.VariantFile(output_path, mode='wz', header=header)
    for rec in records:
        outfile.write(rec)
    outfile.close()

# 2. Main Function
def prepare_vcf(vcf_path):
    # Case 1: Already compressed and indexed?
    if vcf_path.endswith(".gz") and os.path.exists(vcf_path + ".tbi"):
        print("Using existing indexed VCF.")
        return vcf_path
    
    # Case 2: We need to sort it. 
    # We define the target name (e.g. data.vcf -> data.vcf.gz)
    sorted_vcf_path = vcf_path + ".gz"
    
    print(f"Sorting and indexing variants to {sorted_vcf_path}...")
    
    try:
        sort_vcf(vcf_path, sorted_vcf_path)
        
        # STEP B: Index using Pysam
        pysam.tabix_index(sorted_vcf_path, preset="vcf", force=True)
        
        return sorted_vcf_path 
    
    except Exception as e:
        # Cleanup if we failed halfway
        if os.path.exists(sorted_vcf_path):
            os.remove(sorted_vcf_path)
        raise RuntimeError(f"Failed to prepare VCF: {e}")







### FUNCTION TO PARSE ORIENTATION OF BNDs IN A VCF using regex ###


def parse_bnd_orientation(alt_string):
    """
    Parses a VCF 4.2 BND ALT string to extract strand orientation 
    and attachment position.

    Args:
        alt_string (str): The ALT string from the VCF record.

    Returns:
        tuple: (needs_rev_comp: bool, attach_after: bool)
               Returns (None, None) if the string is malformed or not a BND.
    """
    # Uses [^\]\[]+ to match any sequence content before/after brackets
    pattern = re.compile(r"^([^\]\[]+)?([\]\[])([^:]+):(\d+)([\]\[])([^\]\[]+)?$")
    match = pattern.match(alt_string)

    if not match:
        return None, None

    seq_before, bracket_1, _, _, bracket_2, seq_after = match.groups()

    # Case 1: t[p[ -> Forward strand, attached AFTER base t
    if seq_before and bracket_1 == '[' and bracket_2 == '[':
        return False, True

    # Case 2: t]p] -> Reverse strand, attached AFTER base t
    elif seq_before and bracket_1 == ']' and bracket_2 == ']':
        return True, True

    # Case 3: ]p]t -> Forward strand, attached BEFORE base t
    elif seq_after and bracket_1 == ']' and bracket_2 == ']':
        return False, False

    # Case 4: [p[t -> Reverse strand, attached BEFORE base t
    elif seq_after and bracket_1 == '[' and bracket_2 == '[':
        return True, False

    return None, None  # Malformed brackets (e.g., t[p])
