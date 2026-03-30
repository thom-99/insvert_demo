from . import utils_ins_streaming
import pysam

class BufferWriter:
    """Helper to maintain 60-char line width for FASTA output."""
    def __init__(self, file_handle, width=60):
        self.fh = file_handle
        self.width = width
        self.buffer = ""

    def write(self, data):
        self.buffer += data.upper()
        while len(self.buffer) >= self.width:
            self.fh.write(self.buffer[:self.width] + '\n')
            self.buffer = self.buffer[self.width:]
    
    def flush(self):
        if self.buffer:
            self.fh.write(self.buffer + '\n')
            self.buffer = ""

def run(gc_content, ref_fasta, vcf_file, ploidy, output_fasta):
    
    print(f'Streaming variants from {vcf_file}...')
    ref = pysam.FastaFile(ref_fasta)
    vcf = pysam.VariantFile(vcf_file)
    tra_cache = utils_ins_streaming.prefetch_translocations(vcf_file, ref_fasta)
    
    with open(output_fasta, 'w') as out_f:
        writer = BufferWriter(out_f)
        
        for haplotype in range(ploidy):
            # Initialize ONCE per haplotype so inter-chromosomal TRA are tracked
            processed_sources = set()
            processed_sinks = set()

            for chrom in ref.references:
                print(f"Processing {chrom} (Haplotype {haplotype+1})...", end="\r")

                if ploidy==1:
                    out_f.write(f">{chrom}\n")
                else:
                    out_f.write(f">Sample#H{haplotype+1}#{chrom}\n")
                
                try:
                    chrom_variants = list(vcf.fetch(chrom))
                except ValueError:
                    chrom_variants = []

                ref_pos = 0 
                chrom_len = ref.get_reference_length(chrom)
                                
                for var in chrom_variants:

                    # only apply ig the allele at the haplotype is 1 
                    if var.samples[0]['GT'][haplotype] != 1:
                        continue 

                    # 2. Rigorous BND Skip Check
                    # Skip if we already handled the "cut" or "paste" for this event 
                    # it needs to be here so I do not write the reference up to the second BNDs (see next check)
                    svtype = var.info.get("SVTYPE")
                    if svtype == 'BND':
                        event_id = var.info.get('EVENT')
                        role = var.info.get('TRA_ROLE')
                        if (role == 'SOURCE' and event_id in processed_sources) or \
                           (role == 'SINK' and event_id in processed_sinks):
                            continue

                    start = var.pos - 1 #VCF is 1-indexed, python is 0-indexed
                    
                    # Check Overlap
                    if start < ref_pos:
                        continue 

                    # 1. Write Reference up to this SV
                    if start > ref_pos:
                        chunk = ref.fetch(chrom, ref_pos, start)
                        writer.write(chunk)
                        ref_pos = start
                    
                    # 2. Dispatch to utils_ins based on Type

                    if svtype != 'BND':

                        svlen = var.info.get("SVLEN")
                        if isinstance(svlen, tuple): svlen = svlen[0]
                        if svlen is None: svlen = var.stop - var.start
                        length = abs(svlen)
                    

                        if svtype == 'INS':
                            ins_seq = utils_ins_streaming.generate_seq(length, gc_content)
                            utils_ins_streaming.apply_insertion(writer, ins_seq)
                            # ref_pos stays same
                            
                        elif svtype == 'DEL':
                            ref_pos = utils_ins_streaming.apply_deletion(ref_pos, length)
                            
                        elif svtype == 'INV':
                            ref_pos = utils_ins_streaming.apply_inversion(ref, chrom, start, length, writer)
                            
                        elif svtype == 'DUP':
                            # Extract CN logic here or inside utils? 
                            # Keeping extraction here is cleaner for the utility function signature
                            sample_name = list(var.samples.keys())[0] if var.samples else None
                            cn = var.samples[sample_name]['CN'] if sample_name else 2
                            
                            ref_pos = utils_ins_streaming.apply_duplication(ref, chrom, start, length, cn, writer)

                    else:
                        # process TRA (BNDs)
                        role = var.info.get('TRA_ROLE')
                        event_id = var.info.get('EVENT')

                        if role == 'SOURCE' and event_id not in processed_sources:
                            cached_seq = tra_cache.get(event_id, "")
                            length_to_skip = len(cached_seq)
                            if length_to_skip > 0:
                                ref_pos = utils_ins_streaming.apply_deletion(ref_pos, length_to_skip)
                                processed_sources.add(event_id)

                        elif role == 'SINK' and event_id not in processed_sinks:
                            cached_seq = tra_cache.get(event_id, "")
                            if cached_seq:
                                utils_ins_streaming.apply_insertion(writer, cached_seq)
                                processed_sinks.add(event_id)
                            

                # 3. Finish Chromosome
                if ref_pos < chrom_len:
                    chunk = ref.fetch(chrom, ref_pos, chrom_len)
                    writer.write(chunk)
                
                writer.flush()
            
    vcf.close()
    ref.close()
    print(f'\nDone! Output written to {output_fasta}')