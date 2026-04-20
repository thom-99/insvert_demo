from . import utils_ins
import pysam
from rich.progress import Progress

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
    tra_cache = utils_ins.prefetch_translocations(vcf_file, ref_fasta)
    
    #count tot. variants and multiply by ploidy for rich progress bar
    total_variants = sum(1 for _ in vcf)
    total_steps = total_variants * ploidy 

    with open(output_fasta, 'w') as out_f:
        writer = BufferWriter(out_f)
        
        with Progress() as progress:
            task = progress.add_task("[cyan]Inserting SVs...", total=total_steps)

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
                        progress.advance(task)
                        
                        # only apply ig the allele at the haplotype is 1 
                        if var.samples[0]['GT'][haplotype] != 1:
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
                        svtype = var.info.get("SVTYPE")
                        
                        # comoute the length of the variant, BNDs are excluded as they do not have a svlength
                        if svtype != "BND":
                            svlen = var.info.get("SVLEN")
                            if isinstance(svlen, tuple): svlen = svlen[0]
                            if svlen is None: svlen = var.stop - var.start
                            svlen = abs(svlen)

                            # VARIANT PROCESSING
                            if svtype == 'INS':
                                ins_seq = utils_ins.generate_seq(svlen, gc_content)
                                utils_ins.apply_insertion(writer, ins_seq)
                                # ref_pos stays same
                                    
                            elif svtype == 'DEL':
                                ref_pos = utils_ins.apply_deletion(ref_pos, svlen)
                                
                            elif svtype == 'INV':
                                ref_pos = utils_ins.apply_inversion(ref, chrom, start, svlen, writer)
                                    
                            elif svtype == 'DUP':
                                sample_name = list(var.samples.keys())[0] if var.samples else None
                                cn = var.samples[sample_name]['CN'] if sample_name else 2
                                    
                                ref_pos = utils_ins.apply_duplication(ref, chrom, start, svlen, cn, writer)

                        if svtype == 'BND':
                            event_id = var.info.get('EVENT')

                            # ACTION: the 'CUT' i.e. source of cut & paste TRAs
                            del_job = tra_cache["deletions"].get(chrom, {}).get(var.pos)
                            # if there is a deletion job to carry out and I have not carried out before
                            if del_job and event_id not in processed_sources:
                                length, _ = del_job
                                ref_pos = utils_ins.apply_deletion(ref_pos, length)
                                processed_sources.add(event_id)
                                continue

                            # ACTION: the 'PASTE' i.e. source of the cut&paste or copy&paste TRAs
                            ins_job = tra_cache["insertions"].get(chrom, {}).get(var.pos)
                            if ins_job and event_id not in processed_sinks:
                                seq, _ = ins_job
                                utils_ins.apply_insertion(writer, seq)
                                processed_sinks.add(event_id)
                                continue

                                                          

                    # 3. Finish Chromosome
                    if ref_pos < chrom_len:
                        chunk = ref.fetch(chrom, ref_pos, chrom_len)
                        writer.write(chunk)
                    
                    writer.flush()
            
    vcf.close()
    ref.close()
    print(f'\nDone! Output written to {output_fasta}')