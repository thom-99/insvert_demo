from . import utils_ins_streaming
import pysam

class BufferWriter:
    """Helper to maintain 60-char line width for FASTA output."""
    def __init__(self, file_handle, width=60):
        self.fh = file_handle
        self.width = width
        self.buffer = ""

    def write(self, data):
        self.buffer += data
        while len(self.buffer) >= self.width:
            self.fh.write(self.buffer[:self.width] + '\n')
            self.buffer = self.buffer[self.width:]
    
    def flush(self):
        if self.buffer:
            self.fh.write(self.buffer + '\n')
            self.buffer = ""

def run(gc_content, ref_fasta, vcf_file, output_fasta):
    
    print(f'Streaming variants from {vcf_file}...')
    ref = pysam.FastaFile(ref_fasta)
    vcf = pysam.VariantFile(vcf_file)
    
    with open(output_fasta, 'w') as out_f:
        writer = BufferWriter(out_f)
        
        for chrom in ref.references:
            print(f"Processing {chrom}...", end="\r")
            out_f.write(f">{chrom}\n")
            
            try:
                chrom_variants = list(vcf.fetch(chrom))
            except ValueError:
                chrom_variants = []

            ref_pos = 0 
            chrom_len = ref.get_reference_length(chrom)
            
            for var in chrom_variants:
                start = var.pos - 1
                
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

            # 3. Finish Chromosome
            if ref_pos < chrom_len:
                chunk = ref.fetch(chrom, ref_pos, chrom_len)
                writer.write(chunk)
            
            writer.flush()
            
    vcf.close()
    ref.close()
    print(f'\nDone! Output written to {output_fasta}')