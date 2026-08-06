from . import utils_ins
import pysam
import random
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

def run(gc_content, ref_fasta, vcf_file, ploidy, output_fasta, skip_unphased=False):
    random.seed(42)

    print(f'Streaming variants from {vcf_file}...')
    ref = pysam.FastaFile(ref_fasta)
    sorted_vcf_path = utils_ins.prepare_vcf(vcf_file) #vcf preparation ensuring it is properly sorted
    vcf = pysam.VariantFile(sorted_vcf_path)
    tra_cache = utils_ins.prefetch_translocations(sorted_vcf_path, ref_fasta)
    
    # Combined pass: count total variants AND detect unphased heterozygous genotypes
    total_variants = 0
    unphased_count = 0
    for rec in vcf:
        total_variants += 1
        gt = rec.samples[0]['GT']
        if not rec.samples[0].phased and len(set(gt)) > 1:
            unphased_count += 1
    total_steps = total_variants * ploidy 

    if unphased_count > 0:
        action = "Skipping" if skip_unphased else "Randomly assigning to haplotypes"
        print(
            f"\n⚠ WARNING: {unphased_count} variant(s) in the VCF have unphased "
            f"genotypes (e.g., 0/1). Phased genotypes (e.g., 0|1) are required "
            f"for deterministic haplotype placement. {action}.\n"
        )

    with open(output_fasta, 'w') as out_f:
        writer = BufferWriter(out_f)
        
        with Progress() as progress:
            task = progress.add_task("[cyan]Inserting SVs...", total=total_steps)

            for haplotype in range(ploidy):
                # Initialize ONCE per haplotype so inter-chromosomal TRA are tracked
                processed_sources = set()
                processed_sinks = set()
                # Cache random haplotype assignments for unphased BND EVENTs
                # so all mates in the same event get the same assignment
                unphased_event_assignments = {}

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
                        
                        # HANDLING UNPHASED VARIANTS AND ASSIGNING HAPLOTYPES
                        sample = var.samples[0]
                        is_phased = sample.phased

                        # Distinguish truly ambiguous unphased (heterozygous e.g. 0/1)
                        # from deterministic unphased (homozygous e.g. 1/1 or 0/0).
                        # Homozygous genotypes have identical alleles on every haplotype,
                        # so phasing is irrelevant — treat them like phased records.
                        gt = list(sample['GT'])
                        is_heterozygous = len(set(gt)) > 1

                        if not is_phased and is_heterozygous:
                            if skip_unphased:
                                continue

                            # For BNDs: all mates in the same EVENT must get the
                            # same random assignment to avoid partial translocations
                            svtype_check = var.info.get("SVTYPE")
                            event_id = var.info.get("EVENT") if svtype_check == "BND" else None

                            if event_id and event_id in unphased_event_assignments:
                                # Reuse the cached assignment for this EVENT
                                shuffled_gt = unphased_event_assignments[event_id]
                            else:
                                # Make a fresh random assignment
                                shuffled_gt = list(gt)
                                if 1 not in shuffled_gt:
                                    continue  # No ALT allele (e.g., 0/0)
                                random.shuffle(shuffled_gt)

                                # Cache it for BND EVENT consistency
                                if event_id:
                                    unphased_event_assignments[event_id] = shuffled_gt

                            # Warn the user
                            var_id = var.id if var.id else f"{var.chrom}:{var.pos}"
                            sv_label = var.info.get("SVTYPE", "UNKNOWN")
                            assigned_haps = [str(idx + 1) for idx, allele in enumerate(shuffled_gt) if allele == 1]
                            haps_str = ", ".join(assigned_haps)
                            print(
                                f"WARNING: Variant '{var_id}' ({sv_label} at "
                                f"{var.chrom}:{var.pos}) is unphased. Phasing "
                                f"information is needed for correct haplotype "
                                f"assignment. Randomly assigned to Haplotype(s) "
                                f"{haps_str}."
                            )

                            # Apply only if this haplotype got the ALT allele
                            if shuffled_gt[haplotype] != 1:
                                continue
                        else:
                            # Phased: use existing logic
                            if sample['GT'][haplotype] != 1:
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

                            # Write the padding base before standard SVs
                            if ref_pos == start:
                                writer.write(ref.fetch(chrom, start, start + 1))
                                ref_pos = start + 1

                            # VARIANT PROCESSING
                            if svtype == 'INS':
                                explicit_seq, warning = utils_ins.extract_explicit_ins_sequence(var)
                                if warning:
                                    print(f"\n{warning}")
                                
                                if explicit_seq:
                                    ins_seq = explicit_seq
                                else:
                                    ins_seq = utils_ins.generate_seq(svlen, gc_content)
                                
                                utils_ins.apply_insertion(writer, ins_seq)
                                # ref_pos stays same
                                    
                            elif svtype == 'DEL':
                                ref_pos = utils_ins.apply_deletion(ref_pos, svlen)
                                
                            elif svtype == 'INV':
                                ref_pos = utils_ins.apply_inversion(ref, chrom, ref_pos, svlen, writer)
                                    
                            elif svtype == 'DUP':
                                sample_name = list(var.samples.keys())[0] if var.samples else None
                                cn = var.samples[sample_name]['CN'] if sample_name else 2
                                    
                                ref_pos = utils_ins.apply_duplication(ref, chrom, ref_pos, svlen, cn, writer)

                        if svtype == 'BND':
                            event_id = var.info.get('EVENT')

                            # ACTION: the 'CUT' i.e. source of cut & paste TRAs
                            del_job = tra_cache["deletions"].get(chrom, {}).get(var.pos)
                            # if there is a deletion job to carry out and I have not carried out before
                            if del_job and event_id not in processed_sources:
                                length, _ = del_job

                                # Write the padding base before the deletion
                                if ref_pos == start:
                                    writer.write(ref.fetch(chrom, start, start + 1))
                                    ref_pos = start + 1

                                ref_pos = utils_ins.apply_deletion(ref_pos, length)
                                processed_sources.add(event_id)
                                continue

                            # ACTION: the 'PASTE' i.e. sink of the cut&paste or copy&paste TRAs
                            ins_job = tra_cache["insertions"].get(chrom, {}).get(var.pos)
                            if ins_job and event_id not in processed_sinks:
                                src_chr, s_start, s_end, is_inverted, attach_after, _ = ins_job

                                # Write the padding base before the insertion
                                if attach_after and ref_pos == start:
                                    writer.write(ref.fetch(chrom, start, start + 1))
                                    ref_pos = start + 1

                                # Lazy load sequence from reference genome
                                seq = ref.fetch(src_chr, s_start, s_end)
                                if is_inverted:
                                    seq = utils_ins.reverse_complement(seq)

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
