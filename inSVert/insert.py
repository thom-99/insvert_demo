from . import utils_ins
import os
import pysam
import random
from collections import defaultdict
from rich.progress import Progress

class BufferWriter:
    """Helper to maintain 60-char line width for FASTA output."""
    def __init__(self, file_handle, width=60):
        self.fh = file_handle
        self.width = width
        self.buffer = ""

    def write(self, data):
        self.buffer += data.upper()
        pos = 0
        lines = []
        while pos + self.width <= len(self.buffer):
            lines.append(self.buffer[pos:pos + self.width])
            pos += self.width
        if lines:
            self.fh.write('\n'.join(lines) + '\n')
        self.buffer = self.buffer[pos:]
    
    def flush(self):
        if self.buffer:
            self.fh.write(self.buffer + '\n')
            self.buffer = ""

def run(gc_content, ref_fasta, vcf_file, ploidy, output_fasta, skip_unphased=False, sample_name="Sample", split_haplotypes=False, truth_vcf=None):
    random.seed(42)

    print(f'Streaming variants from {vcf_file}...')
    ref = pysam.FastaFile(ref_fasta)
    sorted_vcf_path = utils_ins.prepare_vcf(vcf_file) #vcf preparation ensuring it is properly sorted
    vcf = pysam.VariantFile(sorted_vcf_path)

    # A truth VCF must never replace an input or an output used by this run.
    if truth_vcf:
        truth_path = os.path.realpath(os.path.abspath(truth_vcf))
        protected_paths = {
            os.path.realpath(os.path.abspath(ref_fasta)),
            os.path.realpath(os.path.abspath(vcf_file)),
            os.path.realpath(os.path.abspath(sorted_vcf_path)),
            os.path.realpath(os.path.abspath(output_fasta)),
        }
        if split_haplotypes:
            output_basename, output_extension = os.path.splitext(output_fasta)
            protected_paths.update(
                os.path.realpath(os.path.abspath(f"{output_basename}_hap{haplotype + 1}{output_extension}"))
                for haplotype in range(ploidy)
            )
        if truth_path in protected_paths:
            vcf.close()
            ref.close()
            raise ValueError(
                "Truth VCF output must be different from the input VCF, "
                "reference FASTA, and output FASTA."
            )
    tra_cache = utils_ins.prefetch_translocations(sorted_vcf_path, ref_fasta)
    
    # Combined pass: count variants, unphased genotypes, and DUPs lacking CN.
    total_variants = 0
    unphased_count = 0
    missing_cn_dup_count = 0
    for rec in vcf:
        total_variants += 1
        gt = rec.samples[0]['GT']
        if not rec.samples[0].phased and len(set(gt)) > 1:
            unphased_count += 1
        if rec.info.get("SVTYPE") == "DUP" and rec.samples[0].get("CN") is None:
            missing_cn_dup_count += 1
    total_steps = total_variants * ploidy 

    if unphased_count > 0:
        action = "Skipping" if skip_unphased else "Randomly assigning to haplotypes"
        print(
            f"\n⚠ WARNING: {unphased_count} records in the VCF have unphased "
            f"genotypes (e.g., 0/1). Phased genotypes (e.g., 0|1) are required "
            f"for deterministic haplotype placement. {action}.\n")

    # Cache unphased assignments ONCE across all haplotype passes
    unphased_assignments = {}
    # Cache generated sequences for symbolic <INS> so that they are re-used across haplotypes and not randomly re-generated per haplotype
    symbolic_ins_sequences = {}
    # Track skipped positions as a set of (chrom, pos) so each distinct genomic position is counted exactly once, regardless of ploidy.
    skipped_positions = set()
    supported_svtypes = {'INS', 'DEL', 'INV', 'DUP', 'BND'}
    warned_unsupported_svtypes = set()

    # Track the haplotypes on which each ordinary record was expected and
    # actually inserted. A record is truthful only when these sets match.
    variant_haplotypes = defaultdict(lambda: {"expected": set(), "inserted": set()})

    # BND records describe a whole event, so success is tracked by EVENT and
    # by its required cut/paste actions rather than by individual VCF lines.
    bnd_expected_haplotypes = defaultdict(set)
    bnd_required_actions = defaultdict(set)
    bnd_completed_actions = defaultdict(set)
    if truth_vcf:
        for deletion_jobs in tra_cache["deletions"].values():
            for _, event_id in deletion_jobs.values():
                bnd_required_actions[event_id].add("deletion")
        for insertion_jobs in tra_cache["insertions"].values():
            for *_, event_id in insertion_jobs.values():
                bnd_required_actions[event_id].add("insertion")
    
    # if --split-haplotypes is requested, create the N output files and open them for writing
    output_basename, output_extension = os.path.splitext(output_fasta)
    output_paths = [f"{output_basename}_hap{haplotype + 1}{output_extension}"for haplotype in range(ploidy)] if split_haplotypes else [output_fasta]
    output_files = [open(path, 'w') for path in output_paths]
    
    try:
        with Progress() as progress:
            task = progress.add_task("[cyan]Inserting SVs...", total=total_steps)

            for haplotype in range(ploidy):
                
                # set output fasta file according to haplotype
                out_f = output_files[haplotype] if split_haplotypes else output_files[0]
                writer = BufferWriter(out_f)

                # Initialize ONCE per haplotype so inter-chromosomal TRA are tracked
                processed_sources = set()
                processed_sinks = set()

                for chrom in ref.references:
                    progress.update(task, description=f"Processing {chrom} (Haplotype {haplotype+1})...")

                    if ploidy==1:
                        out_f.write(f">{chrom}\n")
                    else:
                        out_f.write(f">{sample_name}#H{haplotype+1}#{chrom}\n")

                    try:
                        chrom_variants = vcf.fetch(chrom)
                    except ValueError:
                        chrom_variants = []

                    ref_pos = 0
                    chrom_len = ref.get_reference_length(chrom)

                    for variant_index, var in enumerate(chrom_variants):
                        progress.advance(task)

                        svtype = var.info.get("SVTYPE")
                        if svtype is not None and svtype not in supported_svtypes:
                            if svtype not in warned_unsupported_svtypes:
                                progress.console.print(
                                    f"\nWARNING: Unsupported SVTYPE '{svtype}' "
                                    f"detected at {var.chrom}:{var.pos}; records of this type "
                                    "will be skipped."
                                )
                                warned_unsupported_svtypes.add(svtype)
                            continue

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

                            # Universal key: EVENT tag for BNDs (so all mates
                            # in the same event get the same assignment),
                            # (chrom, pos) for everything else.
                            svtype_check = var.info.get("SVTYPE")
                            if svtype_check == "BND":
                                var_key = var.info.get("EVENT")
                            else:
                                var_key = (var.chrom, var.pos)

                            if var_key in unphased_assignments:
                                shuffled_gt = unphased_assignments[var_key]
                            else:
                                # Make a fresh random assignment
                                shuffled_gt = list(gt)
                                if 1 not in shuffled_gt:
                                    continue  # No ALT allele (e.g., 0/0)
                                random.shuffle(shuffled_gt)
                                unphased_assignments[var_key] = shuffled_gt


                            # Apply only if this haplotype got the ALT allele
                            if shuffled_gt[haplotype] != 1:
                                continue
                        else:
                            # Phased: use existing logic
                            if haplotype >= len(sample['GT']) or sample['GT'][haplotype] != 1:
                                continue




                        # Record the intended ALT placement before later checks
                        # can reject it as malformed or overlapping.
                        if truth_vcf:
                            variant_key = (chrom, variant_index)
                            if svtype == "BND":
                                event_id = var.info.get("EVENT")
                                if event_id is not None:
                                    bnd_expected_haplotypes[event_id].add(haplotype)
                            else:
                                variant_haplotypes[variant_key]["expected"].add(haplotype)

                        start = var.pos - 1 # VCF is 1-indexed, python is 0-indexed

                        # Check Overlap only for non-BND variants.
                        # BND companion records intentionally share coordinates with
                        # their active sibling (the one that advances ref_pos). They
                        # are already deduplicated via processed_sources/processed_sinks,
                        # so applying start < ref_pos here would spuriously skip them.
                        if svtype != "BND" and start < ref_pos:
                            skipped_positions.add((haplotype, chrom, var.pos))
                            continue

                        # 1. Write Reference up to this variant
                        if start > ref_pos:
                            chunk = ref.fetch(chrom, ref_pos, start)
                            writer.write(chunk)
                            ref_pos = start

                        # 2. Dispatch to utils_ins based on Type

                        ### SMALL VARIANTS HANDLING (SNPs, MNPs, indels) ###

                        if svtype is None:

                            ref_allele = var.ref or ""
                            alt_str = str(var.alts[0]) if var.alts else ""

                            if not ref_allele or not alt_str:
                                continue # malformed record, skipping

                            # Write ALT str, consume REF lenght from reference
                            writer.write(alt_str)
                            ref_pos = start + len(ref_allele)
                            if truth_vcf:
                                variant_haplotypes[variant_key]["inserted"].add(haplotype)
                            continue


                        ### STRUCTURAL VARIANTS HANDLING (INS, DEL, INV, DUP, BNDs) ###

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

                            # SV PROCESSING
                            if svtype == 'INS':
                                explicit_seq, warning = utils_ins.extract_explicit_ins_sequence(var)
                                if warning:
                                    progress.console.print(f"\n{warning}")

                                if explicit_seq:
                                    ins_seq = explicit_seq
                                else:
                                    insertion_key = (var.chrom, var.pos, var.id)
                                    # generate the sequence for the symbolic <INS> only once and cache it
                                    if insertion_key not in symbolic_ins_sequences:
                                        symbolic_ins_sequences[insertion_key] = utils_ins.generate_seq(svlen, gc_content)

                                    ins_seq = symbolic_ins_sequences[insertion_key]

                                utils_ins.apply_insertion(writer, ins_seq)
                                # ref_pos stays same

                            elif svtype == 'DEL':
                                ref_pos = utils_ins.apply_deletion(ref_pos, svlen)

                            elif svtype == 'INV':
                                ref_pos = utils_ins.apply_inversion(ref, chrom, ref_pos, svlen, writer)

                            elif svtype == 'DUP':
                                vcf_sample_name = list(var.samples.keys())[0] if var.samples else None
                                cn = var.samples[vcf_sample_name].get('CN') if vcf_sample_name else None
                                if cn is None:
                                    cn = 2

                                ref_pos = utils_ins.apply_duplication(ref, chrom, ref_pos, svlen, cn, writer)

                            # Reaching this point means the selected structural
                            # variant was applied without being skipped.
                            if truth_vcf:
                                variant_haplotypes[variant_key]["inserted"].add(haplotype)

                        if svtype == 'BND':
                            event_id = var.info.get('EVENT')

                            # ACTION: the 'CUT' i.e. source of cut & paste TRAs
                            del_job = tra_cache["deletions"].get(chrom, {}).get((var.pos, event_id))
                            # if there is a deletion job to carry out and I have not carried out before
                            if del_job and event_id not in processed_sources:
                                length, _ = del_job

                                # Write the padding base before the deletion
                                if ref_pos == start:
                                    writer.write(ref.fetch(chrom, start, start + 1))
                                    ref_pos = start + 1

                                ref_pos = utils_ins.apply_deletion(ref_pos, length)
                                processed_sources.add(event_id)
                                if truth_vcf:
                                    bnd_completed_actions[(haplotype, event_id)].add("deletion")
                                continue

                            # ACTION: the 'PASTE' i.e. sink of the cut&paste or copy&paste TRAs
                            ins_job = tra_cache["insertions"].get(chrom, {}).get((var.pos, event_id))
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
                                if truth_vcf:
                                    bnd_completed_actions[(haplotype, event_id)].add("insertion")
                                continue



                    # 3. Finish Chromosome
                    if ref_pos < chrom_len:
                        chunk = ref.fetch(chrom, ref_pos, chrom_len)
                        writer.write(chunk)

                    writer.flush()
    
    finally:
        for output_file in output_files:
            output_file.close() # close all output files regardless of insertion errors
    
    if skipped_positions:
        print(f"Warning: skipped {len(skipped_positions)} variant record(s) due to coordinate overlaps.")
    if missing_cn_dup_count:
        print(f"WARNING: {missing_cn_dup_count} DUP record(s) had no copy-number (CN) value and were applied with the default CN=2.")

    vcf.close()
    ref.close()

    if truth_vcf:
        successful_variants = {
            key
            for key, state in variant_haplotypes.items()
            if state["expected"] and state["inserted"] == state["expected"]
        }
        successful_bnd_events = {
            event_id
            for event_id, haplotypes in bnd_expected_haplotypes.items()
            if haplotypes
            and bnd_required_actions[event_id]
            and all(
                bnd_required_actions[event_id]
                <= bnd_completed_actions[(haplotype, event_id)]
                for haplotype in haplotypes
            )
        }

        # Re-read the prepared VCF so records and the complete original header
        # are copied unchanged while unsuccessful records are filtered out.
        per_chrom_index = defaultdict(int)
        output_mode = "wz" if os.fspath(truth_vcf).endswith(".gz") else "w"
        with pysam.VariantFile(sorted_vcf_path) as source_vcf:
            with pysam.VariantFile(
                truth_vcf, mode=output_mode, header=source_vcf.header.copy()
            ) as output_vcf:
                for record in source_vcf:
                    record_index = per_chrom_index[record.chrom]
                    per_chrom_index[record.chrom] += 1

                    if record.info.get("SVTYPE") == "BND":
                        keep = record.info.get("EVENT") in successful_bnd_events
                    else:
                        keep = (record.chrom, record_index) in successful_variants

                    if keep:
                        output_vcf.write(record)

        print(f"Truth VCF written to {truth_vcf}")

    if split_haplotypes:
        print(f'\nDone! Outputs written to {", ".join(output_paths)}')
    else:
        print(f'\nDone! Output written to {output_fasta}')
