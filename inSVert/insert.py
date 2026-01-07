from . import utils_ins
import pysam

'''
insvert module that performs the insertion of a VCF file into a fasta reference
'''


def run(gc_content, ref_fasta, vcf_file, output_fasta):

    print('loading the refetrence genome...')
    genome = utils_ins.parse_fasta(ref_fasta)

    #print(genome.keys()) #printing chromosomes
    
    print('parsing VCF...')
    vcf = pysam.VariantFile(vcf_file)

    print("processing SVs programmatically...")
    for SV in vcf:

        chrom = SV.chrom
        pos = SV.pos - 1
        svtype = SV.info.get("SVTYPE")

        chrom_seq = genome[chrom]['sequence']
        chrom_offset = genome[chrom]['offset']

        # --- DELETION ---
        if svtype == 'DEL':
            svlen = SV.info.get("SVLEN")
            
            # applying deletions and updating offset
            genome[chrom]['offset'] = utils_ins.apply_deletion(chrom_seq, svlen, pos, chrom_offset)

        # --- INSERTION ---
        elif svtype == 'INS':
            svlen =  SV.info.get("SVLEN")
            # generate sequence using the config GC content
            ins_seq = utils_ins.generate_seq(svlen, gc_content) 

            genome[chrom]['offset'] = utils_ins.apply_insertion(chrom_seq, ins_seq, pos, chrom_offset)

        # --- INVERSION ---
        elif svtype == 'INV':
            svlen = SV.info.get("SVLEN")
            end = pos + svlen
            genome[chrom]['offset'] = utils_ins.apply_inversion(chrom_seq, pos, end, chrom_offset)

        # --- DUPLICATION ---
        elif svtype == 'DUP':
            svlen= SV.info.get("SVLEN")

            # EXTRACT COPY NUMBER
            # PySAM stores samples in a dict-like object. 
            # We grab the first sample (index 0) and look up the 'CN' key.
            sample_name = list(SV.samples.keys())[0]
            copy_number = SV.samples[sample_name]['CN']        

            genome[chrom]['offset'] = utils_ins.apply_duplication(chrom_seq, pos, svlen, copy_number, chrom_offset)

    vcf.close()

 
    print(f'writing edited genome to {output_fasta}')

    with open(output_fasta, 'w') as fasta:
        for chrom, data in genome.items():
            
            fasta.write(f'>{chrom}\n')

            # decode bytearray as string
            sequence = data['sequence'].decode('ascii')

            for i in range(0, len(sequence), 60):
                fasta.write(sequence[i:i+60] + '\n')

    print(f'VCF insertion completed. Output written to {output_fasta}.')