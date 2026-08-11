"""
Diagnostic script to identify the exact cause of remaining skipped variants
in inSVert insert after running inSVert simulate.
"""

import pysam
import tempfile
import os
from inSVert import simulate, insert, utils_ins

def run_diagnostic():
    config_path = "config.yaml"
    ref_path = "data/human/test_50mb_2chroms.fa"
    
    if not os.path.exists(config_path) or not os.path.exists(ref_path):
        print(f"Skipping diagnostic test: {config_path} or {ref_path} missing.")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        sim_vcf = os.path.join(tmpdir, "simulated.vcf")
        out_fasta = os.path.join(tmpdir, "output.fa")

        # 1. Run simulate
        simulate.run(
            config_path,
            ref_path,
            sim_vcf,
            seed=42,
            excluded_bed=None,
            non_symbolic=False
        )

        # 2. Inspect sorted VCF and simulate insert manually to catch every skip
        sorted_vcf = utils_ins.prepare_vcf(sim_vcf)
        vcf = pysam.VariantFile(sorted_vcf)
        ref = pysam.FastaFile(ref_path)

        tra_cache = utils_ins.prefetch_translocations(sorted_vcf, ref_path)

        skipped_info = []

        for haplotype in range(2):
            for chrom in ref.references:
                ref_pos = 0
                last_setter = None
                try:
                    chrom_variants = list(vcf.fetch(chrom))
                except ValueError:
                    chrom_variants = []

                for var in chrom_variants:
                    sample = var.samples[0]
                    gt = list(sample['GT'])
                    if gt[haplotype] != 1:
                        continue

                    svtype = var.info.get("SVTYPE")
                    start = var.pos - 1

                    if svtype != "BND" and start < ref_pos:
                        skipped_info.append({
                            "haplotype": haplotype + 1,
                            "chrom": chrom,
                            "pos": var.pos,
                            "start": start,
                            "ref_pos": ref_pos,
                            "svtype": svtype,
                            "id": var.id,
                            "last_setter": last_setter
                        })
                        continue

                    if start > ref_pos:
                        ref_pos = start

                    if svtype != "BND":
                        svlen = var.info.get("SVLEN")
                        if isinstance(svlen, tuple): svlen = svlen[0]
                        if svlen is None: svlen = var.stop - var.start
                        svlen = abs(svlen)

                        if ref_pos == start:
                            ref_pos = start + 1

                        if svtype == "INS":
                            pass
                        elif svtype in ("DEL", "INV", "DUP"):
                            prev_pos = ref_pos
                            ref_pos = ref_pos + svlen
                            last_setter = f"{svtype} {var.id} at {var.pos} len={svlen} (ref_pos {prev_pos} -> {ref_pos})"

                    elif svtype == "BND":
                        event_id = var.info.get("EVENT")
                        del_job = tra_cache["deletions"].get(chrom, {}).get(var.pos)
                        if del_job:
                            length, _ = del_job
                            if ref_pos == start:
                                ref_pos = start + 1
                            prev_pos = ref_pos
                            ref_pos = ref_pos + length
                            last_setter = f"TRA_CUT {event_id} at {var.pos} len={length} (ref_pos {prev_pos} -> {ref_pos})"

        vcf.close()
        ref.close()

        print(f"\n==================================================")
        print(f"DIAGNOSTIC SUMMARY: Found {len(skipped_info)} skipped variants")
        print(f"==================================================")
        for i, info in enumerate(skipped_info[:20], 1):
            print(f"{i}. [{info['svtype']}] {info['chrom']}:{info['pos']} (start={info['start']}, current ref_pos={info['ref_pos']}) ID={info['id']} | Last setter: {info['last_setter']}")

        print("\n=== BREAKENDS FOR TRA_CUT.513 ===")
        vcf_raw = pysam.VariantFile(sim_vcf)
        for r in vcf_raw:
            if r.info.get("EVENT") == "inSVert.TRA_CUT.513":
                print(r.id, r.chrom, r.pos, r.info.get("EVENT"), r.alts)
        vcf_raw.close()

if __name__ == "__main__":
    run_diagnostic()
