import pytest
import os
import random
import pysam
from inSVert.VariantObjects import SNP
from inSVert.utils_sim import parse_config, pick_snp_alt
from inSVert import simulate, insert


def test_snp_format():
    snp = SNP('chr1', 12345, 'inSVert.SNP.1', '0|1', 'A', 'G')
    assert snp.get_end() == 12345
    formatted = snp.format()
    assert formatted == "chr1\t12345\tinSVert.SNP.1\tA\tG\t.\tPASS\tVT=SNP\tGT\t0|1"


def test_pick_snp_alt_transitions():
    random.seed(42)
    # High Ts/Tv ratio ensures transitions
    assert pick_snp_alt('A', 1000.0) == 'G'
    assert pick_snp_alt('G', 1000.0) == 'A'
    assert pick_snp_alt('C', 1000.0) == 'T'
    assert pick_snp_alt('T', 1000.0) == 'C'


def test_pick_snp_alt_transversions():
    random.seed(42)
    # tstv_ratio = 0 forces transversions only
    assert pick_snp_alt('A', 0.0) in ['C', 'T']
    assert pick_snp_alt('G', 0.0) in ['C', 'T']
    assert pick_snp_alt('C', 0.0) in ['A', 'G']
    assert pick_snp_alt('T', 0.0) in ['A', 'G']


def test_pick_snp_alt_non_acgt():
    assert pick_snp_alt('N', 2.0) is None
    assert pick_snp_alt('R', 2.0) is None


def test_parse_config_snp(tmp_path):
    config_file = tmp_path / "config.yaml"
    config_file.write_text("""
genome:
  ploidy: 2
  heterozygosity: 0.5
variants:
  SNP:
    ratio: 0.001
    tstv_ratio: 2.5
""")
    res = parse_config(str(config_file))
    assert 'SNP' in res['variants']
    assert res['variants']['SNP']['ratio'] == 0.001
    assert res['variants']['SNP']['tstv_ratio'] == 2.5


def test_parse_config_snp_missing_ratio(tmp_path):
    config_file = tmp_path / "config.yaml"
    config_file.write_text("""
genome:
  ploidy: 2
  heterozygosity: 0.5
variants:
  SNP:
    tstv_ratio: 2.0
""")
    with pytest.raises(ValueError, match="SNP requires a 'ratio' parameter"):
        parse_config(str(config_file))


def test_parse_config_snp_invalid_ratio(tmp_path):
    config_file = tmp_path / "config.yaml"
    config_file.write_text("""
genome:
  ploidy: 2
  heterozygosity: 0.5
variants:
  SNP:
    ratio: 1.5
""")
    with pytest.raises(ValueError, match="ratio.*must be between 0 and 1"):
        parse_config(str(config_file))


def test_end_to_end_snp_simulation_and_insertion(tmp_path):
    # 1. Create a dummy reference FASTA
    fasta_path = tmp_path / "ref.fa"
    ref_seq = "ACGT" * 2500  # 10,000 bp
    fasta_path.write_text(f">chr1\n{ref_seq}\n")
    pysam.faidx(str(fasta_path))

    # 2. Create config file with SNPs and DEL
    config_file = tmp_path / "config.yaml"
    config_file.write_text("""
genome:
  ploidy: 1
  heterozygosity: 1.0
variants:
  DEL:
    count: 2
    distribution: "normal"
    parameters:
      median_length: 100
      sigma: 10
      min_length: 50
      max_length: 200
  SNP:
    ratio: 0.005  # ~50 SNPs
    tstv_ratio: 2.0
""")

    vcf_path = tmp_path / "out.vcf"
    out_fasta_path = tmp_path / "modified.fa"

    # 3. Simulate
    simulate.run(str(config_file), str(fasta_path), str(vcf_path), seed=123)
    assert os.path.exists(str(vcf_path))

    # Parse simulated VCF to verify SNPs were generated
    vcf_reader = pysam.VariantFile(str(vcf_path))
    snps = [rec for rec in vcf_reader if rec.info.get('VT') == 'SNP']
    assert len(snps) > 0

    # 4. Insert variants into FASTA
    insert.run(0.5, str(fasta_path), str(vcf_path), 1, str(out_fasta_path))
    assert os.path.exists(str(out_fasta_path))

    # Read modified FASTA and verify length change matches deletions
    mod_fa = pysam.FastaFile(str(out_fasta_path))
    mod_seq = mod_fa.fetch("chr1")
    # Output length should be 10000 minus deleted lengths
    assert len(mod_seq) < 10000
