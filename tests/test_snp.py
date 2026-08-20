import random

import pysam
import pytest

from inSVert import insert, simulate
from inSVert.VariantObjects import Polymorphism
from inSVert.utils_sim import (
    MNP_LENGTHS,
    MNP_WEIGHTS,
    parse_config,
    pick_mnp_alt,
    pick_mnp_length,
    pick_snp_alt,
)


def test_polymorphism_formats_snp_and_mnp():
    snp = Polymorphism('chr1', 12345, 'inSVert.SNP.1', '0|1', 'A', 'G')
    assert snp.get_end() == 12345
    assert snp.format() == "chr1\t12345\tinSVert.SNP.1\tA\tG\t.\tPASS\tVT=SNP\tGT\t0|1"

    mnp = Polymorphism('chr1', 12345, 'inSVert.MNP.1', '0|1', 'ACG', 'GTA')
    assert mnp.get_end() == 12347
    assert mnp.format() == "chr1\t12345\tinSVert.MNP.1\tACG\tGTA\t.\tPASS\tVT=MNP\tGT\t0|1"


@pytest.mark.parametrize(('ref', 'alt'), [('', ''), ('A', ''), ('AC', 'T')])
def test_polymorphism_rejects_invalid_alleles(ref, alt):
    with pytest.raises(ValueError):
        Polymorphism('chr1', 1, 'variant.1', '1', ref, alt)


def test_pick_snp_alt_transitions():
    random.seed(42)
    assert pick_snp_alt('A', 1000.0) == 'G'
    assert pick_snp_alt('G', 1000.0) == 'A'
    assert pick_snp_alt('C', 1000.0) == 'T'
    assert pick_snp_alt('T', 1000.0) == 'C'


def test_pick_snp_alt_transversions():
    random.seed(42)
    assert pick_snp_alt('A', 0.0) in ['C', 'T']
    assert pick_snp_alt('G', 0.0) in ['C', 'T']
    assert pick_snp_alt('C', 0.0) in ['A', 'G']
    assert pick_snp_alt('T', 0.0) in ['A', 'G']


def test_pick_snp_alt_non_acgt():
    assert pick_snp_alt('N', 2.0) is None
    assert pick_snp_alt('R', 2.0) is None


def test_pick_mnp_alt():
    random.seed(42)
    ref = 'ACGT'
    alt = pick_mnp_alt(ref, 2.0)
    assert len(alt) == len(ref)
    assert all(ref_base != alt_base for ref_base, alt_base in zip(ref, alt))
    assert pick_mnp_alt('AN', 2.0) is None


def test_pick_mnp_length_uses_internal_distribution(monkeypatch):
    captured = {}

    def fake_choices(population, weights, k):
        captured.update(population=population, weights=weights, k=k)
        return [3]

    monkeypatch.setattr(random, 'choices', fake_choices)
    assert pick_mnp_length() == 3
    assert captured == {
        'population': MNP_LENGTHS,
        'weights': MNP_WEIGHTS,
        'k': 1,
    }


def test_parse_config_polymorphisms(tmp_path):
    config_file = tmp_path / "config.yaml"
    config_file.write_text("""
genome:
  ploidy: 2
  heterozygosity: 0.5
variants:
  SNP:
    count: 10
    tstv_ratio: 2.5
  MNP:
    count: 20
""")
    variants = parse_config(str(config_file))['variants']
    assert variants['SNP'] == {'count': 10, 'tstv_ratio': 2.5}
    assert variants['MNP'] == {'count': 20, 'tstv_ratio': 2.0}


@pytest.mark.parametrize(
    ('variant_type', 'settings', 'message'),
    [
        ('SNP', 'tstv_ratio: 2.0', "SNP requires a 'count' parameter"),
        ('MNP', 'count: 1.5', "MNP requires a 'count' parameter"),
        ('MNP', 'count: 10\n    tstv_ratio: -1', "MNP 'tstv_ratio' must be >= 0"),
    ],
)
def test_parse_config_rejects_invalid_polymorphisms(tmp_path, variant_type, settings, message):
    config_file = tmp_path / "config.yaml"
    config_file.write_text(f"""
genome:
  ploidy: 2
  heterozygosity: 0.5
variants:
  {variant_type}:
    {settings}
""")
    with pytest.raises(ValueError, match=message):
        parse_config(str(config_file))


def test_end_to_end_polymorphism_simulation(tmp_path, capsys):
    fasta_path = tmp_path / "ref.fa"
    ref_seq = "ACGT" * 2500
    fasta_path.write_text(f">chr1\n{ref_seq}\n")
    pysam.faidx(str(fasta_path))

    config_file = tmp_path / "config.yaml"
    config_file.write_text("""
genome:
  ploidy: 1
  heterozygosity: 1.0
variants:
  SNP:
    count: 20
    tstv_ratio: 2.0
  MNP:
    count: 20
    tstv_ratio: 2.0
""")

    vcf_path = tmp_path / "out.vcf"
    simulate.run(str(config_file), str(fasta_path), str(vcf_path), seed=123)
    output = capsys.readouterr().out
    assert "Generating DNA polymorphisms" in output
    assert "Placed " in output and " DNA polymorphisms" in output
    assert "Generating SNPs" not in output

    with pysam.VariantFile(str(vcf_path)) as vcf_reader:
        records = list(vcf_reader)

    snps = [record for record in records if record.info.get('VT') == 'SNP']
    mnps = [record for record in records if record.info.get('VT') == 'MNP']
    assert snps
    assert mnps
    assert all(record.id.startswith('inSVert.SNP.') for record in snps)
    assert all(record.id.startswith('inSVert.MNP.') for record in mnps)
    assert all(len(record.ref) == len(record.alts[0]) == 1 for record in snps)
    assert all(len(record.ref) == len(record.alts[0]) in MNP_LENGTHS for record in mnps)
    assert all(record.pos + len(record.ref) - 1 <= len(ref_seq) for record in records)


def test_insert_applies_mnp_without_length_change(tmp_path):
    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(">chr1\nAACCGGTT\n")
    pysam.faidx(str(fasta_path))

    vcf_path = tmp_path / "mnp.vcf"
    vcf_path.write_text("""##fileformat=VCFv4.2
##contig=<ID=chr1,length=8>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">
##INFO=<ID=VT,Number=1,Type=String,Description="Variant Type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t3\tinSVert.MNP.1\tCC\tTT\t.\tPASS\tVT=MNP\tGT\t1
""")

    output_path = tmp_path / "modified.fa"
    insert.run(0.5, str(fasta_path), str(vcf_path), 1, str(output_path))

    with pysam.FastaFile(str(output_path)) as modified:
        sequence = modified.fetch('chr1')
    assert sequence == 'AATTGGTT'
    assert len(sequence) == 8
