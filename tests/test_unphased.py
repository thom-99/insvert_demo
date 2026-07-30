import pytest
import tempfile
import pathlib
import pysam
import random
from inSVert import insert

@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tmp:
        yield pathlib.Path(tmp)

@pytest.fixture
def sample_ref(temp_dir):
    ref_path = temp_dir / "ref.fa"
    # Create reference FASTA with 1000 bp on chr1
    with open(ref_path, "w") as f:
        f.write(">chr1\n")
        f.write("A" * 1000 + "\n")
    pysam.faidx(str(ref_path))
    return ref_path

def create_vcf(vcf_path, records_lines):
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV Type\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV Length\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n"
        "##INFO=<ID=EVENT,Number=1,Type=String,Description=\"Event ID\">\n"
        "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"Mate ID\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    with open(vcf_path, "w") as f:
        f.write(header)
        for line in records_lines:
            f.write(line + "\n")

def test_unphased_variant_detected_and_warned(sample_ref, temp_dir, capsys):
    vcf_path = temp_dir / "unphased.vcf"
    out_fasta = temp_dir / "out.fa"
    create_vcf(vcf_path, [
        "chr1\t100\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50;END=150\tGT\t0/1"
    ])
    
    insert.run(0.41, str(sample_ref), str(vcf_path), ploidy=2, output_fasta=str(out_fasta), skip_unphased=False)
    captured = capsys.readouterr().out

    assert "1 variant(s) in the VCF have unphased genotypes" in captured
    assert "WARNING: Variant 'del1' (DEL at chr1:100) is unphased" in captured
    assert "Randomly assigned to Haplotype(s)" in captured

def test_unphased_variant_skip_flag(sample_ref, temp_dir, capsys):
    vcf_path = temp_dir / "unphased.vcf"
    out_fasta = temp_dir / "out.fa"
    create_vcf(vcf_path, [
        "chr1\t100\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50;END=150\tGT\t0/1"
    ])

    insert.run(0.41, str(sample_ref), str(vcf_path), ploidy=2, output_fasta=str(out_fasta), skip_unphased=True)
    captured = capsys.readouterr().out

    assert "Skipping" in captured
    # Ensure deleted bases were NOT removed because the unphased variant was skipped
    with open(out_fasta, "r") as f:
        lines = [line.strip() for line in f if not line.startswith(">")]
        genome_seq = "".join(lines)
        # Full chromosome should be unmodified (1000 As)
        assert len(genome_seq) == 2000  # 2 haplotypes of 1000 bp

def test_phased_variant_no_warning(sample_ref, temp_dir, capsys):
    vcf_path = temp_dir / "phased.vcf"
    out_fasta = temp_dir / "out.fa"
    create_vcf(vcf_path, [
        "chr1\t100\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50;END=150\tGT\t0|1"
    ])

    insert.run(0.41, str(sample_ref), str(vcf_path), ploidy=2, output_fasta=str(out_fasta), skip_unphased=False)
    captured = capsys.readouterr().out

    assert "unphased" not in captured

def test_upfront_count_summary(sample_ref, temp_dir, capsys):
    vcf_path = temp_dir / "mixed.vcf"
    out_fasta = temp_dir / "out.fa"
    create_vcf(vcf_path, [
        "chr1\t100\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=110\tGT\t0/1",
        "chr1\t300\tdel2\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=310\tGT\t0/1",
        "chr1\t500\tdel3\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=510\tGT\t0|1"
    ])

    insert.run(0.41, str(sample_ref), str(vcf_path), ploidy=2, output_fasta=str(out_fasta), skip_unphased=False)
    captured = capsys.readouterr().out

    assert "2 variant(s) in the VCF have unphased genotypes" in captured

def test_bnd_event_consistency(sample_ref, temp_dir, capsys):
    vcf_path = temp_dir / "unphased_bnd.vcf"
    out_fasta = temp_dir / "out.fa"
    create_vcf(vcf_path, [
        "chr1\t200\tBND1\tA\tA[chr1:400[\t.\tPASS\tSVTYPE=BND;EVENT=EV1;MATEID=BND2\tGT\t0/1",
        "chr1\t201\tBND4\tA\t]chr1:401]A\t.\tPASS\tSVTYPE=BND;EVENT=EV1;MATEID=BND3\tGT\t0/1",
        "chr1\t400\tBND2\tA\t]chr1:200]A\t.\tPASS\tSVTYPE=BND;EVENT=EV1;MATEID=BND1\tGT\t0/1",
        "chr1\t401\tBND3\tA\tA[chr1:201[\t.\tPASS\tSVTYPE=BND;EVENT=EV1;MATEID=BND4\tGT\t0/1",
    ])

    random.seed(42)
    insert.run(0.41, str(sample_ref), str(vcf_path), ploidy=2, output_fasta=str(out_fasta), skip_unphased=False)
    captured = capsys.readouterr().out

    bnd_warnings = [line for line in captured.splitlines() if "WARNING: Variant 'BND" in line]
    assert len(bnd_warnings) == 8  # 4 BNDs x 2 haplotypes
    
    # Extract the assigned Haplotype for all BND1..BND4 warnings on haplotype 1 run
    # All BND warnings for EV1 must report the exact same assigned haplotype!
    assigned_haps = set()
    for line in bnd_warnings:
        # e.g. "Randomly assigned to Haplotype(s) 1."
        hap_part = line.split("Randomly assigned to Haplotype(s) ")[-1]
        assigned_haps.add(hap_part)
    
    assert len(assigned_haps) == 1, f"Expected consistent haplotype assignment across BND mates, got {assigned_haps}"
