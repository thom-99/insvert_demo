import subprocess
from pathlib import Path

import pysam
import pytest

from inSVert import insert


def create_reference(path, contigs):
    with open(path, "w") as reference:
        for name, length in contigs.items():
            reference.write(f">{name}\n")
            reference.write("A" * length + "\n")
    pysam.faidx(str(path))


def create_vcf(path, records):
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=truth-vcf-test\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV Type\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV Length\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    with open(path, "w") as vcf:
        vcf.write(header)
        for record in records:
            vcf.write(record + "\n")


def test_truth_vcf_keeps_only_fully_inserted_records(tmp_path):
    reference_path = tmp_path / "reference.fa"
    input_vcf_path = tmp_path / "input.vcf"
    output_fasta_path = tmp_path / "output.fa"
    truth_vcf_path = tmp_path / "truth.vcf"
    create_reference(reference_path, {"chr1": 1000})
    create_vcf(input_vcf_path, [
        "chr1\t100\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50;END=150\tGT\t1|0",
        "chr1\t120\tdel_partial\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=130\tGT\t1|1",
        "chr1\t300\tdel_absent\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=310\tGT\t0|0",
        "chr1\t500\tsnp1\tA\tT\t.\tPASS\t.\tGT\t1|1",
    ])

    insert.run(
        0.41,
        str(reference_path),
        str(input_vcf_path),
        ploidy=2,
        output_fasta=str(output_fasta_path),
        truth_vcf=str(truth_vcf_path),
    )

    with pysam.VariantFile(str(input_vcf_path)) as input_vcf:
        input_header = str(input_vcf.header)
    with pysam.VariantFile(str(truth_vcf_path)) as truth_vcf:
        assert str(truth_vcf.header) == input_header
        assert [record.id for record in truth_vcf] == ["del1", "snp1"]


def test_truth_vcf_excludes_unphased_records_when_skipped(tmp_path):
    reference_path = tmp_path / "reference.fa"
    input_vcf_path = tmp_path / "input.vcf"
    output_fasta_path = tmp_path / "output.fa"
    truth_vcf_path = tmp_path / "truth.vcf"
    create_reference(reference_path, {"chr1": 1000})
    create_vcf(input_vcf_path, [
        "chr1\t100\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=110\tGT\t0/1",
    ])

    insert.run(
        0.41,
        str(reference_path),
        str(input_vcf_path),
        ploidy=2,
        output_fasta=str(output_fasta_path),
        skip_unphased=True,
        truth_vcf=str(truth_vcf_path),
    )

    with pysam.VariantFile(str(truth_vcf_path)) as truth_vcf:
        assert list(truth_vcf) == []


def test_truth_vcf_line_count_diff_matches_skipped_variants(tmp_path):
    reference_path = tmp_path / "reference.fa"
    input_vcf_path = tmp_path / "input.vcf"
    output_fasta_path = tmp_path / "output.fa"
    truth_vcf_path = tmp_path / "truth.vcf"
    create_reference(reference_path, {"chr1": 1000})
    create_vcf(input_vcf_path, [
        "chr1\t100\tdel1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50;END=150\tGT\t1",
        "chr1\t120\tdel_skipped\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-10;END=130\tGT\t1",
        "chr1\t300\tsnp1\tA\tT\t.\tPASS\t.\tGT\t1",
    ])

    insert.run(
        0.41,
        str(reference_path),
        str(input_vcf_path),
        ploidy=1,
        output_fasta=str(output_fasta_path),
        truth_vcf=str(truth_vcf_path),
    )

    input_lines = int(
        subprocess.check_output(
            ["wc", "-l", input_vcf_path], text=True
        ).split()[0]
    )
    truth_lines = int(
        subprocess.check_output(
            ["wc", "-l", truth_vcf_path], text=True
        ).split()[0]
    )

    assert input_lines - truth_lines == 1


@pytest.mark.parametrize(
    "fixture_name", ["tra_copy_intra.vcf", "tra_cut_intra.vcf"]
)
def test_truth_vcf_keeps_all_records_of_completed_bnd_event(tmp_path, fixture_name):
    reference_path = tmp_path / "reference.fa"
    input_vcf_path = tmp_path / "bnd.vcf"
    output_fasta_path = tmp_path / "output.fa"
    truth_vcf_path = tmp_path / "truth.vcf"
    create_reference(reference_path, {"chr1": 10000, "chr2": 10000})

    source_vcf = Path(__file__).parent / "test_data" / fixture_name
    input_vcf_path.write_text(source_vcf.read_text())

    insert.run(
        0.41,
        str(reference_path),
        str(input_vcf_path),
        ploidy=2,
        output_fasta=str(output_fasta_path),
        truth_vcf=str(truth_vcf_path),
    )

    with pysam.VariantFile(str(input_vcf_path)) as input_vcf:
        expected_ids = [record.id for record in input_vcf]
    with pysam.VariantFile(str(truth_vcf_path)) as truth_vcf:
        assert [record.id for record in truth_vcf] == expected_ids

