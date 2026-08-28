import pytest

from inSVert import input_validation


def write_vcf(path, sample_names, alt="T", genotypes=None):
    if genotypes is None:
        genotypes = ["0|1"]

    sample_columns = "\t".join(sample_names)
    genotype_columns = "\t".join(genotypes)
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=10000>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_columns}\n"
    )
    record = f"chr1\t10\tvar1\tC\t{alt}\t.\tPASS\t.\tGT\t{genotype_columns}\n"
    path.write_text(header + record, encoding="utf-8")


def test_validate_vcf_accepts_single_sample_biallelic_record(tmp_path, test_ref_path):
    vcf_path = tmp_path / "valid.vcf"
    write_vcf(vcf_path, ["SAMPLE"])

    input_validation.validate_vcf(vcf_path, test_ref_path, target_ploidy=2)


def test_validate_vcf_rejects_multiple_samples(tmp_path, test_ref_path):
    vcf_path = tmp_path / "multi_sample.vcf"
    write_vcf(vcf_path, ["SAMPLE1", "SAMPLE2"], genotypes=["0|1", "0|0"])

    with pytest.raises(ValueError, match="Multi-sample VCFs are not supported"):
        input_validation.validate_vcf(vcf_path, test_ref_path, target_ploidy=2)


def test_validate_vcf_rejects_multiallelic_record(tmp_path, test_ref_path):
    vcf_path = tmp_path / "multiallelic.vcf"
    write_vcf(vcf_path, ["SAMPLE"], alt="T,G")

    with pytest.raises(ValueError, match="Multiallelic record at chr1:10"):
        input_validation.validate_vcf(vcf_path, test_ref_path, target_ploidy=2)


def test_validate_vcf_rejects_gt_above_first_alt(tmp_path, test_ref_path):
    vcf_path = tmp_path / "unsupported_gt.vcf"
    write_vcf(vcf_path, ["SAMPLE"], genotypes=["0|2"])

    with pytest.raises(ValueError, match="Only complete genotypes using reference allele 0"):
        input_validation.validate_vcf(vcf_path, test_ref_path, target_ploidy=2)


def test_validate_vcf_preserves_ploidy_check(tmp_path, test_ref_path):
    vcf_path = tmp_path / "ploidy_mismatch.vcf"
    write_vcf(vcf_path, ["SAMPLE"])

    with pytest.raises(ValueError, match="Ploidy mismatch at chr1:10"):
        input_validation.validate_vcf(vcf_path, test_ref_path, target_ploidy=1)
