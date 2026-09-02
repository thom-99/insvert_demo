import pytest

from inSVert import input_validation


def write_normal_config(
    path, sigma, count=1, median=100, minimum=100, maximum=101
):
    sigma_line = "" if sigma is None else f"      sigma: {sigma}\n"
    path.write_text(
        f"""\
genome:
  ploidy: 2
  heterozygosity: 0.5
variants:
  INS:
    count: {count}
    distribution: normal
    parameters:
      median_length: {median}
{sigma_line}      min_length: {minimum}
      max_length: {maximum}
""",
        encoding="utf-8",
    )


def test_validate_config_rejects_impractical_normal_distribution(tmp_path):
    config_path = tmp_path / "impractical.yaml"
    write_normal_config(config_path, sigma=1_000_000)

    with pytest.raises(ValueError, match="reduce sigma or widen the range"):
        input_validation.validate_config(config_path)


def test_validate_config_accepts_reasonable_normal_distribution(tmp_path):
    config_path = tmp_path / "reasonable.yaml"
    write_normal_config(config_path, sigma=10)

    assert input_validation.validate_config(config_path) is True


def test_validate_config_checks_default_normal_sigma(tmp_path):
    config_path = tmp_path / "impractical_default.yaml"
    write_normal_config(
        config_path,
        sigma=None,
        median=1_000_000,
        minimum=1_000_000,
        maximum=1_000_000,
    )

    with pytest.raises(ValueError, match="reduce sigma or widen the range"):
        input_validation.validate_config(config_path)


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
