from click.testing import CliRunner

from inSVert import cli as cli_module


def test_split_haplotypes_reports_actual_output_paths(tmp_path, monkeypatch):
    reference = tmp_path / "reference.fa"
    reference.touch()
    vcf = tmp_path / "variants.vcf"
    vcf.touch()
    output = tmp_path / "edited.fasta"

    monkeypatch.setattr(cli_module.input_validation, "validate_fasta", lambda *_: True)
    monkeypatch.setattr(cli_module.input_validation, "validate_vcf", lambda *_: True)
    monkeypatch.setattr(cli_module.input_validation, "validate_output_path", lambda *_: True)
    monkeypatch.setattr(cli_module.insert, "run", lambda *_: None)

    result = CliRunner().invoke(
        cli_module.cli,
        [
            "insert",
            str(reference),
            str(vcf),
            "--ploidy",
            "2",
            "--output",
            str(output),
            "--split-haplotypes",
        ],
    )
    
    unwrapped_output = "".join(result.output.splitlines())

    assert result.exit_code == 0, unwrapped_output
    assert str(tmp_path / "edited_hap1.fasta") in unwrapped_output
    assert str(tmp_path / "edited_hap2.fasta") in unwrapped_output
    assert f"Modified genome saved to {output}" not in unwrapped_output
