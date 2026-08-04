import pytest
import pysam
import tempfile
import os
from pathlib import Path
from inSVert import VariantObjects, utils_sim, utils_ins, simulate, insert

def test_insertion_symbolic_default():
    ins = VariantObjects.Insertion('chr1', 100, 50, 'var1', '0|1', 'A')
    assert ins.get_alt() == "<INS>"
    assert ins.ref == "A"

def test_insertion_non_symbolic():
    ins = VariantObjects.Insertion('chr1', 100, 50, 'var1', '0|1', 'A', alt_seq="AGATTACA")
    assert ins.get_alt() == "AGATTACA"
    assert ins.ref == "A"

def test_deletion_non_symbolic():
    del_obj = VariantObjects.Deletion('chr1', 100, -5, 'var2', '0|1', 'A', alt_seq="A", ref_seq="AGATTAC")
    assert del_obj.get_alt() == "A"
    assert del_obj.ref == "AGATTAC"

def test_inversion_non_symbolic():
    inv = VariantObjects.Inversion('chr1', 100, 5, 'var3', '0|1', 'A', alt_seq="ATGTA", ref_seq="ATACA")
    assert inv.get_alt() == "ATGTA"
    assert inv.ref == "ATACA"

def test_duplication_non_symbolic():
    dup = VariantObjects.Duplication('chr1', 100, 5, 'var4', '0|1', 'A', copy_number=2, alt_seq="AGATTAGATTA", ref_seq="AGATTA")
    assert dup.get_alt() == "AGATTAGATTA"
    assert dup.ref == "AGATTA"

class MockRecord:
    def __init__(self, alts, ref="A", svlen=None, var_id="var1", chrom="chr1", pos=100):
        self.alts = alts
        self.ref = ref
        self.id = var_id
        self.chrom = chrom
        self.pos = pos
        self._info = {}
        if svlen is not None:
            self._info['SVLEN'] = svlen
    @property
    def info(self):
        return self._info

def test_extract_explicit_ins_symbolic():
    rec = MockRecord(("<INS>",), svlen=5)
    seq, warning = utils_ins.extract_explicit_ins_sequence(rec)
    assert seq is None
    assert warning is None

def test_extract_explicit_ins_valid():
    rec = MockRecord(("AGATTACA",), ref="A", svlen=7)
    seq, warning = utils_ins.extract_explicit_ins_sequence(rec)
    assert seq == "GATTACA"
    assert warning is None

def test_extract_explicit_ins_length_mismatch():
    rec = MockRecord(("AGATTACA",), ref="A", svlen=10)
    seq, warning = utils_ins.extract_explicit_ins_sequence(rec)
    assert seq == "GATTACA"
    assert warning is not None
    assert "Proceeding using explicit ALT sequence" in warning

def test_extract_explicit_ins_no_svlen():
    rec = MockRecord(("AGATTACA",), ref="A", svlen=None)
    seq, warning = utils_ins.extract_explicit_ins_sequence(rec)
    assert seq == "GATTACA"
    assert warning is None

def test_fetch_ref_span(test_ref_path):
    with pysam.FastaFile(str(test_ref_path)) as ref:
        span = utils_sim.fetch_ref_span("chr1", 1, 10, ref)
        assert len(span) == 10
        assert isinstance(span, str)

def test_end_to_end_non_symbolic_simulation_and_insertion(test_ref_path, tmp_path):
    # 1. Create a minimal YAML config
    config_content = """
genome:
  ploidy: 1
  heterozygosity: 1.0

variants:
  INS:
    count: 2
    distribution: "normal"
    parameters:
      median_length: 20
      sigma: 0
      min_length: 20
      max_length: 20
  DEL:
    count: 2
    distribution: "normal"
    parameters:
      median_length: 20
      sigma: 0
      min_length: 20
      max_length: 20
"""
    config_path = tmp_path / "test_config.yaml"
    config_path.write_text(config_content)
    
    vcf_out = tmp_path / "non_symbolic.vcf"
    
    # 2. Run simulation with non_symbolic=True
    simulate.run(str(config_path), str(test_ref_path), str(vcf_out), seed=42, non_symbolic=True)
    
    assert vcf_out.exists()
    
    # Read VCF and verify non-symbolic ALT entries
    with open(vcf_out, 'r') as f:
        vcf_lines = [line.strip() for line in f if not line.startswith('#')]
    
    assert len(vcf_lines) > 0
    for line in vcf_lines:
        fields = line.split('\t')
        ref_col = fields[3]
        alt_col = fields[4]
        # Should not contain symbolic tags like <INS> or <DEL>
        assert not alt_col.startswith("<")
    
    # 3. Run insertion module using this VCF
    fasta_out = tmp_path / "edited.fasta"
    insert.run(0.41, str(test_ref_path), str(vcf_out), 1, str(fasta_out))
    
    assert fasta_out.exists()
    with pysam.FastaFile(str(fasta_out)) as edited_fa:
        assert len(edited_fa.references) > 0
