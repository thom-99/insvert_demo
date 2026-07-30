import pytest
import pysam
import tempfile
import os
from inSVert import utils_sim, utils_ins

def test_generate_tra_bnds_forward():
    """Verify forward BND generation."""
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        # chr1: A C G T A C G T A C G T A C G T A C G T
        # 1-indexed: 1:A, 2:C, 3:G, 4:T, 5:A, 6:C, 7:G, 8:T, 9:A, 10:C
        # chr2: G G C C G G C C G G C C G G C C G G C C
        # 1-indexed: 1:G, 2:G, 3:C, 4:C, 5:G, 6:G, 7:C, 8:C, 9:G, 10:G, 11:C
        f.write(">chr1\nACGTACGTACGTACGTACGT\n>chr2\nGGCCGGCCGGCCGGCCGGCC\n")
        ref_path = f.name

    try:
        pysam.faidx(ref_path)
        ref_fasta = pysam.FastaFile(ref_path)

        bnds = utils_sim.generate_tra_bnds(
            "TRA_COPY", "chr1", 5, "chr2", 10, 5, "event1", "0|1", ref_fasta, is_reverse=False
        )

        assert len(bnds) == 4
        # P1: t[p[ -> ref_dst at chr2:10 is 'G', pos_src+1 = 6
        assert bnds[0].alt_string == "G[chr1:6["
        rev, attach = utils_ins.parse_bnd_orientation(bnds[0].alt_string)
        assert rev is False
        assert attach is True

        # P2: ]p]t -> ref_src_start at chr1:6 is 'C'
        assert bnds[1].alt_string == "]chr2:10]C"
        rev, attach = utils_ins.parse_bnd_orientation(bnds[1].alt_string)
        assert rev is False
        assert attach is False

        # P3: t[p[ -> ref_src_end at chr1:10 is 'C'
        assert bnds[2].alt_string == "C[chr2:11["
        # P4: ]p]t -> ref_dst at chr2:10 is 'G'
        assert bnds[3].alt_string == "]chr1:10]G"

        ref_fasta.close()
    finally:
        os.remove(ref_path)
        if os.path.exists(ref_path + ".fai"):
            os.remove(ref_path + ".fai")

def test_generate_tra_bnds_reverse():
    """Verify reverse BND generation and compatibility with parse_bnd_orientation."""
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(">chr1\nACGTACGTACGTACGTACGT\n>chr2\nGGCCGGCCGGCCGGCCGGCC\n")
        ref_path = f.name

    try:
        pysam.faidx(ref_path)
        ref_fasta = pysam.FastaFile(ref_path)

        bnds = utils_sim.generate_tra_bnds(
            "TRA_COPY", "chr1", 5, "chr2", 10, 5, "event1", "0|1", ref_fasta, is_reverse=True
        )

        assert len(bnds) == 4
        # P1: t]p] (pointing to pos_src + length = 10), ref_dst at chr2:10 is 'G'
        assert bnds[0].alt_string == "G]chr1:10]"
        rev, attach = utils_ins.parse_bnd_orientation(bnds[0].alt_string)
        assert rev is True
        assert attach is True

        # P2: t]p] -> ref_src_end at chr1:10 is 'C'
        assert bnds[1].alt_string == "C]chr2:10]"
        rev, attach = utils_ins.parse_bnd_orientation(bnds[1].alt_string)
        assert rev is True
        assert attach is True

        # P3: [p[t -> ref_dst_plus1 at chr2:11 is 'C'
        assert bnds[2].alt_string == "[chr1:6[C"
        rev, attach = utils_ins.parse_bnd_orientation(bnds[2].alt_string)
        assert rev is True
        assert attach is False

        # P4: [p[t -> ref_src_start at chr1:6 is 'C'
        assert bnds[3].alt_string == "[chr2:11[C"
        rev, attach = utils_ins.parse_bnd_orientation(bnds[3].alt_string)
        assert rev is True
        assert attach is False

        ref_fasta.close()
    finally:
        os.remove(ref_path)
        if os.path.exists(ref_path + ".fai"):
            os.remove(ref_path + ".fai")

def test_parse_config_reverse_ratio(tmp_path):
    config_content = """
genome:
  ploidy: 2
  heterozygosity: 0.5

variants:
  TRA_COPY:
    count: 10
    distribution: "normal"
    reverse_ratio: 0.25
    parameters:
      median_length: 50000
      sigma: 1000
      min_length: 10000
      max_length: 1000000

  TRA_CUT:
    count: 5
    distribution: "normal"
    parameters:
      reverse_ratio: 0.75
      median_length: 50000
      sigma: 1000
      min_length: 10000
      max_length: 1000000
"""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    res = utils_sim.parse_config(str(config_file))
    assert res['variants']['TRA_COPY']['reverse_ratio'] == 0.25
    assert res['variants']['TRA_CUT']['reverse_ratio'] == 0.75

def test_parse_config_invalid_reverse_ratio(tmp_path):
    config_content = """
genome:
  ploidy: 2
  heterozygosity: 0.5

variants:
  TRA_COPY:
    count: 10
    distribution: "normal"
    reverse_ratio: 1.5
    parameters:
      median_length: 50000
      sigma: 1000
      min_length: 10000
      max_length: 1000000
"""
    config_file = tmp_path / "config.yaml"
    config_file.write_text(config_content)

    with pytest.raises(ValueError, match="reverse_ratio' must be between 0.0 and 1.0"):
        utils_sim.parse_config(str(config_file))
