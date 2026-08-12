import pytest
import pysam
from inSVert.utils_ins import prefetch_translocations, reverse_complement

def test_prefetch_tra_copy_inter(test_data_dir, test_ref_path):
    vcf_path = str(test_data_dir / "tra_copy_inter.vcf")
    tra_map = prefetch_translocations(vcf_path, str(test_ref_path))
    
    assert 'chr1' in tra_map['insertions']
    assert (500, 'COPY1') in tra_map['insertions']['chr1']
    src_chr, s_start, s_end, is_inverted, attach_after, event_id = tra_map['insertions']['chr1'][(500, 'COPY1')]
    assert event_id == 'COPY1'
    ref = pysam.FastaFile(str(test_ref_path))
    sequence = ref.fetch(src_chr, s_start, s_end)
    if is_inverted:
        sequence = reverse_complement(sequence)
    assert sequence == 'TGCA' * 250  # 1000bp
    
    assert not tra_map['deletions']

def test_prefetch_tra_copy_intra(test_data_dir, test_ref_path):
    vcf_path = str(test_data_dir / "tra_copy_intra.vcf")
    tra_map = prefetch_translocations(vcf_path, str(test_ref_path))
    
    assert 'chr1' in tra_map['insertions']
    assert (500, 'COPY2') in tra_map['insertions']['chr1']
    src_chr, s_start, s_end, is_inverted, attach_after, event_id = tra_map['insertions']['chr1'][(500, 'COPY2')]
    assert event_id == 'COPY2'
    ref = pysam.FastaFile(str(test_ref_path))
    sequence = ref.fetch(src_chr, s_start, s_end)
    if is_inverted:
        sequence = reverse_complement(sequence)
    assert sequence == 'ACGT' * 250  # 1000bp
    
    assert not tra_map['deletions']

def test_prefetch_tra_cut_inter(test_data_dir, test_ref_path):
    vcf_path = str(test_data_dir / "tra_cut_inter.vcf")
    tra_map = prefetch_translocations(vcf_path, str(test_ref_path))
    
    assert 'chr1' in tra_map['insertions']
    assert (500, 'CUT1') in tra_map['insertions']['chr1']
    src_chr, s_start, s_end, is_inverted, attach_after, event_id = tra_map['insertions']['chr1'][(500, 'CUT1')]
    assert event_id == 'CUT1'
    ref = pysam.FastaFile(str(test_ref_path))
    sequence = ref.fetch(src_chr, s_start, s_end)
    if is_inverted:
        sequence = reverse_complement(sequence)
    assert sequence == 'TGCA' * 250  # 1000bp
    
    assert 'chr2' in tra_map['deletions']
    assert (1000, 'CUT1') in tra_map['deletions']['chr2']
    assert tra_map['deletions']['chr2'][(1000, 'CUT1')] == (1000, 'CUT1')

def test_prefetch_tra_cut_intra(test_data_dir, test_ref_path):
    vcf_path = str(test_data_dir / "tra_cut_intra.vcf")
    tra_map = prefetch_translocations(vcf_path, str(test_ref_path))
    
    assert 'chr1' in tra_map['insertions']
    assert (7000, 'CUT2') in tra_map['insertions']['chr1']
    src_chr, s_start, s_end, is_inverted, attach_after, event_id = tra_map['insertions']['chr1'][(7000, 'CUT2')]
    assert event_id == 'CUT2'
    ref = pysam.FastaFile(str(test_ref_path))
    sequence = ref.fetch(src_chr, s_start, s_end)
    if is_inverted:
        sequence = reverse_complement(sequence)
    assert sequence == 'ACGT' * 250  # 1000bp
    
    assert 'chr1' in tra_map['deletions']
    assert (3000, 'CUT2') in tra_map['deletions']['chr1']
    assert tra_map['deletions']['chr1'][(3000, 'CUT2')] == (1000, 'CUT2')

def test_prefetch_no_event_skips_bnds(test_data_dir, test_ref_path, capsys):
    vcf_path = str(test_data_dir / "no_event_bnd.vcf")
    tra_map = prefetch_translocations(vcf_path, str(test_ref_path))
    
    assert not tra_map['insertions']
    assert not tra_map['deletions']
    
    captured = capsys.readouterr()
    assert "2 BND record(s) found without EVENT tags" in captured.out

def test_prefetch_mixed_events(test_data_dir, test_ref_path, capsys):
    vcf_path = str(test_data_dir / "mixed_events.vcf")
    tra_map = prefetch_translocations(vcf_path, str(test_ref_path))
    
    assert 'chr1' in tra_map['insertions']
    assert (500, 'COPY1') in tra_map['insertions']['chr1']
    src_chr, s_start, s_end, is_inverted, attach_after, event_id = tra_map['insertions']['chr1'][(500, 'COPY1')]
    assert event_id == 'COPY1'
    
    captured = capsys.readouterr()
    assert "2 BND record(s) found without EVENT tags" in captured.out

def test_prefetch_sequence_content(test_data_dir, test_ref_path):
    vcf_path = str(test_data_dir / "tra_copy_inter.vcf")
    tra_map = prefetch_translocations(vcf_path, str(test_ref_path))
    src_chr, s_start, s_end, is_inverted, attach_after, event_id = tra_map['insertions']['chr1'][(500, 'COPY1')]
    ref = pysam.FastaFile(str(test_ref_path))
    sequence = ref.fetch(src_chr, s_start, s_end)
    if is_inverted:
        sequence = reverse_complement(sequence)
    assert sequence == 'TGCA' * 250
