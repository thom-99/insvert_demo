import pytest
import io
from inSVert.utils_ins import reverse_complement, generate_seq, apply_insertion, apply_deletion

def test_reverse_complement_basic():
    assert reverse_complement('ACGT') == 'ACGT'

def test_reverse_complement_single():
    assert reverse_complement('A') == 'T'

def test_reverse_complement_polyA():
    assert reverse_complement('AAAA') == 'TTTT'

def test_reverse_complement_with_N():
    assert reverse_complement('ACNGT') == 'ACNGT'

def test_reverse_complement_lowercase():
    assert reverse_complement('acgt') == 'acgt'

def test_reverse_complement_empty():
    assert reverse_complement('') == ''

def test_reverse_complement_palindrome():
    assert reverse_complement('AATT') == 'AATT'

def test_generate_seq_length():
    seq = generate_seq(100, 0.5)
    assert len(seq) == 100

def test_generate_seq_zero_length():
    # Wait, original generate_seq raises ValueError for 0?
    # Let's check: if length < 0: raise ValueError
    # If length == 0 it should return '' 
    seq = generate_seq(0, 0.5)
    assert seq == ''

def test_generate_seq_negative_length():
    with pytest.raises(ValueError):
        generate_seq(-1, 0.5)

def test_generate_seq_only_dna_bases():
    seq = generate_seq(1000, 0.5)
    assert set(seq).issubset(set('ACGT'))

def test_apply_insertion():
    buffer = io.StringIO()
    apply_insertion(buffer, 'ACGT')
    assert buffer.getvalue() == 'ACGT'

def test_apply_deletion():
    new_pos = apply_deletion(100, 50)
    assert new_pos == 150

def test_apply_deletion_negative_length():
    new_pos = apply_deletion(100, -50)
    assert new_pos == 150
