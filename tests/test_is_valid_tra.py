import pytest
from inSVert.utils_ins import is_valid_tra

# --- TRA_COPY (2 adjacencies) ---

def test_tra_copy_inter_valid(make_bnd):
    P1 = make_bnd('chr1', 500, 'COPY1.P1', 'T[chr2:1001[', 'COPY1.P2', 'COPY1')
    P2 = make_bnd('chr2', 1001, 'COPY1.P2', ']chr1:500]T', 'COPY1.P1', 'COPY1')
    P3 = make_bnd('chr2', 2000, 'COPY1.P3', 'A[chr1:501[', 'COPY1.P4', 'COPY1')
    P4 = make_bnd('chr1', 501, 'COPY1.P4', ']chr2:2000]A', 'COPY1.P3', 'COPY1')
    adjacencies = [(P1, P2), (P3, P4)]
    assert is_valid_tra('COPY1', adjacencies) == ('TRA_COPY', 'chr2', (1000, 2000), 'chr1', 500, None, False, True)

def test_tra_copy_intra_valid(make_bnd):
    P1 = make_bnd('chr1', 500, 'COPY2.P1', 'T[chr1:3001[', 'COPY2.P2', 'COPY2')
    P2 = make_bnd('chr1', 3001, 'COPY2.P2', ']chr1:500]A', 'COPY2.P1', 'COPY2')
    P3 = make_bnd('chr1', 4000, 'COPY2.P3', 'T[chr1:501[', 'COPY2.P4', 'COPY2')
    P4 = make_bnd('chr1', 501, 'COPY2.P4', ']chr1:4000]A', 'COPY2.P3', 'COPY2')
    adjacencies = [(P1, P2), (P3, P4)]
    assert is_valid_tra('COPY2', adjacencies) == ('TRA_COPY', 'chr1', (3000, 4000), 'chr1', 500, None, False, True)

def test_tra_copy_both_sides_large_distance(make_bnd):
    P1 = make_bnd('chr1', 100, 'ID1', 'A[chr2:500[', 'ID2')
    P2 = make_bnd('chr2', 500, 'ID2', ']chr1:100]A', 'ID1')
    P3 = make_bnd('chr2', 1500, 'ID3', 'A[chr1:101[', 'ID4')
    P4 = make_bnd('chr1', 101, 'ID4', ']chr2:1500]A', 'ID3')
    adjacencies = [(P1, P2), (P3, P4)]
    assert is_valid_tra('TEST', adjacencies) == ('TRA_COPY', 'chr2', (499, 1500), 'chr1', 100, None, False, True)

def test_tra_copy_single_adjacency_returns_none(make_bnd):
    P1 = make_bnd('chr1', 500, 'COPY1.P1', 'T[chr2:1001[', 'COPY1.P2')
    P2 = make_bnd('chr2', 1001, 'COPY1.P2', ']chr1:500]T', 'COPY1.P1')
    assert is_valid_tra('TEST', [(P1, P2)]) is None

def test_tra_copy_empty_adjacencies():
    assert is_valid_tra('TEST', []) is None

# --- TRA_CUT (3 adjacencies) ---

def test_tra_cut_inter_valid(make_bnd):
    P1 = make_bnd('chr1', 500, 'CUT1.P1', 'T[chr2:1001[', 'CUT1.P2', 'CUT1')
    P2 = make_bnd('chr2', 1001, 'CUT1.P2', ']chr1:500]T', 'CUT1.P1', 'CUT1')
    P3 = make_bnd('chr2', 2000, 'CUT1.P3', 'A[chr1:501[', 'CUT1.P4', 'CUT1')
    P4 = make_bnd('chr1', 501, 'CUT1.P4', ']chr2:2000]A', 'CUT1.P3', 'CUT1')
    H1 = make_bnd('chr2', 1000, 'CUT1.H1', 'A[chr2:2001[', 'CUT1.H2', 'CUT1')
    H2 = make_bnd('chr2', 2001, 'CUT1.H2', ']chr2:1000]T', 'CUT1.H1', 'CUT1')
    adjacencies = [(P1, P2), (P3, P4), (H1, H2)]
    assert is_valid_tra('CUT1', adjacencies) == ('TRA_CUT', 'chr2', (1000, 2000), 'chr1', 500, 1000, False, True)

def test_tra_cut_intra_valid(make_bnd):
    P1 = make_bnd('chr1', 7000, 'CUT2.P1', 'T[chr1:3001[', 'CUT2.P2', 'CUT2')
    P2 = make_bnd('chr1', 3001, 'CUT2.P2', ']chr1:7000]A', 'CUT2.P1', 'CUT2')
    P3 = make_bnd('chr1', 4000, 'CUT2.P3', 'T[chr1:7001[', 'CUT2.P4', 'CUT2')
    P4 = make_bnd('chr1', 7001, 'CUT2.P4', ']chr1:4000]A', 'CUT2.P3', 'CUT2')
    H1 = make_bnd('chr1', 3000, 'CUT2.H1', 'T[chr1:4001[', 'CUT2.H2', 'CUT2')
    H2 = make_bnd('chr1', 4001, 'CUT2.H2', ']chr1:3000]A', 'CUT2.H1', 'CUT2')
    adjacencies = [(H1, H2), (P1, P2), (P3, P4)]
    assert is_valid_tra('CUT2', adjacencies) == ('TRA_CUT', 'chr1', (3000, 4000), 'chr1', 7000, 3000, False, True)

def test_tra_cut_no_heal_returns_none(make_bnd):
    P1 = make_bnd('chr1', 500, 'P1', 'T[chr2:1001[', 'P2')
    P2 = make_bnd('chr2', 1001, 'P2', ']chr1:500]T', 'P1')
    P3 = make_bnd('chr2', 2000, 'P3', 'A[chr1:501[', 'P4')
    P4 = make_bnd('chr1', 501, 'P4', ']chr2:2000]A', 'P3')
    P5 = make_bnd('chr1', 600, 'P5', 'C[chr3:100[', 'P6')
    P6 = make_bnd('chr3', 100, 'P6', ']chr1:600]C', 'P5')
    adjacencies = [(P1, P2), (P3, P4), (P5, P6)]
    assert is_valid_tra('TEST', adjacencies) is None

# --- Return type validation ---
def test_return_type_tra_copy(make_bnd):
    P1 = make_bnd('chr1', 500, 'COPY1.P1', 'T[chr2:1001[', 'COPY1.P2')
    P2 = make_bnd('chr2', 1001, 'COPY1.P2', ']chr1:500]T', 'COPY1.P1')
    P3 = make_bnd('chr2', 2000, 'COPY1.P3', 'A[chr1:501[', 'COPY1.P4')
    P4 = make_bnd('chr1', 501, 'COPY1.P4', ']chr2:2000]A', 'COPY1.P3')
    res = is_valid_tra('TEST', [(P1, P2), (P3, P4)])
    assert len(res) == 8
    assert isinstance(res[0], str)
    assert isinstance(res[1], str)
    assert isinstance(res[2], tuple)
    assert len(res[2]) == 2
    assert isinstance(res[3], str)
    assert isinstance(res[4], int)
    assert res[5] is None
    assert isinstance(res[6], bool)
    assert isinstance(res[7], bool)

def test_return_type_tra_cut(make_bnd):
    P1 = make_bnd('chr1', 500, 'CUT1.P1', 'T[chr2:1001[', 'CUT1.P2')
    P2 = make_bnd('chr2', 1001, 'CUT1.P2', ']chr1:500]T', 'CUT1.P1')
    P3 = make_bnd('chr2', 2000, 'CUT1.P3', 'A[chr1:501[', 'CUT1.P4')
    P4 = make_bnd('chr1', 501, 'CUT1.P4', ']chr2:2000]A', 'CUT1.P3')
    H1 = make_bnd('chr2', 1000, 'CUT1.H1', 'A[chr2:2001[', 'CUT1.H2')
    H2 = make_bnd('chr2', 2001, 'CUT1.H2', ']chr2:1000]T', 'CUT1.H1')
    res = is_valid_tra('TEST', [(P1, P2), (P3, P4), (H1, H2)])
    assert len(res) == 8
    assert res[0] == 'TRA_CUT'
    assert isinstance(res[5], int)

def test_is_inverted_false_for_forward_brackets(make_bnd):
    P1 = make_bnd('chr1', 500, 'COPY1.P1', 'T[chr2:1001[', 'COPY1.P2')
    P2 = make_bnd('chr2', 1001, 'COPY1.P2', ']chr1:500]T', 'COPY1.P1')
    P3 = make_bnd('chr2', 2000, 'COPY1.P3', 'A[chr1:501[', 'COPY1.P4')
    P4 = make_bnd('chr1', 501, 'COPY1.P4', ']chr2:2000]A', 'COPY1.P3')
    res = is_valid_tra('TEST', [(P1, P2), (P3, P4)])
    assert res[6] is False

# --- Edge cases ---
def test_four_adjacencies_returns_none(make_bnd):
    P1 = make_bnd('chr1', 1, 'P1', 'A[chr2:2[', 'P2')
    P2 = make_bnd('chr2', 2, 'P2', ']chr1:1]A', 'P1')
    adjacencies = [(P1, P2)] * 4
    assert is_valid_tra('TEST', adjacencies) is None
