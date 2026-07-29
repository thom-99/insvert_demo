import pytest
from inSVert.utils_ins import parse_bnd_orientation

@pytest.mark.parametrize("alt_string, expected", [
    # 1. All 4 standard VCF 4.2 patterns with single-base sequence
    ('A[chr2:1000[', (False, True)),
    ('A]chr2:1000]', (True, True)),
    (']chr2:1000]A', (False, False)),
    ('[chr2:1000[A', (True, False)),

    # 2. Multi-base sequences
    ('ACGT[chr2:1000[', (False, True)),
    (']chr2:1000]ACGT', (False, False)),

    # 3. Complex chromosome names
    ('A[chr2_random:1000[', (False, True)),
    ('A[chrUn_gl000220:5000[', (False, True)),

    # 4. Malformed inputs
    ('A[chr2:1000]', (None, None)),
    ('A]chr2:1000[', (None, None)),
    ('<DEL>', (None, None)),
    ('', (None, None)),
    ('ACGT', (None, None)),
    ('A[chr2[', (None, None)),

    # 5. Edge cases
    (' A[chr2:1000[ ', (False, True)), # Does not strip internally
    ('A[chr1:999999999[', (False, True)),
])
def test_parse_bnd_orientation(alt_string, expected):
    assert parse_bnd_orientation(alt_string) == expected
