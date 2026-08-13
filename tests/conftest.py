import pytest
from pathlib import Path

@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data" 

@pytest.fixture
def test_ref_path(test_data_dir):
    return test_data_dir / "test_ref.fa"

class MockInfo:
    def __init__(self, data):
        self._data = data
    def get(self, key, default=None):
        return self._data.get(key, default)

class MockBND:
    def __init__(self, chrom, pos, id, alt_string, mateid, event=None):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.alts = (alt_string,) if alt_string else None
        self.info = MockInfo({'SVTYPE': 'BND', 'MATEID': (mateid,), 'EVENT': event})

@pytest.fixture
def make_bnd():
    return MockBND
