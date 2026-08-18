import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "circos_translocations.py"
SPEC = importlib.util.spec_from_file_location("circos_translocations", SCRIPT)
circos = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(circos)


def paf_line(
    query, qlen, qstart, qend, target, tlen, tstart, tend, mapq=60,
    alignment_type="P", cigar=None,
):
    block = qend - qstart
    return "\t".join(
        map(
            str,
            [query, qlen, qstart, qend, "+", target, tlen, tstart, tend, block, block, mapq],
        )
    ) + f"\ttp:A:{alignment_type}" + (f"\tcg:Z:{cigar}" if cigar else "") + "\n"


def test_read_fasta_lengths(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1 description\nAAAA\nAA\n>chr2\nTTT\n")
    assert circos.read_fasta_lengths(fasta) == {"chr1": 6, "chr2": 3}


def test_parse_and_infer_foreign_block(tmp_path):
    paf = tmp_path / "alignment.paf"
    paf.write_text(
        paf_line("sample#H1#chr1", 11000, 0, 5000, "chr1", 10000, 0, 5000)
        + paf_line("sample#H1#chr1", 11000, 5000, 6000, "chr2", 10000, 2000, 3000)
        + paf_line("sample#H1#chr1", 11000, 6000, 11000, "chr1", 10000, 5000, 10000)
    )
    alignments = circos.parse_paf(paf, min_alignment_length=500)
    events = circos.infer_translocations(alignments)
    assert events == [
        circos.Translocation(
            "sample#H1#chr1", "chr2", 2000, 3000, "chr1", 5000, 1.0, 60
        )
    ]


def test_secondary_alignment_can_explain_primary_cigar_insertion(tmp_path):
    paf = tmp_path / "alignment.paf"
    paf.write_text(
        paf_line(
            "q", 11000, 0, 11000, "chr1", 10000, 0, 10000,
            cigar="5000M1000I5000M",
        )
        + paf_line(
            "q", 11000, 5000, 6000, "chr2", 10000, 2000, 3000,
            mapq=0, alignment_type="S",
        )
    )
    events = circos.infer_translocations(
        circos.parse_paf(paf, min_alignment_length=500, min_mapq=20)
    )
    assert len(events) == 1
    assert events[0].source_chromosome == "chr2"
    assert events[0].destination_chromosome == "chr1"
    assert events[0].destination_position == 5000


def test_filter_small_contigs_removes_linked_events():
    event = circos.Translocation("q", "tiny", 0, 100, "chr1", 500, 1.0, 60)
    lengths, events, message = circos.filter_small_contigs(
        {"chr1": 10000, "tiny": 5}, [event], min_contig_fraction=0.01
    )
    assert lengths == {"chr1": 10000}
    assert events == []
    assert "Excluded 1 contig" in message
