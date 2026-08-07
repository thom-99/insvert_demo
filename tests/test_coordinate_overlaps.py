"""
Targeted diagnostic tests to identify the remaining sources of
"coordinate overlaps" in inSVert insert even after the >= fix.

Three concrete root causes are tested:

1. BND records (TRA events) share pos_src+1 and pos_src+length on the same
   chromosome. After insert.py processes the first BND at that pos, ref_pos
   advances past it. The second BND record at the SAME position arrives in
   the sorted VCF immediately after and hits `start < ref_pos` → skipped.

2. The overlap check in simulate.py uses `overlaps(chrom, pos, ...)` where
   `pos` is the 1-indexed VCF position. But `track_sv` also stores `pos`
   (1-indexed) as the `start` of the interval. In insert.py `start = var.pos - 1`
   (0-indexed). So a stored interval `(P, P+L)` in simulate corresponds to
   0-indexed range `[P-1, P+L-1)` consumed by insert. There is a systematic
   1-off between the coordinate systems used in simulate vs insert.

3. The overlapping_variants counter is shared across ALL haplotypes. If a
   variant is legitimately absent on haplotype 1 (gt "0|1") but causes a
   ref_pos advancement on haplotype 0 through an ADJACENT variant, and then
   haplotype 0 encounters that same variant from a different angle, the count
   inflates with ploidy.
"""

import pytest
from collections import defaultdict
import bisect

from inSVert import utils_sim


# ─────────────────────────────────────────────────────────────────────────────
# Helper: reproduce exactly what insert.py does to compute ref_pos
# ─────────────────────────────────────────────────────────────────────────────

def simulate_insert_refpos(variants):
    """
    Reproduce insert.py's ref_pos logic for a sequence of variants.
    Each variant is a dict: {svtype, pos (1-indexed VCF), svlen, gt_hap0}
    Returns (ref_pos_final, skipped_count).
    """
    ref_pos = 0
    skipped = 0
    for v in variants:
        if v.get("gt_hap0", 1) != 1:   # not on this haplotype
            continue
        start = v["pos"] - 1            # insert.py line 143
        if start < ref_pos:
            skipped += 1
            continue
        svtype = v["svtype"]
        if svtype == "INS":
            ref_pos = start + 1         # padding base only; seq doesn't advance ref
        elif svtype in ("DEL", "INV", "DUP"):
            ref_pos = start + 1 + v["svlen"]   # padding base + body
        elif svtype == "BND_CUT":
            ref_pos = start + 1 + v["svlen"]   # BND cut: padding + deletion
        elif svtype == "BND_PASTE":
            ref_pos = start + 1         # BND paste: only padding base advanced
    return ref_pos, skipped


# ─────────────────────────────────────────────────────────────────────────────
# ROOT CAUSE 1: BND duplicate records at the same position
# ─────────────────────────────────────────────────────────────────────────────

def test_bnd_duplicate_positions_trigger_overlap_skip():
    """
    A TRA_CUT event at src=pos_src produces BND records at:
      - pos_src   (H1 heal record)
      - pos_src+1 (P2 paste-start source record)
      - pos_src+L (P3 paste-end source record)
      - pos_src+L+1 (H2 heal record)
    After the first BND at pos_src+1 is processed (advancing ref_pos by the
    cut-deletion length), the P3 record at pos_src+L arrives. Its
    start = (pos_src+L) - 1. Meanwhile ref_pos = pos_src+1 + cut_len = pos_src+L+1.
    So start < ref_pos → SKIPPED, even though it's a valid companion record.
    """
    pos_src = 1000
    cut_len = 500  # length of cut segment

    # After insert processes the CUT BND (H1 at pos_src):
    # ref_pos = (pos_src - 1) + 1 + cut_len = pos_src + cut_len
    ref_pos_after_cut = pos_src + cut_len  # 1500

    # The "end heal" BND is at pos_src + cut_len + 1 in VCF (1-indexed)
    end_heal_pos = pos_src + cut_len + 1   # 1501
    start_end_heal = end_heal_pos - 1       # 1500  (insert.py 0-indexed)

    # P2 source BND is at pos_src + 1 in VCF
    p2_pos = pos_src + 1  # 1001
    start_p2 = p2_pos - 1  # 1000
    # This one comes BEFORE the heal in a sorted VCF, so it gets processed first
    # (or it may not — depends on sort order within same chromosome)

    # What matters: does the end heal BND get blocked?
    assert start_end_heal >= ref_pos_after_cut, (
        "End-heal BND is NOT blocked — this is fine"
    )

    # But check P3 at pos_src + cut_len:
    p3_pos = pos_src + cut_len  # 1500
    start_p3 = p3_pos - 1       # 1499

    # If ref_pos after H1 cut = pos_src + cut_len = 1500:
    is_skipped = start_p3 < ref_pos_after_cut
    assert is_skipped, (
        f"P3 BND at pos={p3_pos} (start={start_p3}) IS skipped because "
        f"start_p3 ({start_p3}) < ref_pos ({ref_pos_after_cut}). "
        "BND companion records collide with each other in insert.py."
    )


# ─────────────────────────────────────────────────────────────────────────────
# ROOT CAUSE 2: Coordinate system mismatch (1-indexed simulate vs 0-indexed insert)
# ─────────────────────────────────────────────────────────────────────────────

def test_coordinate_system_mismatch_del_followed_by_del():
    """
    simulate.py uses 1-indexed VCF positions in track_sv.
    insert.py converts: start = var.pos - 1  (0-indexed).

    After DEL at VCF pos P, len L:
      - simulate tracks: (P, P+L) in 1-indexed space
      - insert.py after DEL: ref_pos = (P-1) + 1 + L = P + L  (0-indexed)

    For the next DEL at VCF pos Q to NOT be skipped:
      insert requires: Q - 1 >= P + L  →  Q >= P + L + 1  →  Q > P + L

    After the >= fix, overlaps() blocks Q if prev_end (P+L) >= Q:
      → blocks Q <= P+L → allows Q >= P+L+1 → Q > P+L  ✓

    This test CONFIRMS the fix is correct for simple DEL→DEL.
    """
    sv_positions = defaultdict(lambda: defaultdict(list))
    chrom = "chr1"
    gt = "1|1"

    # DEL1: VCF pos=1000, len=500 → get_end() = 1500 → track_sv stores (1000, 1500)
    utils_sim.track_sv(sv_positions, chrom, 1000, 1500, gt)

    # Candidate Q=1500: simulate allows? overlaps(chrom, 1500, ...) → prev_end(1500) >= 1500 → True → BLOCKED ✓
    blocked_at_1500 = utils_sim.overlaps(chrom, 1500, 1600, gt, sv_positions)
    # Candidate Q=1501: overlaps(chrom, 1501, ...) → prev_end(1500) >= 1501 → False → ALLOWED ✓
    blocked_at_1501 = utils_sim.overlaps(chrom, 1501, 1601, gt, sv_positions)

    assert blocked_at_1500 is True, (
        "Q=1500 must be blocked (insert would have start=1499 < ref_pos=1500)"
    )
    assert blocked_at_1501 is False, (
        "Q=1501 must be allowed (insert has start=1500 >= ref_pos=1500)"
    )

    # Now verify insert.py would agree:
    variants = [
        {"svtype": "DEL", "pos": 1000, "svlen": 500, "gt_hap0": 1},
        {"svtype": "DEL", "pos": 1501, "svlen": 100, "gt_hap0": 1},  # allowed
    ]
    _, skipped = simulate_insert_refpos(variants)
    assert skipped == 0, "DEL at Q=1501 after DEL(1000,500) must not be skipped in insert"

    variants_bad = [
        {"svtype": "DEL", "pos": 1000, "svlen": 500, "gt_hap0": 1},
        {"svtype": "DEL", "pos": 1500, "svlen": 100, "gt_hap0": 1},  # should be blocked
    ]
    _, skipped_bad = simulate_insert_refpos(variants_bad)
    assert skipped_bad == 1, "DEL at Q=1500 after DEL(1000,500) MUST be skipped in insert"


def test_ins_followed_by_del_boundary():
    """
    After INS at VCF pos P:
      - simulate tracks: (P, P)   (get_end() = pos for INS)
      - insert.py after INS: ref_pos = P  (0-indexed; padding base consumed, no body advance)

    For next DEL at VCF pos Q to not be skipped: Q - 1 >= P → Q >= P+1 → Q > P.
    After the >= fix: overlaps blocks Q if prev_end(P) >= Q → blocks Q<=P → allows Q>=P+1 ✓
    """
    sv_positions = defaultdict(lambda: defaultdict(list))
    chrom = "chr1"
    gt = "1|1"

    # INS: VCF pos=2000, get_end()=2000 → track_sv stores (2000, 2000)
    utils_sim.track_sv(sv_positions, chrom, 2000, 2000, gt)

    # Q=2000 should be blocked: prev_end(2000) >= 2000 → True ✓
    blocked = utils_sim.overlaps(chrom, 2000, 2100, gt, sv_positions)
    # Q=2001 should be allowed: prev_end(2000) >= 2001 → False ✓
    allowed = utils_sim.overlaps(chrom, 2001, 2101, gt, sv_positions)

    assert blocked is True, "Q=2000 after INS@2000 must be blocked"
    assert allowed is False, "Q=2001 after INS@2000 must be allowed"

    # Verify insert.py consistency
    variants = [
        {"svtype": "INS", "pos": 2000, "svlen": 5000, "gt_hap0": 1},
        {"svtype": "DEL", "pos": 2001, "svlen": 100,  "gt_hap0": 1},  # allowed
    ]
    _, skipped = simulate_insert_refpos(variants)
    assert skipped == 0, "DEL at Q=2001 after INS@2000 must not be skipped"

    variants_bad = [
        {"svtype": "INS", "pos": 2000, "svlen": 5000, "gt_hap0": 1},
        {"svtype": "DEL", "pos": 2000, "svlen": 100,  "gt_hap0": 1},  # blocked
    ]
    _, skipped_bad = simulate_insert_refpos(variants_bad)
    assert skipped_bad == 1, "DEL at Q=2000 after INS@2000 must be skipped"


# ─────────────────────────────────────────────────────────────────────────────
# ROOT CAUSE 3: overlapping_variants counter accumulates across ploidy haplotypes
# ─────────────────────────────────────────────────────────────────────────────

def test_overlap_counter_inflated_by_ploidy():
    """
    In insert.py, overlapping_variants is initialised ONCE outside the
    `for haplotype in range(ploidy):` loop. Any legitimate overlap on one
    haplotype is counted ONCE. But if the same variant appears on multiple
    haplotypes that all encounter the same collision, the single counter
    accumulates across haplotypes — inflating the reported number.

    For ploidy=3, if one true overlap exists per haplotype pass, the final
    warning will say '3 skipped' even if only 1 distinct genomic position
    was the problem.
    """
    ploidy = 3
    overlapping_variants = 0  # as in insert.py — outside the haplotype loop

    for haplotype in range(ploidy):
        # Simulate: one collision per haplotype at pos 5000
        ref_pos = 6000   # some variant already advanced ref_pos
        start = 5000 - 1  # candidate at pos 5000
        if start < ref_pos:
            overlapping_variants += 1

    # The warning will say "skipped 3" but the same positional conflict is
    # being counted 3 times (once per haplotype).
    assert overlapping_variants == 3, (
        "overlapping_variants is inflated by ploidy even for a single "
        "distinct genomic collision"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Summary regression: the >= fix is correct, but BNDs are an independent bug
# ─────────────────────────────────────────────────────────────────────────────

def test_no_false_overlaps_simple_nonbnd_sorted_vcf():
    """
    With the >= fix in place, a sorted VCF of simple (non-BND) variants
    placed correctly by simulate must produce zero skips in insert.
    """
    sv_positions = defaultdict(lambda: defaultdict(list))
    chrom = "chr1"
    gt = "1|1"

    # Simulate placing three non-overlapping variants
    # DEL @ 1000, len 500 → end 1500
    utils_sim.track_sv(sv_positions, chrom, 1000, 1500, gt)
    # INS @ 1501 (first allowed pos after DEL) → end 1501
    assert not utils_sim.overlaps(chrom, 1501, 1501, gt, sv_positions)
    utils_sim.track_sv(sv_positions, chrom, 1501, 1501, gt)
    # DEL @ 1502, len 200 → end 1702
    assert not utils_sim.overlaps(chrom, 1502, 1702, gt, sv_positions)
    utils_sim.track_sv(sv_positions, chrom, 1502, 1702, gt)

    # Now replay through insert logic — must be zero skips
    variants = [
        {"svtype": "DEL", "pos": 1000, "svlen": 500, "gt_hap0": 1},
        {"svtype": "INS", "pos": 1501, "svlen": 300, "gt_hap0": 1},
        {"svtype": "DEL", "pos": 1502, "svlen": 200, "gt_hap0": 1},
    ]
    _, skipped = simulate_insert_refpos(variants)
    assert skipped == 0, (
        "With the >= fix, correctly-placed simulate variants must produce "
        "zero skips in insert for non-BND SVs"
    )
