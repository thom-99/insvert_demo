# Fixing the Remaining Coordinate Clashes in inSVert

## Problem Statement

After running `inSVert simulate` followed by `inSVert insert`, the insert module reports coordinate overlap warnings (e.g., "Warning: skipped 31 variant records due to coordinate overlaps"). This should be impossible: `simulate` uses an `overlaps()` function to guarantee that no two variants share genomic coordinates. Yet `insert` encounters them anyway.

Two independent bugs cause this. One is critical (causing the vast majority of clashes), one is minor (causing rare edge-case collisions).

---

## Bug 1 (Critical): `is_valid_tra()` Misidentifies the HEAL Adjacency

### 1.1 Context: How `TRA_CUT` Works in the Pipeline

A Cut & Paste translocation (`TRA_CUT`) does two things:

1. **Cuts** a segment of length `L` out of `chrom_src` at `pos_src`.
2. **Pastes** that segment into `chrom_dst` at `pos_dst`.

In `simulate.py`, `generate_tra_bnds()` encodes this as **6 BND records** (3 mate pairs / adjacencies):

```text
Source chromosome (chrom_src):

     ... ─── pos_src ─── pos_src+1 ─── ... ─── pos_src+L ─── pos_src+L+1 ─── ...
              │  H1  │   │   P2   │             │   P3   │   │    H2     │
              └──────────────────────────────────────────────────┘
                           excised segment (length L)
```

- **Paste Start adjacency:** `(P1, P2)` — joins `chrom_dst:pos_dst` to `chrom_src:pos_src+1`
- **Paste End adjacency:** `(P3, P4)` — joins `chrom_src:pos_src+L` to `chrom_dst:pos_dst+1`
- **Heal adjacency:** `(H1, H2)` — joins `chrom_src:pos_src` to `chrom_src:pos_src+L+1`

During insertion, `prefetch_translocations()` calls `is_valid_tra()` to determine which adjacency is the **Heal pair**. The Heal pair defines the deletion span on the source chromosome. Getting this wrong means `insert.py` will delete the wrong region.

### 1.2 Why the Current Algorithm Fails (Intra-Chromosomal Case)

For **inter-chromosomal** events (`chrom_src ≠ chrom_dst`), identification is trivial: the Heal pair is the only adjacency where both breakends share the same chromosome. The two Paste pairs each span two different chromosomes.

For **intra-chromosomal** events (`chrom_src == chrom_dst`), **all three adjacencies** have both breakends on the same chromosome. The current algorithm in `is_valid_tra()` (lines 263–292 of `utils_ins.py`) tries each adjacency as a candidate Heal pair and checks:

> "If I treat this pair as Heal, do exactly 2 of the remaining 4 breakends fall inside the interval, and 2 fall outside with distance ≤ 1?"

This check is **ambiguous**. Consider a concrete example:

```text
Event: inSVert.TRA_CUT.513
pos_src = 16,019,896    (source segment starts here)
pos_dst =  1,266,852    (destination paste site)
L       =     50,852    (segment length)
```

The 6 breakends on `chr1` are:

| Breakend | Position      | Role                          |
|----------|---------------|-------------------------------|
| P1       |  1,266,852    | Paste sink start              |
| P4       |  1,266,853    | Paste sink end                |
| H1       | 16,019,896    | Heal left (source gap start)  |
| P2       | 16,019,897    | Paste source start            |
| P3       | 16,070,748    | Paste source end              |
| H2       | 16,070,749    | Heal right (source gap end)   |

Now test the **Paste Start adjacency `(P1, P2)`** as a Heal candidate:

- Interval: `[1,266,852 ... 16,019,897]` (span = 14.75 MB!)
- Inside this interval: `P4` (1,266,853) and `H1` (16,019,896) → **2 breakends inside** ✓
- Outside this interval: `P3` (16,070,748) and `H2` (16,070,749) → **2 breakends outside** ✓
- Distance between outside pair: `|16,070,749 − 16,070,748| = 1 ≤ 1` ✓

**All checks pass.** The algorithm accepts `(P1, P2)` as the Heal pair and declares a 14.75 MB deletion. This is catastrophically wrong — the real deletion is only 50.8 KB between `H1` and `H2`.

The consequence: `ref_pos` jumps by 14.75 MB, and every variant placed by `simulate` within that 14.75 MB window is flagged as a coordinate overlap and skipped. In our test run, this caused **473 false skips**.

### 1.3 The Mathematical Invariant That Solves This

The key insight is that the 6 breakend positions of a `TRA_CUT` have a rigid geometric structure. On the **source chromosome**, 4 breakends are placed at positions that satisfy:

$$x_1, \quad x_2 = x_1 + 1, \quad x_3, \quad x_4 = x_3 + 1$$

where $x_1 < x_2 < x_3 < x_4$ and $x_3 - x_2 = L - 1$.

On the **destination** side, 2 breakends sit at:

$$y_1, \quad y_2 = y_1 + 1$$

The Heal adjacency connects $x_1$ to $x_4$ (the outermost source breakends). The Paste adjacencies connect source breakends to destination breakends.

### 1.4 Proposed Algorithm: Adjacency Span Minimization

Among the 3 adjacencies of a `TRA_CUT`, each pair `(r, m)` has a **span** defined as `|r.pos − m.pos|`:

| Adjacency    | Span                                      |
|-------------|-------------------------------------------|
| Paste Start `(P1, P2)` | `|pos_dst − (pos_src + 1)|`      |
| Paste End `(P3, P4)`   | `|(pos_src + L) − (pos_dst + 1)|` |
| Heal `(H1, H2)`        | `pos_src + L + 1 − pos_src = L + 1` |

Each Paste adjacency spans from the destination site to the source site. The Heal adjacency spans exactly `L + 1` bases (the cut length plus the two flanking reference bases).

**Critical observation:** The Paste adjacencies must cross from `pos_dst` to `pos_src`. For the overlap checker in `simulate` to have allowed both positions, `pos_src` and `pos_dst` must be non-overlapping regions. Since `simulate` tracks `[pos_src, pos_src + L]` for the source and `[pos_dst, pos_dst + 1]` for the destination, these regions cannot overlap. Therefore:

- If `pos_dst < pos_src`: Paste Start spans at least `pos_src + 1 − pos_dst > L` (because `pos_dst + 1 < pos_src`).
- If `pos_dst > pos_src + L`: Paste End spans at least `pos_dst + 1 − (pos_src + L) > 1`, and Paste Start spans `pos_dst − (pos_src + 1) > L − 1`.

**In all valid configurations, the Heal adjacency has the smallest span among the intrachromosomal adjacencies.** This is because the Heal pair connects two positions that flank the excised segment (`pos_src` and `pos_src + L + 1`), while each Paste pair must bridge the gap between `pos_dst` and `pos_src`.

> [!NOTE]
> **Edge case:** If `pos_dst` falls inside `[pos_src, pos_src + L]`, the source and destination regions would overlap — but `simulate`'s `overlaps()` function prevents this. So for any VCF produced by `inSVert`, the minimum-span rule is guaranteed to be correct.
>
> For **3rd-party VCFs** where the caller may not enforce non-overlapping coordinates, the minimum-span heuristic should be combined with a validation step (checking the `x2 = x1 + 1` and `x4 = x3 + 1` quadruplet pattern on the remaining breakends).

### 1.5 Implementation

In `utils_ins.py`, replace lines 263–292 of `is_valid_tra()` (the `elif count == 3:` branch) with:

```python
elif count == 3:
    # STEP 1: Filter to intrachromosomal adjacencies only.
    # For inter-chromosomal TRA_CUT, only the Heal pair is intrachromosomal.
    # For intra-chromosomal TRA_CUT, all 3 are intrachromosomal.
    intra_adjs = [(i, r, m) for i, (r, m) in enumerate(adjacencies) if r.chrom == m.chrom]
    
    if not intra_adjs:
        return None
    
    # STEP 2: If only 1 intrachromosomal adjacency exists → it must be Heal.
    # (This is the inter-chromosomal case.)
    if len(intra_adjs) == 1:
        _, r, m = intra_adjs[0]
        heal_adj = (r, m)
    else:
        # STEP 3: Intra-chromosomal case (2 or 3 intrachromosomal adjacencies).
        # The Heal adjacency has the MINIMUM span among intrachromosomal pairs.
        # This is guaranteed because simulate enforces non-overlapping source/dest
        # regions, making Paste spans always larger than L+1.
        _, r, m = min(intra_adjs, key=lambda t: abs(t[1].pos - t[2].pos))
        heal_adj = (r, m)
    
    # VALIDATION: Verify the Heal span is consistent (span > 1, i.e., L >= 1).
    h_min = min(heal_adj[0].pos, heal_adj[1].pos)
    h_max = max(heal_adj[0].pos, heal_adj[1].pos)
    if h_max - h_min <= 1:
        return None  # Degenerate event, skip
```

The rest of the function (lines 297 onward: computing `s_start_0idx`, `s_end_0idx`, finding sink breakends, etc.) remains unchanged.

---

## Bug 2 (Minor): Source Region Tracked 1 bp Too Short in `simulate.py`

### 2.1 The Problem

In `simulate.py` line 240, the source region of a `TRA_CUT` is tracked as:

```python
utils_sim.track_sv(sv_positions, chrom_src, pos_src, pos_src + l, gt)
```

But `generate_tra_bnds()` places breakend `H2` at position `pos_src + l + 1`. This means position `pos_src + l + 1` is **not protected** by the overlap checker. A subsequent variant placed at exactly that position would pass the `overlaps()` check but would collide with the H2 breakend during insertion.

### 2.2 The Fix

In `simulate.py` line 240, extend the tracked interval by 1:

```diff
- utils_sim.track_sv(sv_positions, chrom_src, pos_src, pos_src + l, gt)
+ utils_sim.track_sv(sv_positions, chrom_src, pos_src, pos_src + l + 1, gt)
```

This ensures the full footprint of the `TRA_CUT` on the source chromosome (from `H1` at `pos_src` to `H2` at `pos_src + l + 1`) is protected from overlap.

### 2.3 Consistency Check: Destination Tracking

Line 238 tracks the destination as:

```python
utils_sim.track_sv(sv_positions, chrom_dst, pos_dst, pos_dst + 1, gt)
```

The destination breakends `P1` and `P4` sit at `pos_dst` and `pos_dst + 1`. The tracked interval `[pos_dst, pos_dst + 1]` covers both positions. ✓ No fix needed here.

### 2.4 Consistency Check: `find_valid_tra_coords` overlap check

Line 400 checks source overlaps as:

```python
overlaps(chrom_src, pos_src, pos_src + sv_length, gt, sv_positions)
```

This should also be updated to match the new tracking interval:

```diff
- overlaps(chrom_src, pos_src, pos_src + sv_length, gt, sv_positions)
+ overlaps(chrom_src, pos_src, pos_src + sv_length + 1, gt, sv_positions)
```

---

## Implementation Checklist

### Fix 1: `is_valid_tra()` Heal Pair Identification
- [ ] **File:** `inSVert/utils_ins.py`
- [ ] **Location:** Lines 263–292 (inside the `elif count == 3:` branch)
- [ ] **Action:** Replace the geometric nesting loop with the minimum-span selection algorithm described in Section 1.5 above.
- [ ] **Test:** Run `inSVert simulate` + `inSVert insert` on the 2-chromosome human subset (`test_50mb_2chroms.fa`) with `config.yaml` at ploidy 2. Verify that the "skipped variant records" warning drops from ~31 to near-zero.
- [ ] **Test:** Verify all 74 existing unit tests still pass (`pytest tests/ -v`).
- [ ] **Test:** Create a dedicated test for intra-chromosomal `TRA_CUT` where `pos_dst < pos_src` to confirm the Heal pair is correctly identified.

### Fix 2: Source Tracking Interval Off-by-One
- [ ] **File:** `inSVert/simulate.py`
- [ ] **Location:** Line 240
- [ ] **Action:** Change `pos_src + l` → `pos_src + l + 1`
- [ ] **File:** `inSVert/utils_sim.py`
- [ ] **Location:** Line 400 (inside `find_valid_tra_coords`)
- [ ] **Action:** Change `pos_src + sv_length` → `pos_src + sv_length + 1`
- [ ] **Test:** Run simulation and verify that no variant is placed at `pos_src + l + 1` for any TRA_CUT event.
