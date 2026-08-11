# Pure Mathematical Topology of Cut & Paste Translocations (`TRA_CUT`)

## 1. The Goal: 100% ID-Independent Resolution for 3rd-Party VCFs

To support **any 3rd-party VCF file** (from callers like SvABA, Manta, Gridss, or custom simulation pipelines), `inSVert` cannot rely on specific breakend ID naming conventions (such as `.H1` or `.P1`). 

We must identify the **HEAL adjacency $(H_1, H_2)$** using **pure coordinate geometry and VCF 4.2 topology invariants**.

---

## 2. The Fundamental Invariant of Cut & Paste Translocations

A Cut & Paste translocation excises a segment of length $L$ from a source chromosome (`chrom_src`) between $p_{\text{src}}$ and $p_{\text{src}} + L + 1$, and pastes it into a destination position $p_{\text{dst}}$ on `chrom_dst`.

Regardless of whether $p_{\text{dst}}$ is on a different chromosome (inter-chromosomal) or the same chromosome (intra-chromosomal), the 6 breakend records in VCF 4.2 represent 3 pairs of mates:

```text
Source Chromosome Coordinates:
--------------------------------------------------------------------------------
Base:        ... [ p_src ] [ p_src + 1 ] ... [ p_src + L ] [ p_src + L + 1 ] ...
Role:            [   H1  ] [    P2     ]     [    P3     ] [       H2      ]
Position:            x1         x2                x3               x4
--------------------------------------------------------------------------------
Relation:                x2 = x1 + 1                  x4 = x3 + 1
```

### The 4 Source Breakends ($x_1 < x_2 < x_3 < x_4$)
On `chrom_src`, 4 breakends are created around the cut site:
1. $x_1 = H_1.\text{pos} = p_{\text{src}}$ (reference base immediately preceding the cut)
2. $x_2 = P_2.\text{pos} = p_{\text{src}} + 1 = x_1 + 1$ (first base of the excised segment)
3. $x_3 = P_3.\text{pos} = p_{\text{src}} + L$ (last base of the excised segment)
4. $x_4 = H_2.\text{pos} = p_{\text{src}} + L + 1 = x_3 + 1$ (reference base immediately following the cut)

### The 2 Destination Breakends ($y_1, y_2$)
On `chrom_dst`, 2 breakends are created at the paste insertion site:
- $y_1 = P_1.\text{pos} = p_{\text{dst}}$
- $y_2 = P_4.\text{pos} = p_{\text{dst}} + 1 = y_1 + 1$

### The Mate Connections (Adjacencies)
- **Heal Pair:** $(x_1, x_4)$ — connects $H_1$ to $H_2$, bridging the source deletion.
- **Paste Start Pair:** $(y_1, x_2)$ — connects destination sink $P_1$ to cut start $P_2$.
- **Paste End Pair:** $(x_3, y_2)$ — connects cut end $P_3$ to destination sink $P_4$.

---

## 3. Mathematical Proof: Unambiguous Identification of $(H_1, H_2)$

Given **any 3 BND mate pairs** (6 records total) for a `TRA_CUT` event:

### Step 1. Separate Source and Destination Breakends
Count how many breakends lie on each chromosome:
- The **source chromosome** `chrom_src` contains **4 breakends**.
- The **destination chromosome** `chrom_dst` contains **2 breakends** (or if intra-chromosomal, all 6 breakends are on the same chromosome).

### Step 2. Topological Classification (Intra-Chromosomal Case)
When all 6 breakends are on the same chromosome, sort their coordinates:

The 4 source breakends $x_1, x_2, x_3, x_4$ satisfy:
$$x_2 = x_1 + 1 \quad \text{and} \quad x_4 = x_3 + 1$$

And the **HEAL adjacency is uniquely and unconditionally $(x_1, x_4)$** — the outermost pair of the 4 source breakends!

The remaining two breakends $(y_1, y_2)$ at $p_{\text{dst}}$ and $p_{\text{dst}} + 1$ form the destination insertion site.

---

## 4. Why the Previous Algorithm Failed on Intra-Chromosomal Events

The previous implementation attempted to find `heal_adj` by picking an arbitrary pair $(r, m)$ and checking if the other 4 breakends were split: 2 inside $(r.\text{pos}, m.\text{pos})$ and 2 outside.

For intra-chromosomal events where $p_{\text{dst}} < p_{\text{src}}$:
- The pair $(y_1, x_2)$ spanned $[p_{\text{dst}}, x_2]$.
- Inside $[y_1, x_2]$ lay $y_2$ and $x_1$ (2 breakends!).
- Outside $[y_1, x_2]$ lay $x_3$ and $x_4$ (2 breakends!).
- $x_3$ and $x_4$ satisfied $|x_4 - x_3| = 1 \le 1$.

Because $(y_1, x_2)$ accidentally satisfied the loose geometric criteria, the algorithm misidentified $(y_1, x_2)$ (spanning from $1.26\text{ MB}$ to $16.01\text{ MB}$) as the HEAL deletion pair instead of $(x_1, x_4)$ (spanning $16.01\text{ MB}$ to $16.07\text{ MB}$).

---

## 5. The Pure Coordinate Algorithm for `is_valid_tra()`

```python
def identify_tra_cut_heal_pair(adjacencies):
    """
    100% ID-independent topological identification of the HEAL pair for TRA_CUT.
    Works for any 3rd-party VCF caller (VCF 4.2).
    """
    # 1. Group breakends by chromosome
    chrom_bnds = defaultdict(list)
    for r, m in adjacencies:
        chrom_bnds[r.chrom].append(r)
        chrom_bnds[m.chrom].append(m)

    # 2. Inter-chromosomal case: source_chrom has 4 breakends, sink_chrom has 2
    for chrom, bnds in chrom_bnds.items():
        if len(bnds) == 4:
            # Sort the 4 source breakends by position: x1 < x2 < x3 < x4
            bnds.sort(key=lambda b: b.pos)
            x1, x2, x3, x4 = bnds[0], bnds[1], bnds[2], bnds[3]

            # The HEAL pair connects x1 (p_src) and x4 (p_src + L + 1)
            # Find the adjacency in adjacencies containing x1 and x4
            for r, m in adjacencies:
                if (r == x1 and m == x4) or (r == x4 and m == x1):
                    return (r, m)

    # 3. Intra-chromosomal case (all 6 breakends on same chromosome):
    # The 6 breakends consist of 4 source breakends (x1, x2, x3, x4 with x2=x1+1 and x4=x3+1)
    # and 2 sink breakends (y1, y2 with y2=y1+1).
    all_bnds = [b for pair in adjacencies for b in pair]
    all_bnds.sort(key=lambda b: b.pos)

    # Find the 4 source breakends matching the x2 = x1 + 1 and x4 = x3 + 1 pattern
    for i in range(len(all_bnds)):
        for j in range(i + 3, len(all_bnds)):
            x1, x4 = all_bnds[i], all_bnds[j]
            # Check if there exist x2 (x1.pos + 1) and x3 (x4.pos - 1) among all_bnds
            x2 = next((b for b in all_bnds if b.pos == x1.pos + 1 and b != x1), None)
            x3 = next((b for b in all_bnds if b.pos == x4.pos - 1 and b != x4), None)
            if x2 and x3 and x2.pos < x3.pos:
                # Found the source quadruplet (x1, x2, x3, x4)!
                # (x1, x4) is the HEAL pair!
                for r, m in adjacencies:
                    if (r == x1 and m == x4) or (r == x4 and m == x1):
                        return (r, m)

    return None
```

This mathematical algorithm is **100% independent of VCF record IDs**, relies purely on physical coordinate invariants of VCF 4.2 breakends, and works universally for both `inSVert` and any 3rd-party VCF caller.
