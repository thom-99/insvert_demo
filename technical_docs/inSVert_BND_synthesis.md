# Comprehensive Guide: How inSVert Writes and Expects Breakends (BNDs)

This document synthesizes how **inSVert** generates (writes) and parses (expects to receive) translocations and dispersed duplications represented in VCF 4.2 Breakend (`BND`) format, based on the codebase implementation in [`inSVert/utils_sim.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/utils_sim.py) and [`inSVert/utils_ins.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/utils_ins.py).

---

## 1. Conceptual Overview & Classification

inSVert models complex structural variants using **VCF 4.2 Breakend (`BND`)** specifications. It categorizes translocations into two distinct biological event types:

1. **Dispersed Duplication (`TRA_COPY` / Interspersed Duplication / Copy & Paste)**:
   - A segment from a source chromosome region is duplicated and inserted into a target sink location.
   - Represented by **2 adjacencies (4 `BND` records)**.
2. **Translocation (`TRA_CUT` / Transposition / Cut & Paste)**:
   - A segment is excised from a source chromosome region and pasted into a target sink location. The genomic gap left on the source chromosome is repaired ("healed").
   - Represented by **3 adjacencies (6 `BND` records)**: 
     - 2 adjacencies (4 `BND`s) for the paste junction at the destination.
     - 1 HEAL adjacency (2 `BND`s) bridging the cut site on the source chromosome.

---

## 2. Mandatory VCF Attributes Required by inSVert

For inSVert to successfully group, validate, and parse BND records, input VCFs must satisfy:

* **`SVTYPE=BND`**: Required in the `INFO` column for all breakend records.
* **`EVENT=<event_id>`**: **Crucial.** All breakends belonging to the same biological event must share the exact same `EVENT` ID tag in `INFO`.
* **`MATEID=<mate_id>`**: Specifies the reciprocal breakend ID for each pair in an adjacency.
* **Phased Genotypes (`GT`)**: Phased genotypes (e.g., `0|1` or `1|0`) ensure deterministic haplotype placement. (Unphased `0/1` records sharing an `EVENT` ID are assigned consistently across all mates).

---

## 3. ALT Bracket Notation & Strand Orientations

inSVert parses ALT strings using regular expressions ([`parse_bnd_orientation`](file:///home/thomas/Desktop/insvert_demo/inSVert/utils_ins.py#L428-L466)) to determine insertion orientation and attachment position:

| ALT Pattern | Bracket Type | Meaning | Insertion Orientation |
| :--- | :--- | :--- | :--- |
| `t[p[` | Forward | Joined sequence extending right of `p` attaches *after* base `t` | **Forward** |
| `]p]t` | Forward | Joined sequence extending left of `p` attaches *before* base `t` | **Forward** |
| `t]p]` | Reverse | Reverse-complemented sequence extending left of `p` attaches *after* base `t` | **Inverted (Reverse Complement)** |
| `[p[t` | Reverse | Reverse-complemented sequence extending right of `p` attaches *before* base `t` | **Inverted (Reverse Complement)** |

---

## 4. How inSVert Writes BNDs (`simulate.py` / `utils_sim.py`)

When generating simulated variants ([`generate_tra_bnds`](file:///home/thomas/Desktop/insvert_demo/inSVert/utils_sim.py#L417-L493)):

1. **PASTE Breakends (P1, P2, P3, P4)**:
   - `P1` (at `chrom_dst:pos_dst`) links to `P2` (at `chrom_src:pos_src+1`).
   - `P3` (at `chrom_src:pos_src+length`) links to `P4` (at `chrom_dst:pos_dst+1`).
   - For **forward** orientation, bracket formats `t[p[` and `]p]t` are constructed.
   - For **reverse** (inverted) orientation, bracket formats `t]p]` and `[p[t` are constructed.
2. **HEAL Breakends (H1, H2)** *(Used only for `TRA_CUT`)*:
   - `H1` (at `chrom_src:pos_src`) links to `H2` (at `chrom_src:pos_src+length+1`).
   - Re-joins the genomic reference flanking the excised segment (always forward orientation).

---

## 5. How inSVert Expects to Receive & Parse BNDs (`utils_ins.py`)

During genome reconstruction ([`prefetch_translocations`](file:///home/thomas/Desktop/insvert_demo/inSVert/utils_ins.py#L275-L346) & [`is_valid_tra`](file:///home/thomas/Desktop/insvert_demo/inSVert/utils_ins.py#L93-L272)):

1. **Group by `EVENT`**: Reads all VCF records with `SVTYPE=BND` and groups them by `EVENT`.
2. **Topology & Distance Check**:
   - **2 Adjacencies (`TRA_COPY`)**:
     - Evaluates the distance between paired breakends on both sides.
     - The side with distance $> 1\text{ bp}$ is recognized as the **source segment** (`chrom_src:[start, end]`).
     - The side with distance $\le 1\text{ bp}$ is recognized as the **target sink site** (`chrom_dst:pos_dst`).
     - Registers an **insertion job** at `chrom_dst:pos_dst` to copy `ref.fetch(chrom_src, start-1, end)`.
   - **3 Adjacencies (`TRA_CUT`)**:
     - Identifies the **HEAL adjacency `(H1, H2)`** on `chrom_src`, which encloses the 2 source paste breakends inside the interval `(h_min, h_max)`.
     - Registers a **deletion job** at `h_min` of length `(h_max - 1) - h_min` (skipping the excised reference sequence).
     - Registers an **insertion job** at the sink position `pos_dst` outside the HEAL interval.

---

## 6. Concrete VCF Examples

### Example A: Dispersed Duplication (`TRA_COPY` — 4 BND lines)
*Copies `chr2:1001-2000` and inserts it into `chr1:500` in forward orientation.*

```vcf
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated with breakend">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	500	DUP_EVENT1.P1	T	T[chr2:1001[	.	PASS	SVTYPE=BND;EVENT=DUP_EVENT1;MATEID=DUP_EVENT1.P2	GT	0|1
chr2	1001	DUP_EVENT1.P2	A	]chr1:500]A	.	PASS	SVTYPE=BND;EVENT=DUP_EVENT1;MATEID=DUP_EVENT1.P1	GT	0|1
chr2	2000	DUP_EVENT1.P3	C	C[chr1:501[	.	PASS	SVTYPE=BND;EVENT=DUP_EVENT1;MATEID=DUP_EVENT1.P4	GT	0|1
chr1	501	DUP_EVENT1.P4	G	]chr2:2000]G	.	PASS	SVTYPE=BND;EVENT=DUP_EVENT1;MATEID=DUP_EVENT1.P3	GT	0|1
```

---

### Example B: Cut & Paste Translocation (`TRA_CUT` — 6 BND lines)
*Excises `chr2:1001-2000` (heals `chr2:1000` to `chr2:2001`) and pastes it into `chr1:500`.*

```vcf
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated with breakend">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	500	CUT_EVENT1.P1	T	T[chr2:1001[	.	PASS	SVTYPE=BND;EVENT=CUT_EVENT1;MATEID=CUT_EVENT1.P2	GT	0|1
chr2	1001	CUT_EVENT1.P2	A	]chr1:500]A	.	PASS	SVTYPE=BND;EVENT=CUT_EVENT1;MATEID=CUT_EVENT1.P1	GT	0|1
chr2	2000	CUT_EVENT1.P3	C	C[chr1:501[	.	PASS	SVTYPE=BND;EVENT=CUT_EVENT1;MATEID=CUT_EVENT1.P4	GT	0|1
chr1	501	CUT_EVENT1.P4	G	]chr2:2000]G	.	PASS	SVTYPE=BND;EVENT=CUT_EVENT1;MATEID=CUT_EVENT1.P3	GT	0|1
chr2	1000	CUT_EVENT1.H1	A	A[chr2:2001[	.	PASS	SVTYPE=BND;EVENT=CUT_EVENT1;MATEID=CUT_EVENT1.H2	GT	0|1
chr2	2001	CUT_EVENT1.H2	T	]chr2:1000]T	.	PASS	SVTYPE=BND;EVENT=CUT_EVENT1;MATEID=CUT_EVENT1.H1	GT	0|1
```

---

## 7. Summary Reference Table

| Event Type | Adjacency Count | Total BND Records | Source Chromosome Action | Target Chromosome Action |
| :--- | :---: | :---: | :--- | :--- |
| **`TRA_COPY`** | 2 | 4 | Remains unmodified | Sequence copied & inserted at `pos_dst` |
| **`TRA_CUT`** | 3 | 6 | Sequence deleted across HEAL interval | Sequence pasted at `pos_dst` |
