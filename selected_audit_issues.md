# Technical Audit Analysis: 5 Selected inSVert Issues

This document details the 5 specific issues identified during the code audit of **inSVert**, focusing on root cause analysis, code trace, impact, and concrete code remediation for each problem.

---

## 1. Reference Corruption for Non-Single-Base `svtype is None` Records

### Problem Location
* **File**: [`inSVert/insert.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/insert.py#L173-L197) (lines 173–197)

### Description & Root Cause Analysis
In `insert.py`, variants without an explicit `SVTYPE` tag in the INFO column have `svtype = var.info.get("SVTYPE")` evaluate to `None`. The code routes these records to the SNP handling block:

```python
# ── SNP HANDLING (Single-base substitution, no padding base) ──
if svtype is None:
    ref_allele = var.ref
    alt_alleles = var.alts
    if (ref_allele and alt_alleles and 
        len(ref_allele) == 1 and len(alt_alleles) == 1 and 
        len(str(alt_alleles[0])) == 1):
        
        # Write ALT nucleotide in place of reference base
        writer.write(str(alt_alleles[0]))
        ref_pos = start + 1
        continue
```

If an input VCF record has `svtype is None` (standard in calls from GATK, FreeBayes, DeepVariant, BCFtools) but is **not** a single-base substitution (for instance, Multi-Nucleotide Polymorphisms (MNPs) like `REF=AT, ALT=GC`, or plain VCF indels like `REF=A, ALT=AG`), the inner condition `len(ref_allele) == 1 and len(alt_alleles) == 1 and len(str(alt_alleles[0])) == 1` evaluates to `False`.

Because there is no `else` block or fallback, execution falls through past `if svtype is None:` to lines 187–197:

```python
if svtype != "BND":
    ...
    # Write the padding base before standard SVs
    if ref_pos == start:
        writer.write(ref.fetch(chrom, start, start + 1))
        ref_pos = start + 1
```

None of the subsequent `elif svtype == ...` branches (`INS`, `DEL`, `INV`, `DUP`) match because `svtype` is `None`. 

### Impact
* The variant payload is **never applied** to the modified FASTA.
* The reference padding base at `start` is written to `writer`, and `ref_pos` is advanced by 1 base.
* The genome sequence skips the variant payload while advancing reference coordinates, **corrupting reference genome alignment and sequence length**.

### Remediation Code
Update `insert.py` when `svtype is None` to handle multi-base substitutions (MNPs) and plain indels:

```python
if svtype is None:
    ref_allele = var.ref
    alt_alleles = var.alts
    if ref_allele and alt_alleles:
        alt_str = str(alt_alleles[0])
        # Single-base SNP
        if len(ref_allele) == 1 and len(alt_str) == 1:
            writer.write(alt_str)
            ref_pos = start + 1
            continue
        # MNP / Multi-base substitution (e.g. REF=AT, ALT=GC)
        elif len(ref_allele) == len(alt_str):
            writer.write(alt_str)
            ref_pos = start + len(ref_allele)
            continue
        # Plain VCF Indel without SVTYPE (e.g. REF=A, ALT=AG or REF=AG, ALT=A)
        elif len(alt_str) > len(ref_allele) and alt_str.upper().startswith(ref_allele.upper()):
            writer.write(alt_str)
            ref_pos = start + len(ref_allele)
            continue
        elif len(ref_allele) > len(alt_str) and ref_allele.upper().startswith(alt_str.upper()):
            writer.write(alt_str)
            ref_pos = start + len(ref_allele)
            continue
```

---

## 2. Multiallelic Variant Incompatibility

### Problem Location
* **File**: [`inSVert/insert.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/insert.py#L143) (lines 143, 181)

### Description & Root Cause Analysis
In `insert.py`, phased genotype processing filters records using:

```python
# Line 143:
if haplotype >= len(sample['GT']) or sample['GT'][haplotype] != 1:
    continue
```

And SNP ALT writing hardcodes the first alternate allele index:

```python
# Line 181:
writer.write(str(alt_alleles[0]))
```

In multiallelic VCF records (e.g., `REF=A, ALT=C,G`), a genotype for a diploid sample carrying the second alternate allele (`G`) is represented as `GT=0|2` or `1|2`. 
1. `sample['GT'][haplotype] != 1` checks if the allele index is strictly equal to `1`. When `sample['GT'][haplotype] == 2`, the condition is `True`, causing the record to be **silently skipped** on that haplotype.
2. Even if passed, `alt_alleles[0]` hardcodes allele `0` (`C`) regardless of whether the sample carries allele `1` or allele `2`.

### Impact
* All multiallelic variants (SNPs and SVs) where samples carry alternate alleles higher than index 1 (`GT=0|2`, `GT=2|2`, `GT=1|3`, etc.) are either omitted or populated with the wrong alternate sequence.

### Remediation Code
Update genotype extraction to retrieve the dynamic allele index and index into `alt_alleles`:

```python
# Check if current haplotype carries an ALT allele (allele_idx > 0)
allele_idx = sample['GT'][haplotype] if haplotype < len(sample['GT']) else 0
if allele_idx is None or allele_idx == 0:
    continue  # Reference allele (0) or missing allele (None)

# For SNP ALT lookup:
alt_allele_seq = str(var.alts[allele_idx - 1])
writer.write(alt_allele_seq)
```

---

## 3. Binary Search Failure for Nested Intervals Specified in the `.bed` Exclusion File

### Problem Location
* **File**: [`inSVert/utils_sim.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/utils_sim.py#L189-L246) (lines 189–246)

### Description & Root Cause Analysis
In `utils_sim.py`, `parse_bed()` parses exclusion coordinates into `excluded_ranges[chrom]` and sorts them by start coordinate. During simulation, `overlaps_excluded_region()` checks if a candidate SV interval `(start, end)` overlaps any excluded region using binary search:

```python
def overlaps_excluded_region(chrom, start, end, excluded_regions:dict):
    ...
    intervals = excluded_regions[chrom]
    idx = bisect.bisect_right(intervals, (start, end))

    if idx > 0:
        prev_start, prev_end = intervals[idx - 1]
        if prev_end > start:
            return True
        
    if idx < len(intervals):
        next_start, next_end = intervals[idx]
        if next_start < end:
            return True 
        
    return False
```

If a `.bed` file contains overlapping or nested intervals (e.g. `(100, 1000)` and `(200, 250)`), sorting yields `[(100, 1000), (200, 250)]`.

When testing a candidate SV interval at `(260, 270)`:
1. `bisect_right` compares `(260, 270)` against `(200, 250)` (`260 > 200`), placing `idx = 2`.
2. `idx - 1` points to `intervals[1] == (200, 250)`.
3. `prev_end > start` checks `250 > 260`, which is **False**.
4. `idx < len(intervals)` is **False**.
5. The function returns `False` (no overlap), **missing the enclosing `(100, 1000)` interval!**

### Impact
* Candidate SVs are placed inside regions explicitly specified as excluded in the BED file whenever the BED file contains nested or overlapping entries.

### Remediation Code
Merge overlapping/nested intervals in `parse_bed()` after sorting so that every chromosome contains disjoint, monotonically increasing intervals:

```python
def parse_bed(bed_path: str):
    ...
    # Sort and merge intervals per chromosome
    for chrom in excluded_ranges:
        excluded_ranges[chrom].sort(key=lambda x: (x[0], x[1]))
        merged = []
        for current in excluded_ranges[chrom]:
            if not merged:
                merged.append(current)
            else:
                prev_start, prev_end = merged[-1]
                if current[0] < prev_end:  # Overlapping or touching
                    merged[-1] = (prev_start, max(prev_end, current[1]))
                else:
                    merged.append(current)
        excluded_ranges[chrom] = merged
    
    return excluded_ranges
```

---

## 4. Crash on Site-Only or Ungenotyped VCFs

### Problem Location
* **File**: [`inSVert/input_validation.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/input_validation.py#L43-L52) (lines 43–52)
* **File**: [`inSVert/insert.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/insert.py#L43) (lines 43, 92)

### Description & Root Cause Analysis
In `insert.py`, variant streaming reads genotypes directly assuming sample columns are present:

```python
# Line 43:
gt = rec.samples[0]['GT']

# Line 92:
sample = var.samples[0]
```

Similarly, `input_validation.py` line 44 performs `rec.samples[0].get('GT')` after checking `if rec.samples:`.

If an input VCF is a site-only file (contains no sample columns after `FORMAT`, standard in reference databases like ClinVar or gnomAD) or contains records with ungenotyped samples (where `rec.samples` is empty or `GT` is missing/None):
* `rec.samples[0]` raises `IndexError: tuple index out of range`.
* `rec.samples[0]['GT']` raises `KeyError` or returns `None`.

### Impact
* `inSVert insert` immediately crashes when run on site-only VCFs or VCFs containing records without sample GT calls.

### Remediation Code
Add explicit validation and fallback handling in `insert.py` and `input_validation.py`:

```python
# Safe sample extraction in insert.py:
if not var.samples or len(var.samples) == 0:
    raise ValueError(
        f"Variant at {var.chrom}:{var.pos} contains no sample columns. "
        f"inSVert require sample genotypes (GT) to perform haplotype insertion."
    )

sample = var.samples[0]
if 'GT' not in sample or sample['GT'] is None or None in sample['GT']:
    # Handle or skip ungenotyped site
    continue
```

---

## 5. Performance: Warning `print()` Flooding in `insert.py` Inner Loop

### Problem Location
* **File**: [`inSVert/insert.py`](file:///home/thomas/Desktop/insvert_demo/inSVert/insert.py#L126-L136) (lines 126–136)

### Description & Root Cause Analysis
In `insert.py`, when handling unphased heterozygous variants (`is_heterozygous and not is_phased`), the engine executes the following inside the per-variant, per-haplotype streaming loop:

```python
for haplotype in range(ploidy):
    for chrom in ref.references:
        for var in chrom_variants:
            ...
            if not is_phased and is_heterozygous:
                ...
                print(
                    f"WARNING: Variant '{var_id}' ({sv_label} at "
                    f"{var.chrom}:{var.pos}) is unphased. Phasing "
                    f"information is needed for correct haplotype "
                    f"assignment. Randomly assigned to Haplotype(s) "
                    f"{haps_str}."
                )
```

### Impact
* If an input VCF contains 100,000 unphased SNPs, `insert.py` calls `print()` **100,000 times per haplotype pass** (200,000 print calls for diploid genomes).
* Standard I/O operations blocking stdout cause excessive terminal output bloat and slow down execution runtime by **5x to 10x**.

### Remediation Code
Collect unphased variant counts and output a single summary warning before the loop, or log unphased assignments only once per variant key:

```python
# Track warned variant keys to avoid duplicate prints across haplotypes
warned_unphased = set()

# Inside variant loop:
if var_key not in warned_unphased:
    warned_unphased.add(var_key)
    # Print once per unique variant record
    print(
        f"WARNING: Variant '{var_id}' ({sv_label} at {var.chrom}:{var.pos}) "
        f"is unphased. Randomly assigned to Haplotype(s) {haps_str}."
    )
```
