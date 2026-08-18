# Variant simulation



## Overlap checking

When inSVert proposes a variant, it must first check whether that variant would conflict with one that has already been accepted. Repeatedly comparing the new variant with every previous variant would become slow as the VCF grows. Instead, inSVert keeps accepted intervals in sorted lists and uses binary search to find the only nearby intervals that could overlap.

### Variants are tracked separately for each haplotype

The interval tracker has the following structure:

```python
# {chromosome: {haplotype_index: [(start, end), ...]}}
sv_positions = defaultdict(lambda: defaultdict(list))
```

For example, in a diploid simulation, chromosome 1 has one sorted interval list for haplotype 1 and another for haplotype 2. A heterozygous variant with genotype `0|1` is stored only in the second list.

This distinction is important. Two variants can occupy the same reference coordinates without conflicting if they occur on different haplotypes. inSVert therefore checks only the haplotypes whose genotype contains allele `1`

### Coordinates stored by the tracker

The tracker uses 1-based reference coordinates, matching VCF `POS`. Its intervals are treated as closed at both ends. This means that an interval `(100, 200)` reserves positions 100 through 200.

| Variant | Interval reserved by the simulator | Why |
| --- | --- | --- |
| `INS` | `(POS, POS)` | The insertion occurs at one anchored breakpoint. Its inserted length does not consume reference bases. |
| `DEL` | `(POS, POS + abs(SVLEN))` | The VCF anchor and the deleted reference span are reserved. |
| `INV` | `(POS, POS + SVLEN)` | The VCF anchor and the inverted reference span are reserved. |
| Tandem `DUP` | `(POS, POS + SVLEN)` | The VCF anchor and duplicated reference unit are reserved. |
| `TRA_COPY` destination | `(POS, POS + 1)` | The two sides of the insertion junction are reserved. The source remains unchanged. |
| `TRA_CUT` destination | `(POS, POS + 1)` | The two sides of the insertion junction are reserved. |
| `TRA_CUT` source | `(POS, POS + SVLEN + 1)` | The source segment and the two source junctions are reserved. |
| SNP or MNP | The replaced reference bases | A SNP consumes one base. An MNP consumes the length of its REF allele. |



### Finding a possible overlap with `bisect`

Python's [`bisect`](https://docs.python.org/3/library/bisect.html) module searches an already sorted list without scanning it from the beginning. `bisect_right()` returns the index at which the proposed interval could be inserted while preserving the sort order.

Because accepted intervals in each list do not overlap, only two existing intervals can conflict with the proposal:

1. the interval immediately before the insertion point
2. the interval immediately after the insertion point

The core check is therefore small:

```python
idx = bisect.bisect_right(intervals, (start, end))

if idx > 0:
    _, prev_end = intervals[idx - 1]
    if prev_end >= start:
        return True

if idx < len(intervals):
    next_start, _ = intervals[idx]
    if next_start <= end:
        return True
```

If neither neighbor overlaps, the candidate is safe for the haplotypes being checked. Once accepted, `bisect.insort()` places it at the correct position in the same sorted list:

```python
bisect.insort(sv_positions[chrom][hap_idx], (start, end))
```



### Excluding regions from a BED file

BED exclusion intervals are stored by chromosome, sorted, and merged when they overlap or touch. inSVert then uses the same binary-search idea to compare each candidate with the two neighboring BED intervals.

BED coordinates are 0-based and half-open. For example, `[99, 100)` identifies the single base reported as position 100 in a VCF. Any VCF-derived coordinate must be converted to this convention before it is compared with a BED interval. SNPs and MNPs currently perform this conversion explicitly:

```python
bed_start = pos - 1
bed_end = bed_start + allele_length
```

For `DEL`, `INV`, and tandem `DUP`, the simulator passes `(POS, END)` to the half-open BED check. This excludes the VCF anchor at `POS` and covers the reference span after it, which is the span consumed or rearranged by the operation. 

## Building and sampling length distributions

For most structural variants, the configuration provides a count, a distribution, and length limits. inSVert samples the requested number of lengths before it starts choosing genomic positions.

```yaml
DEL:
  count: 100
  distribution: normal
  parameters:
    median_length: 5000
    sigma: 500
    min_length: 50
    max_length: 100000
```

### Normal distribution

For a normal distribution, `median_length` is used as the centre of the distribution. Since a normal distribution is symmetric, its mean and median are the same. `sigma` controls how widely lengths vary around that centre. If `sigma` is omitted, inSVert uses 10 percent of `median_length`.

```python
sample = random.gauss(mu, sigma)
```

### Pareto distribution

A Pareto distribution produces many shorter variants and a smaller number of much longer variants, resembling real-world variant lenghts. The configuration supplies the desired median and minimum. inSVert converts them into the Pareto shape parameter `alpha`:

$$
\alpha = \frac{\log(2)}{\log(\text{median length} / \text{minimum length})}
$$

It then samples a standard Pareto value and scales it by the configured minimum:

```python
sample = random.paretovariate(alpha) * min_length
```


