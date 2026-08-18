# Variant insertion

inSVert applies a VCF to a reference FASTA while keeping all VCF positions tied to the original reference. Two ideas make this possible:

1. variants are processed in reference-coordinate order
2. a cursor records how far inSVert has read from the original reference

Together, these ideas avoid recalculating every later coordinate after each insertion or deletion.

## Keeping track of inserted variants

### The coordinate problem

Suppose an insertion adds 1,000 bases near the beginning of a chromosome. Every base after it now appears 1,000 positions later in the output sequence and we need to account for this in some way. If inSVert rewrote all later VCF coordinates after every edit, coordinate handling would quickly become complicated.

Instead, VCF coordinates always continue to refer to the unchanged input reference genome. inSVert tracks its advancements through that reference with `ref_pos`.

`ref_pos` is a 0-based cursor pointing to the first reference base that has not yet been written, replaced, or skipped. Essentially it points to just before a variant is applied. 

### The algorithm

Before insertion begins, inSVert makes sure the VCF is coordinate-sorted, BGZF-compressed, and indexed. It can then visit each chromosome's variants from 'left to right'.

For every variant, the 1-based VCF position is converted to a 0-based Python coordinate:

```python
start = var.pos - 1
```

If there is unchanged reference sequence between `ref_pos` and the next variant, that sequence is copied first:

```python
if start > ref_pos:
    chunk = ref.fetch(chrom, ref_pos, start)
    writer.write(chunk)
    ref_pos = start
```

inSVert then applies the variant and updates `ref_pos` according to how many reference bases the operation consumed.

####  each variant type edits the genome in its own way, therefore the ref_pos has to be modified accordingly
| Operation | Sequence written | Effect on `ref_pos` |
| --- | --- | --- |
| SNP or MNP | The explicit ALT allele | Advance by the length of REF |
| Insertion | The anchor base, then inserted sequence | Advance past the anchor only |
| Deletion | The anchor base | Skip the deleted reference span |
| Inversion | The anchor base, then the reverse-complemented span | Advance past the inverted span |
| Tandem duplication | The anchor base, then the reference unit repeated according to `CN` | Advance past the original unit once |
| Copy translocation at destination | Sequence fetched from the source | Do not consume source bases at the destination |
| Cut translocation at source | No sequence for the cut span | Skip the removed source span |

The key rule is simple: `ref_pos` moves according to the input reference, not according to how much output was produced.

An insertion can write thousands of new bases while leaving the cursor immediately after its reference anchor. A deletion can move the cursor forward while writing none of the deleted bases. This is why changes in output length do not shift the coordinates of variants that appear later in the VCF.

### A small example

Consider a reference with a 3-base deletion followed by an insertion:

```text
reference:     ABCDEFGHIJ
delete:           DEF
insert after I:   XYZ
```

inSVert processes it as follows:

1. Copy `ABC`, the unchanged sequence before the deletion.
2. Move `ref_pos` past `DEF` without writing those bases.
3. Copy `GHI`, the unchanged sequence before the insertion.
4. Write `XYZ` without moving further through the reference.
5. Copy the remaining `J`.

The output is `ABCGHIXYZJ`. Both variants were interpreted using positions on the original `ABCDEFGHIJ` reference.



## Memory-efficient variant insertion

### The central idea

A simple way to modify a genome would be to load an entire chromosome into a string, apply every change to that string, and write it when finished. For large genomes and multiple haplotypes, those chromosome-sized strings can consume substantial memory.

inSVert instead behaves like someone copying a book while following a sorted list of corrections:

1. copy text up to the next correction
2. apply that correction
3. continue from the next unread position in the original book

At no point does inSVert need to build the complete modified genome in memory.

### Reading the reference, in chunks

The reference FASTA is indexed with `pysam.FastaFile`. This allows inSVert to fetch a specific interval without reading every earlier base:

```python
chunk = ref.fetch(chrom, ref_pos, start)
writer.write(chunk)
```

The interval before a variant is written immediately. After the final variant on a chromosome, the remaining reference interval is fetched and written in the same way.

Sequence needed for an inversion, duplication, or translocation is also fetched only when that operation is reached. For example, a translocated segment is loaded from its source coordinates when inSVert processes the destination:

```python
seq = ref.fetch(src_chr, s_start, s_end)
if is_inverted:
    seq = reverse_complement(seq)

writer.write(seq)
```

### Writing the output, incrementally

Output sequence is sent to `BufferWriter`, a class used for convenience to write DNA to the output fasta while wrapping it in fixed-width 60 char lines, as soon as it is available. 

Conceptually, the data flow is:

```text
indexed FASTA interval -> apply one edit -> FASTA line writer -> output file
```

The output file grows continuously on disk. It is not assembled as one genome-sized Python string.

