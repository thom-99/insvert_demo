# inSVert

inSVert is a command-line tool for simulating genomic variants and applying them to a reference FASTA.

It can generate a VCF containing simulated variants or apply variants from an existing VCF to produce modified haplotype sequences.

inSVert is open source and freely available on [GitHub](https://github.com/thom-99/inSVert).

## What can inSVert do?

inSVert provides two main commands:

### Simulate variants

The `simulate` command generates variants across a reference genome according to a YAML configuration file.

The simulation can control:

- the number and length distribution of structural variants
- SNP and MNP densities
- ploidy and heterozygosity
- duplication copy numbers
- translocation orientation
- regions excluded through a BED file
- reproducibility through a random seed

### Apply variants to a reference

The `insert` command applies variants from a VCF to a reference FASTA.

Variants are applied separately to each haplotype according to their genotypes. The resulting haplotypes can be written to one combined FASTA or to separate files.

## Supported variants

inSVert supports:

- single-nucleotide polymorphisms
- multi-nucleotide polymorphisms
- insertions
- deletions
- inversions
- tandem duplications
- interspersed duplications
- cut-and-paste translocations

Some variants, particularly interspersed duplications and translocations, require a specific breakend representation in the input VCF.


## Under the hood

For readers interested in how inSVert works internally:

- [Variant simulation](under-the-hood/variant-simulation.md) explains variant overlap checking, haplotype-aware coordinate tracking, BED exclusions, and variant length distributions.
- [Variant insertion](under-the-hood/variant-insertion.md) explains reference-coordinate tracking as variant are applied and memory-efficient genome construction.

