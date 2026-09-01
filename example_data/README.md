# inSVert example datasets

Each directory contains a small reference FASTA and a `config.yaml` that can
be used to try the `simulate` and `insert` commands without downloading a full
genome. The references are representative subsets, not complete assemblies.

Always insert a VCF into the exact reference FASTA used to create it:

```bash
inSVert simulate config.yaml reference.fa --seed 42 -o variants.vcf
inSVert insert reference.fa variants.vcf --ploidy <PLOIDY> -o edited.fa
```

## `human/`

Contains a small multi-contig human reference subset and a diploid (`ploidy: 2`)
configuration. It is designed to demonstrate many SNPs and MNPs together with
translocations, while keeping deletions, insertions, duplications, and
inversions comparatively rare. The `TRA_COPY` and `TRA_CUT` settings use
30 kb events and include a 20% probability of reverse orientation.

if you have space, try downloading the full human genome
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip hg38.fa.gz
```
and going crazy with the configfile!

## `yeast/`

Contains a compact *Saccharomyces cerevisiae* reference and a haploid
(`ploidy: 1`) configuration. It is a smaller structural-variant example with
modest SNP/MNP density, short deletions and insertions, duplications,
inversions, and a few translocations. Its shorter event lengths are chosen to
fit within the relatively small yeast chromosomes.
