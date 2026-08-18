<p align="center">
  <img src="img/small_logo.png" alt="inSVert Logo" width="300">
</p>

# inSVert 
inSVert is a bioinformatic utility deisgned to perform two distinct, yet interconnected, tasks:
- simulate genomic variants, encoding them into a VCF file
- insert variants from a VCF file into a reference genome

### so, what can I do with it?
On principle, inSVert main purpose lies on producing a synthetic ground-truth set of structural variants and inserting them into a reference genome, which will be further processed with the ultimate goal of benchmarking other tools such as variant callers and read mappers. 

# Installation
when installing inSVert, using a virtual enviroment is highly recommended

ensure your packaging tools are updated:
```bash
python3 -m pip install --upgrade pip setuptools wheel
```
install inSVert with:
```bash
git clone https://github.com/thom-99/inSVert
cd inSVert
pip install .
```
# Usage

## inSVert simulate
The first module simulates a custom set of variants according to the user instructions provided in the config.yaml file.
You can use the configfile.yaml available in this repo or generate a template one with 

```
inSVert generate-configfile
```
then edit it based on what you want to simulate

The user can choose to simulate variants according to a [pareto distribution](https://en.wikipedia.org/wiki/Pareto_distribution), which more closely reflects the natural distribution of variants (with fewer long variants and more short variants), or a [normal distribution](https://en.wikipedia.org/wiki/Normal_distribution).  

inSVert also takes into account polyploid organisms: the user uses the 'ploidy' and 'heterozygousity' parameters to instruct the simulate module about how many genome copies he intends to simulate the variants on (most likely this corresponds to the ploidy number of the organism of interest) and the probability of variants of being heterozygous (present only on one genome copy). 

to simulate variants, simply type 
```
inSVert simulate config.yaml reference.fasta -o simulated.vcf
```
where the first argument is the path to the config.yaml file and the second one the path to your reference genome in fasta format.

optional arguments:
```
-o, --output PATH
    Path for the output VCF file. Default: simulated.vcf.

--seed INTEGER
    Random seed for reproducible simulations.

--exclude PATH
    BED file containing genomic regions to exclude from simulation, such as
    centromeres or mitochondrial DNA.

--non-symbolic
     Write explicit REF and ALT sequences instead of symbolic alleles such as
    <INS> and <DEL>. This can substantially increase the output VCF size
```

## inSVert insert
given a VCF file , either produced by *inSVert simulate* or provided by the user, the  variants contained in the file will be programmatically inserted into a specified reference genome in fasta format. Although it may seem trivial, this is by far the most complex step as it requires careful tracking of the inserted variants to avoid indexing problems and to avoid placing variants one on top of the other. 

For this reason it is a strict requirement that the VCF file is produced from the same reference in which we are trying to insert the variants and that the VCF file is sorted. Therefore, inSVert check both requiremenrs and take care of sorting the VCF file if not already sorted.


to insert variants from a sorted VCF to a reference genome, simply type 
```
inSVert insert reference.fasta simulated.vcf --ploidy 2 -o simulated.fasta
```
where the first argument is the path to the reference genome and the second one the path to the VCF chosen by the user; the --ploidy argument is not optional and requires to specify how many copies of the genome to simulate. If you are using inSVert simulate to produce a VCF, it has to match the ploidy argument of the config.yaml. In any case, the genotype string of your variants in the VCF should be informative about the ploidy number you need to insert here.  

optional arguments:

```
-o, --output PATH
    Path for the output FASTA file.

--truth-vcf PATH
    Write a ground-truth VCF containing only variants successfully inserted
    into the output FASTA. Recommended for third-party VCFs, where malformed,
    unsupported, or conflicting records may be skipped.

--gc FLOAT
    GC fraction used when generating insertion sequences. Must be between
    0 and 1. Default: 0.41.

--sample-name TEXT
    Sample name used in multi-haplotype FASTA headers. Default: Sample.

--skip-unphased-variants
    Skip unphased variants instead of assigning them randomly to a haplotype.

--split-haplotypes
    Write each haplotype to a separate FASTA file. For example, an output
    named output.fa with ploidy 2 produces output_hap1.fa and output_hap2
```


# Architecture 

![Alt text](img/inSVert_whitebackground.png)

inSVert has a decoupled architecture, designed so that its modules can be used as a standalone bioinformatic utility. 

