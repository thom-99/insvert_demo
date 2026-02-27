# demo specififcations 

inSVert is a toolkit for the simulation of structural variants and for the insertion of structural variants into a reference genome. 

inSVert main utility lies in benchmarking different read mappers and variant callers against a ground thruth set of structural variants. The software is composed by two modules: simulate & insert. 

![Alt text](img/benchmarking_workflow.png)

### inSVert simulate
The first module simulates a custom set of structural variants such as Deletions, Insertions, Inversion and Duplications according to the user instructions provided in the config.yaml file. The default options models the SV length distribution as a lognormal distribution providing a pattern that more closely resemles real variants (fewer very long variants and many short ones), however there are other distributions to choose from. 
(to implement) inSVert also aims to be the first structural variantion simulator that is fleible in terms of ploidy, so that It can work also with genomes having a ploidy number of 3 or above, like many plant genomes do. For this reason the user needs to specify the ploidy number and a measure of heterozygousity.    

to simulate structural variants, simply type 
```
inSVert simulate config.yaml reference.fasta -o simulated.vcf
```
where the first argument is the path to the config.yaml file and the second one the path to your reference genome in fasta format, you can specify in which file you want your simulated SVs after the -o option. 



### inSVert insert
given a VCF file , either produced by *inSVert simulate* or provided by the user, the Structural Variants contained in the file will be programmatically inserted into a specified reference genome in fasta format. Although it may seem trivial, this is by far the most complex step as it requires careful tracking of the inserted variants to avoid indexing problems and to avoid placing variants one on top of the other. 

It is a strict requirement that the VCF file is produced from the same reference in which we are trying to insert the variants. 
This can be easily checked by inspecting the first few lines of the VCF
simply type 
```
head myfile.vcf 
```
and check for a correspondance between the reference of the VCF and the one you want to put the variants in. 

to insert Structural Variants from a sorted VCF to a reference genome, simply type 
```
inSVert insert reference.fasta simulated.vcf -gc 0.41 -o simulated.fasta
```
where the first argument is the path to the reference genome and the second one the path to the VCF chosen by the user. 
The optional argument -gc allows the user to specify the GC ratio of their reference genome. This is used when building insertions, in order to make DNA sequences more realistic. The default is set to the human genome GC content (0.41). 


---





### practical considerations

- due to limited computing power, the organism of choice will be yeast.
- checked validity of the VCF file with vcftools (vcf-validator)
- checked the effects of inserting SVs by pairing and plotting the un-edited and edited fasta references with Gepard. 

---

# TO DO


for the final version:


- allow to simulate based on other distributions (student and normal) 
- add Translocations, movements of DNA from one chromosome to another one. (if it does not mess with polyploids) 
- handle polidy number (hardest task yet)
specifically:
- In the run function, generate this GT string for every SV and pass it to the VariantObjects constructor.
- Update StructuralVariant.__init__ to accept and store a genotype string
- Modify the format() method to replace the hardcoded 1/1 with the instance's self.genotype attribute
- loop restructuring in insert_streaming.py : wrap the entire chromsome processing loop inside a new loop (for h_idx in range(ploidy))
- Update the FASTA header writer to use the Sample#H{index}#Chrom format (e.g., Sample#H1#chr1) to ensure the output is pangenome-compatible and easy to distinguish. Add an optional setting to be able to replace the default 'Sample' with the organism name. 
- Inside the variant loop, use pysam to access the GT field for the current record
- Allele Check: Only execute the apply_insertion/deletion/etc. logic if the genotype bit at the current h_idx is 1. If it is 0, skip the variant and treat that region as reference for that specific haplotype pass.
- Update the run function to either accept a ploidy argument or automatically determine it by inspecting the first record of the VCF file. Implement a test that if it finds a genotype with more copies than plody, it gives an error. 
CLI Updates (cli.py):Update the simulate command to pass the new config parameters to the simulation engine.Update the insert command to handle the increased output file size (which will be original reference size $\times$ ploidy).


- containerize in docker image 
- write a nextflow benchmarking pipeline 



