# demo specififcations 

inSVert is a toolkit for the simulation of structural variants and for the insertion of structural variants into a reference genome. 

inSVert main utility lies in benchmarking different read mappers and variant callers against a ground thruth set of structural variants. The software is composed by two modules: simulate & insert. 

### inSVert simulate
The first module simulates a custom set of structural variants such as Deletions, Insertions, Inversion and Duplications according to the user instructions provided in the config.yaml file. The default options models the SV length distribution as a lognormal distribution providing a pattern that more closely resemles real variants (fewer very long variants and many short ones), however there are other distributions to choose from. 
(to implement) inSVert also aims to be the first structural variantion simulator that is fleible in terms of ploidy, so that It can work also with genomes having a ploidy number of 3 or above, like many plant genomes do. For this reason the user needs to specify the ploidy number and a measure of heterozygousity.    

to simulate structural variants, simply type 
```
inSVert simulate config.yaml reference.fasta -o simulated.vcf
```
where the first argument is the path to the config.yaml file and the second one the path to your reference genome in fasta format, you can specify in which file you want your simulated SVs after the -o option. 



### inSVert insert
given an sorted VCF file, either produced by *inSVert simulate* or provided by the user, the Structural Variants contained in the file will be programmatically inserted into a specified reference genome in fasta format. Although it may seem trivial, this is by far the most complex step as it requires careful tracking of the inserted variants to avoid indexing problems and to avoid placing variants one on top of the other. 

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

bottlenecks:
- the whole reference genome has to be loaded into memory, sucking up a huge amount of RAM


### practical considerations

- due to limited computing power, the organism of choice will be yeast.
- checked validity of the VCF file with vcftools (vcf-validator)
- checked the effects of inserting SVs by pairing and plotting the un-edited and edited fasta references with Gepard. 

---

# TO DO


for the final version:

- allow to simulate based on other distributions (student and normal) 
- use interval tree to compute the overlap wich is a more efficient solution O(logN) rather than the current O(N^2) 
- add a forced sorting step for the VCF, otherwise if the VCF to insert is not sorted the program will break.
- add Translocations, movements of DNA from one chromosome to another one. (if it does not mess with polyploids) 
- handle polidy number (hardest task yet)
- reduce memory usage by implementing some form of lazy-loading (reading chromsomes??)



