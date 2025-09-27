# demo specififcations 


### command-line interface design
```
inSVert <command> [options]

Commands:
    simulate    -analyze input VCF and produce a realistic fake
    insert      -inserts the SVs into the reference genome
    pipeline    -combines the above in a single process
```

**simulate** 

input : input.vcf

output : simulated.vcf 

1. parse the input VCF file seprate Svs by SV type (DEL,INS,INV,DUP,DUP:TANDEM)
2. for each SV type, get the lengths and fit that data to a power-law distribution
3. based on user input or based on the instances of each SV type in the original VCF, sample a length randomly for a SV
4. build a SV event using that length into a structured VCF file
5. repeat steps 3 and 4 until the target amount of SVs are produced

**insert**

input : simulated.vcf + input.fa 

output : simulated.fa 


------
------

### practical considerations

due to limited computing power, the organism of choice will be yeast.

1. create a new vcf file: 
for this purpose I used Saccharomyces [kudriavzevii](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR7517606&display=download) ONT reads mapped to [Saccharomyces cervisiae](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/) reference genome, in this fashion i am sure to find a rich profile of structural variants when producing a VCF.

    mapping is done using [minimap2](https://github.com/lh3/minimap2), which is particularly suited for mapping long reads. A preliminary variant calling on the resulting bam file is performed with [sniffles](https://github.com/fritzsedlazeck/Sniffles) 

