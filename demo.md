# demo specififcations 

due to limited computing power, the organism of choice will be yeast.

1. create a new vcf file 
for this purpose I used Saccharomyces [kudriavzevii](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR7517606&display=download) ONT reads mapped to [Saccharomyces cervisiae](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/) reference genome, in this fashion i am sure to find a rich profile of structural variants when producing a VCF.

mapping is done using [minimap2](https://github.com/lh3/minimap2), which is particularly suited for mapping long reads. A preliminary variant calling on the resulting bam file is performed with [sniffles](https://github.com/fritzsedlazeck/Sniffles) 

