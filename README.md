TGSRICEPAN
-

# Introduction
Rice is one of the most important crops for human. The third-generation sequencing (TGS, also called long-read sequencing, LRS) helps us assemble more high-quality genomes and construct more complete pan-genomes. Here are codes of the article "**Long-read sequencing of 111 rice genomes reveals significantly larger pan-genomes**".

# Codes
The main pipelines are as following and the main self-writen scripts are in "scripts" directory.

### Reads2Genomes.sh 
The pipeline of preprocessing, assembling, polishing and scaffolding from raw long/short reads to assemblies.

### Reads2SV.sh
The pipeline of SV calling, filtering and merging from long reads to SVs.

### Gapfill.sh
The pipeline of filling gaps in Nipponbare genomes with corrected reads and assembled contigs.

### Genomes2Pan.sh
The pipeline of constructing pan-genomes from genomes.
