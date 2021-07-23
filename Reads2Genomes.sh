#### De novo assembly, polishing, scaffolding and evaluation
# need tools
# NanoPlot
# nextdenovo
# minimap2
# racon
# medaka
# fastqc
# trimmomatic
# bowtie2
# samtools
# plion
# ragoo(ragtag)
# nucmer
# quast.py

ROOT=rice_project

# mkdir raw_ont
# mkdir raw_sgs
# mkdir ref

# mkdir qc_ont
# mkdir qc_sgs
# mkdir trim_ont
# mkdir trim_sgs
# mkdir assemble
# mkdir racon
# mkdir medaka
# mkdir plion
# mkdir ragoo

#######
# QC long reads (NanoPlot)
#######
cd $ROOT
NanoPlot -t 4 \
	--fastq raw_ont/${sample}.fq.gz \
	-o qc_ont/${sample} \
	-p ${sample} \
	--plots hex dot


#######
# Trim long reads with quality >= 7 , length >=1000 (NanoFilt)
#######
cd $ROOT
zcat raw_ont/${sample}.fq.gz | NanoFilt -q 7 -l 1000 | gzip > trim_ont/${sample}.fq.gz



#######
# Assemble (Nextdenovo)
#######
cd $ROOT
mdkir assemble/${sample}

## prepare input files
cd $ROOT
ls trim_ont/${sample}.fq.gz > assemble/${sample}/${sample}.input.fofn
cd assemble/${sample}

## find suitable seed_cutoff length
# seq_stat -f 1000 -g 390Mb ${sample}.input.fofn

## prepare nextdenovo parameters 
cp ../nextdenovo_template.cfg ./${sample}.nextdenovo.cfg

# nextdenovo_template.cfg 
:<<EOF
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes
deltmp = yes
rerun = 3
parallel_jobs = 10
input_type = raw
input_fofn = input.fofn
workdir = rundir
[correct_option]
read_cutoff = 1000
seed_cutoff = 20000 # 15000 
blocksize = 16g
pa_correction = 10
seed_cutfiles = 10
sort_options = -m 16g -t 4 -k 40
minimap2_options_raw = -x ava-ont -t 4
correction_options = -p 10

[assemble_option]
random_round = 20
minimap2_options_cns = -x ava-ont -t 40 -k17 -w17
nextgraph_options = -a 1
EOF

## main assemble program
nextDenovo ${sample}.nextdenovo.cfg
ass=assemble/${sample}/rundir/03.ctg_graph/01.ctg_graph.sh.work/ctg_graph00/nextgraph.assembly.contig.fasta



########
# Polish with long reads (Racon & Medaka)
########

# Racon
cd $ROOT
mdkir racon/${sample}

## map
minimap2 -d racon/${sample}/${sample}_frag.mmi ${ass}
minimap2 -x map-ont \
	-t ${th} \
	--MD ${sample}/${sample}_frag.mmi \
	-o ${sample}/${sample}_readtoass.paf \
	trim_ont/${sample}.fq.gz

## recommend parameter for medaka
racon -t ${th} \
	-m 8 \
	-x -6 \
	-g -8 \
	-w 500 \
	trim_ont/${sample}.fq.gz ${sample}/${sample}_readtoass.paf ${ass} > ${sample}/${sample}_racon.fa
ass_racon=racon/${sample}/${sample}_racon.fa


# Medaka
cd $ROOT
mkdir medaka/${sample}
medaka_consensus -i trim_ont/${sample}.fq.gz \
	-d $ass_racon \
	-o $sample \
	-t $th \
	-m r941_min_high_g303
ass_medaka=medaka/${sample}/${sample}_medaka.fa


#######
# QC short reads (FastQC)
#######
cd $ROOT
mdkir qc_sgs/${sample}
file1=raw_sgs/${sample}_R1.fq.gz
file2=raw_sgs/${sample}_R2.fq.gz
fastqc -t 4 \
	-o qc_sgs/${sample} \
	${file1}

fastqc -t 4 \
	-o qc_sgs/${sample} \
	${file1}

########
# Trim short reads, clip adapter, cut head 10bp (Trimmomatic)
########
cd $ROOT
mdkir trim_sgs/${sample}
file1=raw_sgs/${sample}_R1.fq.gz
file2=raw_sgs/${sample}_R2.fq.gz
upd=unpair
adir=~/tool/Trimmomatic-0.39/adapters
java -jar ~/tool/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${file1} ${file2} \
	trim_sgs/${sample}/${sample}_trim_R1.fq.gz trim_sgs/${sample}/${upd}/${sample}_unpaired_R1.fq.gz \
	trim_sgs/${sample}/${sample}_trim_R2.fq.gz trim_sgs/${sample}/${upd}/${sample}_unpaired_R2.fq.gz \
	ILLUMINACLIP:${adir}/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:10



########
# Polish with short reads (Pilon)
########
cd $ROOT
mkdir pilon/${sample}
th=8
ngs1=trim_sgs/${sample}_R1.fq.gz
ngs2=trim_sgs/${sample}_R2.fq.gz

## map, bowtie2
bowtie2-build --threads ${th} ${ass_medaka} pilon/${sample}/idx
bowtie2 -p ${th} \
	-x pilon/${sample}/idx \
	-S pilon/${sample}/${sample}_readtoassmedaka_ngs.sam \
	-1 ${ngs1} \
	-2 ${ngs2} \
	2>pilon/${sample}/${sample}.mapping
# sort 
samtools sort -@ ${th} \
	-o pilon/${sample}/${sample}_ngs.bam \
	pilon/${sample}/${sample}_ngs.sam
samtools index pilon/${sample}/${sample}_ngs.bam

java -Xmx180G \
	-jar ~/tool/pilon-1.23.jar \
	--threads ${th} \
	--genome ${ass_medaka} \
	--frags pilon/${sample}/${sample}_readtoassmedaka_ngs.bam \
	--outdir pilon/${sample}
ass_plion=pilon/${sample}/${sample}_pilon.fa



########
# Scaffold, Reference-guided (RaGOO RagTag)
########
cd $ROOT
mkdir ragoo/${sample}

## get corrected reads
cat assemble/${sample}/rundir/02.cns_align/01.get_cns.sh.work/get_cns*/cns.fasta > ragoo/${sample}/${sample}.corr.fa
cor_ont_reads=ragoo/${sample}/${sample}.corr.fa
ref=ref/msu7/chrs.fa
query=$ass_plion
th=8
nucmer=~/tool/mummer-3.9.4alpha/nucmer

## break misassemblies in assembly with corrected reads
ragtag.py correct \
	-t ${th} \ 
	-u \
	--aligner ${nucmer} \
	-o ragoo/${sample}/${sample}_correct \
	-T corr \
	-R ${cor_ont_reads} \
	${ref} ${query}

## scaffold the contigs 
ragtag.py scaffold \
	-t ${th} \
	-u \
	--aligner ${nucmer} \
	-o ragoo/${sample}/${sample}_scaffold \
	${ref} ragoo/${sample}/${sample}_correct/${sample}_pilon.corrected.fasta
chromosome_level_scaffold=ragoo/${sample}/${sample}_scaffold/ragtag.scaffolds.fasta



########
# Quast
########
# cd $ROOT
# mkdir quast/${sample}
ref=ref/msu7/chrs.fa
th=4
reggff=ef/msu7/chrs.gff
quast.py -t ${th} \
	--large ${chromosome_level_scaffold} \
    -r ${ref} \
    -g ${refgff} \
    -o quast/${sample} \
    --no-icarus \
    --no-snps \
    --min-identity 90.0

