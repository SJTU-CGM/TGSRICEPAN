#### SV calling, merging and filtering
# need tools
# minimap2
# sniffle
# SURVIVOR

ROOT=rice_project
#######
# SV calling (minimap2, sniffles)
#######
# map
th=16
refidx=ref/msu7/chrs
mkdir mapmsu7
minimap2 -t ${th} \
	--MD \
	-a ${refidx}.mmi \
	-o mapmsu7/${sample}_ont.sam \
	trim_ont/${sample}.fq.gz

# sort
samtools sort -@ ${th} \
	-o mapmsu7/${sample}_ont.bam \
	mapmsu7/${sample}_ont.sam

mkdir sv_sniffles
# call SV
mkdir sv_sniffles/raw
sniffles -m mapmsu7/${sample}_ont.bam \
	--tmp_file tmp \
	-v sv_sniffles/raw/${sample}.vcf

#######
# SV filter (SURVIVOR)
#######
mkdir sv_sniffles/filter 
ls sv_sniffles/raw/*vcf > sv_sniffles/raw_vcf_files.txt
while read line;
do
	outfile=${line##*/}
	SURVIVOR filter $line NA 50 -1 0.05 10 filter/${outfile}
done < raw_vcf_files.txt

#######
# SV merge (SURVIVOR)
#######
ls filter/*vcf > filter_vcf_files.txt
SURVIVOR merge filter_vcf_files.txt 1000 0 1 1 0 50 merged_SURVIVOR_allsample.vcf


