# needs tools
# quast.py
# minimap2
# faSomeRecords
# gclust
# blastn
# eupan

# needs scripts
# unalign_drop_direct.py
# paf_select2drop_cov.py
# getTax.py
# filter_exclude_contig.py
# cluster2fasta.py
# elongate_seq.py (optional)

ROOT=rice_project

########
# Quast
########
# cd $ROOT
# mkdir quast/${sample}
ref=ref/msu7/chrs.fa
th=4
reggff=ef/msu7/chrs.gff
#quast.py -t ${th} \
#	--large ${chromosome_level_scaffold} \
#    -r ${ref} \
#    -g ${refgff} \
#    -o quast/${sample} \
#    --no-icarus \
#    --no-snps \
#    --min-identity 90.0

########
# Get unalign regions from quast
########
ref=ref/msu7/chrs.fa
mkdir unalign
python3 unalign_drop_direct.py quast/${sample}/contigs_reports/contigs_report_ragtag-scaffolds.unaligned.info \
	${chromosome_level_scaffold} \
	unalign/${sample}_unalign_500.fa \
	500 \
	${sample}

# remap to msu7 chromosomes & mitochondrion/plastid
ref=ref/msu7/chrs.fa
minimap2 -x asm10 ref/chrs.fa unalign/${sample}_unalign_500.fa > unalign/${sample}_unalign_tomsu7.paf
python3 ~/100rice/code/paf_select2drop_cov.py drop unalign/${sample}_tomsu7.paf 0.8 > unalign/${sample}_unalign_inmsu7.list
~/tool/faSomeRecords unalign/${sample}_unalign_500.fa -exclude unalign/${sample}_unalign_inmsu7.list unalign/${sample}_nomsu7.fa

ref=ref/msu7/chrmipl.fa
minimap2 -x asm10 ref/chrs.fa unalign/${sample}_nomsu7.fa > unalign/${sample}_unalign_tomipl.paf
python3 ~/100rice/code/paf_select2drop_cov.py drop unalign/${sample}_tomipl.paf 0.8 > unalign/${sample}_unalign_inmipl.list
~/tool/faSomeRecords unalign/${sample}_unalign_nomsu7.fa -exclude unalign/${sample}_unalign_inmipl.list unalign/${sample}_nomsu7_nomipl.fa

#######
# Merge
#######
for file in unalign/*_nomsu7_nomipl.fa;
do
	cat $file >> unalign/unalign_no2_merge.fa
done


#######
# Cluster sequences (Gclust & Blast)
#######
# gclust
mkdir gclust
~/tool/gclust/script/sortgenome.pl --genomes-file unalign/unalign_no2_merge.fa \
	--sortedgenomes-file gclust/unalign_no2_merge.sorted.fa
~/tool/gclust/gclust -loadall \
	-minlen 20 \
	-both \
	-nuc \
	-threads 40 \
	-ext 1 \
	-sparse 2 \
	-memiden 90 \
	gclust/unalign_no2_merge.sorted.fa > gclust/unalign.gclust.clu
python3 cluster2fasta.py gclust/unalign_no2_merge.sorted.fa gclust/unalign.gclust.clu gclust/unalign.gclust.fa


# blast
th=40
~/tool/EUPAN-v0.44/bin/eupan rmRedundant blastCluster -c 0.9 \
	-t ${th} \
	gclust/unalign.gclust.fa \
	selfblast \
	~/tool/ncbi-blast/bin
unalign_cluster_fa=selfblast/non-redundant.fa

#######
# Blast nt database and remove contaminations (blastn)
#######
mkdir blastnt
query=$unalign_cluster_fa
th=10
blastn -db ~/data/db/nt_20200708/nt \
        -query ${query} \
        -out blastnt/unalign_clsutered.blastnt \
        -evalue 1e-5 \
        -best_hit_overhang 0.25 \
        -perc_identity 0.5 \
        -max_target_seqs 5 \
        -num_threads $th \
        -outfmt "6 qacc sacc evalue pident length qstart qend sstart send"

python3 getTax.py blastnt/unalign_clsutered.blastnt 2 ~/data/db/taxonomy_20200715/nucl_gb.accession2taxid ~/data/db/taxonomy_20200715/rankedlineage.dmp > blastnt/unalign_clsutered.tax
python3 filter_exclude_contig.py blastnt/unalign_clsutered.blastnt blastnt/uunalign_clsutered.tax 0 0 0 > blastnt/unalign_clsutered.polution
~/tool/faSomeRecords ${query} -exclude blastnt/unalign_clsutered.polution blastnt/unalign_clsutered_clean.fa



#######
# Sequences elongation (for gene prediction, optional)
#######
mkdir elongated
# cat ragoo/*/*_scaffold/ragtag.scaffolds.fasta > elongated/temp.raw_scaffolds.fa
python elongate_seq.py blastnt/unalign_clsutered_clean.fa 5000 elongated/temp.raw_scaffolds.fa elongated/elongated.fa > seq_elongated_table.txt

