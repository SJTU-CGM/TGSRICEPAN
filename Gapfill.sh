# needs tools
# TGS-GapCloser

# needs scripts
# gapfill.py

####### For samples for GJ

#######
# Corrected reads
#######
corr_reads=${sample}.corr_reads.fa
th=1
~/tool/TGS-GapCloser_v1.0.1/TGS-GapCloser.sh \
  --scaff ref/msu7/chrs.fa \
  --reads ${corr_reads} \
  --output msu7_${sample}_corr_reads_gapfill \
  --tgstype ont \
  --min_idy 0.3 \
  --min_match 300 \
  --ne \
  --thread ${th} >pipe.log 2>pipe.err

######
# Assembled contigs
######
ass_contigs=${sample}.contigs.fa
th=1
~/tool/TGS-GapCloser_v1.0.1/TGS-GapCloser.sh \
--scaff ref/msu7/chrs.fa \
--reads ${ass_contigs} \
--output msu7_${sample}_ass_contigs_gapfill \
--tgstype ont \
--min_idy 0.9 \
--min_match 1000 \
--ne \
--thread ${th} >pipe.log 2>pipe.err

######
# Merge
######
cat */*updated_scaff_infos > merged.updated_scaff_infos
# chrs.contigs.fa: break all N of chr and to every contig,
# exist in output of TGS-GapCloser, named msu7_${sample}_ass_contigs_gapfill.contig
mv msu7_${sample}_ass_contigs_gapfill.contig ref/msu7/chrs.contigs.fa
cat ref/msu7/chrs.contigs.fa
python gapfill.py merged.updated_scaff_infos min > msu7.chrs.min.fa