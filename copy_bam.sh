cp /data/Project_BROOKLYN/accepted_hits_combined/accepted_hits_01-05_compare/accepted_hits* /media/jaxi/differential_expression/stranded/

cp /data/seqfiles_SW/accepted_hits_sw/accepted_hits_1.bam /data/seqfiles_SW/accepted_hits_sw/accepted_hits[2-5].bam /media/jaxi/differential_expression/unstranded/

bam_dir1="/media/jaxi/differential_expression/stranded"
bam_dir2="/media/jaxi/differential_expression/unstranded"

FILE1=$(ls $bam_dir1 | egrep "accepted")
FILE2=$(ls $bam_dir1 | egrep "accepted")

for f1 in $FILE1; do samtools sort -n $bam_dir1/$f1 "sorted."${i%".bam"}; done

for f2 in $FILE2; do samtools sort -n $bam_dir2/$f2 "sorted."${i%".bam"}; done


