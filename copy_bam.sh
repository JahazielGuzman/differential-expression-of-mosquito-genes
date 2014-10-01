#copy the stranded, paired end files 01-05 to the external drive
cp /data/Project_BROOKLYN/accepted_hits_combined/accepted_hits_01-05_compare/accepted_hits* /media/jaxi/differential_expression/stranded/

#copy the unstranded, paired end files 1-5 to the external drive
cp /data/seqfiles_SW/accepted_hits_sw/accepted_hits_1.bam /data/seqfiles_SW/accepted_hits_sw/accepted_hits[2-5].bam /media/jaxi/differential_expression/unstranded/

#set the paths where stranded and unstranded bam files will be
bam_dir1="/media/jaxi/differential_expression/stranded"
bam_dir2="/media/jaxi/differential_expression/unstranded"

# get the name of each file on a seperate line for stranded and unstranded
FILE1=$(ls $bam_dir1 | egrep "accepted")
FILE2=$(ls $bam_dir1 | egrep "accepted")

# for each stranded bam file i, sort the file using sametools and output sorted.i.bam
for f1 in $FILE1; do samtools sort -n $bam_dir1/$f1 "sorted."${i%".bam"}; done

# for each unstranded bam file i, sort the using samtools and output sorted.i.bam
for f2 in $FILE2; do samtools sort -n $bam_dir2/$f2 "sorted."${i%".bam"}; done


