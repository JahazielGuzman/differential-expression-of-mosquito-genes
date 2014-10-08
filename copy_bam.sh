# use this statement to identify rRNA sequences from the aegypti.gff3 annotation file
filter_rRNA=$(egrep "VectorBase.rRNA" /data/temp_bam/aegypti.gff3 | cut -f 9 |cut -d' ' -f 1 | cut -d'=' -f 2)
RNAs=$(for RNA in $filter_rRNA; do echo $RNA; done)
# extra formatting may be necessary, remove semicolons
echo $RNAs | sed "s/;//g"


#copy the stranded, paired end files 01-05 to the external drive
cp /data/Project_BROOKLYN/accepted_hits_combined/accepted_hits_01-05_compare/accepted_hits* /media/jaxi/differential_expression/stranded/

#copy the unstranded, paired end files 1-5 to the external drive
cp /data/seqfiles_SW/accepted_hits_sw/accepted_hits_1.bam /data/seqfiles_SW/accepted_hits_sw/accepted_hits[2-5].bam /media/jaxi/differential_expression/unstranded/

#set the paths where stranded and unstranded bam files will be
bam_dir1="/media/jaxi/differential_expression/stranded"
bam_dir2="/media/jaxi/differential_expression/unstranded"
ANNOT="/media/jaxi/differential_expression/noriboaegypti.gff3"

# get the name of each file on a seperate line for stranded and unstranded
FILE1=$(ls $bam_dir1 | egrep "accepted")
FILE2=$(ls $bam_dir2 | egrep "accepted")

# for each stranded bam file i, sort the file using sametools and output sorted.i.bam
# and show progress

for f1 in $FILE1; do samtools sort -@ 11 -n $bam_dir1/$f1 "sorted."${f1%".bam"}; done

# for each unstranded bam file i, sort the usi ng samtools and output sorted.i.bam
# record progressecho "done with stranded file. NOW do unstranded files"

for f2 in $FILE2; do samtools sort -@ 11 -n $bam_dir2/$f2 "sorted."${f2%".bam"}; done

#featureCounts read summarization

for f1 in $FILE1; do featureCounts $bam_dir1/$f1 -a $ANNOT -F -g -f -t 'exon' -O -s 1 -M -T 12 -p -o ${f1%".bam"}".txt"; done

for f2 in $FILE2; do featureCounts $bam_dir2/$f2 -a $ANNOT -F -g -f -t 'exon' -O -s 0 -M -T 12 -p -o ${f2%".bam"}".txt"; done
