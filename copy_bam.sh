# use these statement to identify rRNA sequences from the aegypti.gff3 annotation file
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

sed -ri "s/ID=/gene_id=/g" $ANNOT 

# get the name of each file on a seperate line for stranded and unstranded
FILE1=$(ls $bam_dir1 | egrep "accepted")
FILE2=$(ls $bam_dir2 | egrep "accepted")

# for each stranded bam file i, sort the file using sametools and output sorted.i.bam
# and show progress

for f1 in $FILE1; do samtools sort -@ 11 -n $bam_dir1/$f1 "sorted."${f1%".bam"}; done

# for each unstranded bam file i, sort the usi ng samtools and output sorted.i.bam
# record progressecho "done with stranded file. NOW do unstranded files"

for f2 in $FILE2; do samtools sort -@ 11 -n $bam_dir2/$f2 "sorted."${f2%".bam"}; done

mv sorted*0* stranded/
mv sorted* unstranded/

###################################################################################
############ only run these lines again if you have the files and #################
############  only want to do read summarization ################

#set the paths where stranded and unstranded bam files will be
bam_dir1="/media/jaxi/differential_expression/stranded"
bam_dir2="/media/jaxi/differential_expression/unstranded"
ANNOT="/media/jaxi/differential_expression/noriboaegypti.gff3"

# get the name of each file on a seperate line for stranded and unstranded
FILE1=$(ls $bam_dir1 | egrep "accepted")
FILE2=$(ls $bam_dir2 | egrep "accepted")

#########################################################################

#featureCounts read summarization

for f1 in $FILE1; do featureCounts $bam_dir1/$f1 -a $ANNOT -F -g -f -t 'exon' -O -s 1 -M -T 12 -p -o ${f1%".bam"}".txt"; done

for f2 in $FILE2; do featureCounts $bam_dir2/$f2 -a $ANNOT -F -g -f -t 'exon' -O -s 0 -M -T 12 -p -o ${f2%".bam"}".txt"; done

mv accepted*0* stranded/diff_expf/
mv accepted* unstranded/diff_expf/

###### create a subset of the annotated genome, to do read summarization only for hypothetical proteins ############
###### the hypothetical proteins are labeled as mRNA;s #############################################################

egrep ".*(mRNA).*hypothetical" noriboaegypti.gff3 > hypothetical_aegypti.gff3

ANNOT="/media/jaxi/differential_expression/hypothetical_aegypti.gff3"

################# do read summarization considering only these sequences ############################

FILE1=$(ls $bam_dir1 | egrep "^accepted_hits0[1-5]\.")
FILE2=$(ls $bam_dir2 | egrep "^accepted_hits(_1|[2-4])\.")

for f1 in $FILE1; do featureCounts $bam_dir1/$f1 -a $ANNOT -F -g -f -t 'mRNA' -O -s 1 -M -T 12 -p -o $bam_dir1/"diff_expf/mRNA."${f1%".bam"}".txt"; done

for f2 in $FILE2; do featureCounts $bam_dir2/$f2 -a $ANNOT -F -g -f -t 'mRNA' -O -s 0 -M -T 12 -p -o  $bam_dir2/"diff_expf/mRNA."${f2%".bam"}".txt"; done

##### FeatureCounts does not recognize the gene_id of one of the hypothetical proteins ##############
##### we now must place this gene's id in the count file ############################################

# obtain the id which was not properly detected
MISSING_HP=$(egrep ".*mRNA.1413190.1421790" hypothetical_aegypti.gff3 | cut -f 9 | cut -d ';' -f 1 | cut -d '=' -f 2
)

# we will now insert the id into the counts files
# because it is missing

# the count files are in the diff_expf folder of the stranded and unstranded folder
de_dir1=$bam_dir1/"diff_expf"
de_dir2=$bam_dir2/"diff_expf"

# get all the files which correspond to mRNA counts
FILE1=$(ls $de_dir1 | egrep "mRNA.*txt$")
FILE2=$(ls $de_dir2 | egrep "mRNA.*txt$")

# for each file, insert the missing id into the line it is missing from and write it to a new file
# mRNA.i.txt where i corresponds to file accepted_hits[i].txt
for f1 in $FILE1; do sed -E "s/(^.Supercontig)/"$MISSING_HP"\1/g" $de_dir1/$f1 > $de_dir1/${f1/"accepted_hits"/""}; done
for f2 in $FILE2; do sed -E "s/(^.Supercontig)/"$MISSING_HP"\1/g" $de_dir2/$f2 > $de_dir2/${f2/"accepted_hits"/""}; done

# obtain the bam files for sample 6-10
FILE3=$(ls $bam_dir1 | egrep "^accepted_hits0([6-9]|10)\.")

# output hypothetical protein counts for each bam file
for f3 in $FILE3; do featureCounts $bam_dir1/$f3 -a $ANNOT -F -g -f -t 'mRNA' -O -s 1 -M -T 12 -p -o $bam_dir1/"diff_expf/mRNA."${f3%".bam"}".txt"; done

# add the missing id to the count files
FILE3=$(ls $de_dir1 | egrep "mRNA.*([6-9]|10)\.txt$")
for f3 in $FILE3; do sed -E "s/(^.Supercontig)/"$MISSING_HP"\1/g" $de_dir1/$f3 > $de_dir1""/${f3/"accepted_hits"/""}; done

FILE4=$(ls $bam_dir1 | egrep "^accepted_hits(1|2)[0-9]\.")

for f4 in $FILE4; do featureCounts $bam_dir1/$f4 -a $ANNOT -F -g -f -t 'mRNA' -O -s 1 -M -T 12 -p -o $bam_dir1/"diff_expf/mRNA."${f4%".bam"}".txt"; done

ANNOT="/media/jaxi/differential_expression/noriboaegypti.gff3"

for f4 in $FILE4; do featureCounts $bam_dir1/$f4 -a $ANNOT -F -g -f -t 'exon' -O -s 1 -M -T 12 -p -o $bam_dir1/"diff_expf/"${f4%".bam"}".txt"; done

FILE4=$(ls $de_dir1 | egrep "mRNA.*(1|2)[1-9]\.txt$")

for f4 in $FILE4; do sed -E "s/(^.Supercontig)/"$MISSING_HP"\1/g" $de_dir1/$f4 > $de_dir1""/${f4/"accepted_hits"/""}; done
