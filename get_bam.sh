# this file assumes that the current directory contains subdirectories named stranded which contains bam files for reads
# obtained through a stranded rna-seq experiment and a subdirectory named unstranded which contains bam files with reads
# that are unstranded. Also assumes there is a directory named bam_files to store the txt converted bam files

bam_dir1="stranded"
bam_dir2="unstranded"

# get the names of the files to view using samtools
BAM_FILES=$(ls $bam_dir1/accepted_hits*)
BAM_FILES=$BAM_FILES" "$(ls $bam_dir2/accepted_hits*)

# output each accepted_hits bam file to a text file with the same name
for bf in $BAM_FILES; do if expr ${bf:0:1} = "u"; then bf=${bf/"unstranded/"/""}; else bf=${bf/"stranded/"/""}; fi;  bf=${bf%".bam"}".txt"; samtools view $bf "bam_files/"$bf; bf="bam_files"/$bf; gzip $bf; done

bams_2_gzip=$(ls bam_files)
for b in $bams_2_gzip; do pigz/pigz "bam_files/"$b; done
