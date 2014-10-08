# remove all transcripts from reference genome aegypti.gff3 where the 
# transcript is a ribosomal RNA (rRNA) sequence, in other words its 
# ID or Parent match the id of a known ribosomal sequence
# use the command aegypti_ref.seek(0) to rewind the file

import csv

# read in reference genome
aegypti_ref = open("/data/temp_bam/aegypti.gff3",'rb')

# write to a new file which will be the updated reference genome
noribo_ref = open("/media/jaxi/differential_expression/noriboaegypti.gff3", 'w')

# read in the file that contains the rRNA id's
rRNA = open("/media/jaxi/differential_expression/rRNA_contigs.txt","rb")

# we will read the reference genome into aegypti_reader and
# write lines of interest to aegypti_writer output stream
aegypti_reader = csv.reader(aegypti_ref, delimiter='\t')
aegypti_writer = csv.writer(noribo_ref, delimiter='\t')

# read the rRNA id's and then divide each id
# into an element into a list which will be stored in
# the list rRNA
rRNA = rRNA.read()
rRNA = rRNA.split(" ")

# there is a newline in the last string element
# this removes the newline
rRNA[-1] = rRNA[:-1]

# for each row i in the reference genome that is not
# a comment (i is a list with 9 elements)
# check for each rRNA id j
# whether j appears in row i
# and if it does not, print row i to the file noriboaegypti.gff3
for gene in aegypti_reader:
	if len(gene) == 9: 	
		b = 0
		for rnas in rRNA:
			b = b or (str(rnas) in str(gene[-1])) # ID/Parent fields are in element 9
		if not b:
			aegypti_writer.writerow(gene)
