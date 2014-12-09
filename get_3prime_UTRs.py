# this script will open the noriboaegypti.gff3 annotation file, it will find sequences which are
# labeled as three_prime_utr and will obtain an identifier, the chromosome, and the starting and ending
# chromosomal coordinates for each these sequences and will
# output it into the file 3_prime_utr_list.txt

import csv

# read lines from the annotated genome
noribo_ref = open("noriboaegypti.gff3", 'rb')

# write the 3_prime_utr data to 
UTR_file = open("3_prime_utr_list.txt","w")

# we will read the reference genome into aegypti_reader and
# write lines of interest to aegypti_writer output stream
aegypti_reader = csv.reader(noribo_ref, delimiter='\t')

# for sequence in the annotation, if it is a three prime utrthree_prime_utr:
# print the tab delimited sequence identifier in gene[9], the chromosome 
# in gene[0], and the starting and ending coordinates in gene[3] and gene[4]

utr_string = ""
for gene in aegypti_reader:
	if len(gene) == 9 and gene[2] == "three_prime_utr":
		utr_string += "\t".join(map(str, [gene[-1].split("=")[-1][:-1]
			,gene[0],gene[3],gene[4]])) + "\n"

UTR_file.write(utr_string)
