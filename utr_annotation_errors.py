"""
this file must be run like this: python utr_annotation_errors.py sample_name
where sample_name is the name of some read file whose contents are contained 
in a tab-delimited gzipped text file  stored in a folder named bam_files. 
A directory tests is also assumed to exist and this is where the output will be
stored.
"""

import csv
import gzip
import sys

### read in the name of an rna-seq sample 

if len(sys.argv) < 2:
	exit_string = "please enter the file name of the accepted_hits txt file without extension"
	print exit_string
	sys.exit(1)


read_file = "bam_files/" + sys.argv[1] + ".txt.gz"
read_file = gzip.open(read_file,"rb")

UTRs = open("3_prime_utr_list.txt","rb").read()
UTR_results = open("tests/utr_results_for_" + str(sys.argv[1]) + ".txt", "a")
UTRs = UTRs.split("\n")

rnaseq_reads = csv.reader(read_file, delimiter="\t")

"""
   utr has 4 fields, Its identifier in utr[0], its chromosome in utr[1]
   its starting coordinate in utr[2] and its ending coordinate in utr[3]
   each read has 12 fields, but those of interest are its chromosome in read[2]
   its start position in read[3], and its sequence in read[9].

   for each 3'-UTR utr[p...r] where p and r are integer chromosomal coordinates in the annotation file, 
   we will determine the number of reads read[i...j] corresponding to two scenarios: scenario A
   where p <= i <= j <= r - v and v is an integer variable. And scenario B where
   r <= j, note we do not concern ourselves with i and p in scenario B. After computing the counts
   both these scenarios, the name of each UTR as well as the counter for A and B respectively will be
   output to the file stored in UTR_results in the tests directory.
"""

v = 50
utr_dict = {}
utr_test_string = ""
for utr in UTRs:
	utr = utr.split("\t")
	utr_dict[utr[0]] = [0,0]
	for read in rnaseq_reads:
		if read[2] == utr[1]:
			if int(read[3]) < int(utr[3]):
				read_end = int(read[3]) + len(read[9]) - 1
				if int(utr[2]) <= int(read[3]) and read_end <= (int(utr[3]) - v):
					utr_dict[utr[0]][0] += 1
				elif utr[3] <= read_end:
					utr_dict[utr[0]][1] += 1 
			else:
				break
		elif read[2] < utr[1]:
			continue
		else:
			break
	
	##### print the results to the console as well as to the output file
	print "finished " + str(utr[0]) + " " + str(utr_dict[utr[0]][0]) + " " + str(utr_dict[utr[0]][1])
	utr_test_string = "\t".join(map(str,[utr[0],utr_dict[utr[0]][0],utr_dict[utr[0]][1]])) + "\n"
	UTR_results.write(utr_test_string)
	#### reset the read file to the beginning
	read_file.seek(0)
	rnaseq_reads = csv.reader(read_file, delimiter="\t")
