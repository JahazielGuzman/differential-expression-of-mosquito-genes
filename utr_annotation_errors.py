import csv
import gzip
import sys

if len(sys.argv) < 2:
	exit_string = "please enter the file name of the accepted_hits txt file without extension"
	print exit_string
	sys.exit(1)

read_file = "bam_files/" + sys.argv[1] + ".txt.gz"
read_file = gzip.open(read_file,"rb")
UTRs = open("3_prime_utr_list.txt","rb").read()
UTR_results = open("tests/utr_results_for_" + str(sys.argv[0]) + ".txt", "w")
UTRs = UTRs.split("\n")

rnaseq_reads = csv.reader(read_file, delimiter="\t")

"""
   utr has 4 fields, Its identifier in utr[0], its chromosome in utr[1]
   its starting coordinate in utr[2] and its ending coordinate in utr[3]
   each read has 12 fields, but those of interest are its chromosome in read[2]
   its start position in read[3], and its sequence in read[9]
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
	print "finished " + str(utr[0])
	read_file.seek(0)
	rnaseq_reads = csv.reader(read_file, delimiter="\t")

for utr in utr_dict:
	utr_test_string += "\t".join(map(str,[utr,utr_dict[utr][0],utr_dict[utr][1]])) + "\n"

UTR_results.write(utr_test_string)
