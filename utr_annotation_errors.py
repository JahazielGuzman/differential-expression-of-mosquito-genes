import csv
import gzip

if len(sys.argv) < 1:
exit_string = "please enter the file name of the accepted_hits file, the file should be a txt file\n please leave out the file extension"
print exit_string
sys.exit(1)

read_file = sys.argv[0] + ".txt"
read_file = gzip.open(read_file,"rb")
UTRs = open("3_prime_utr_list.txt","rb").read()
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
			if read[3] < utr[3]:
				read_end = read[3] + len(read[9]) - 1
				if utr[2] <= read[3] and read_end <= utr[3] - v:
					utr_dict[utr[0]][0] += 1
				elif utr[3] <= read_end:
					utr_dict[utr[0]][1] += 1 
			else:
				break
		elif read[2] < utr[1]:
			continue
		else:
			break
	rnaseq_reads.seek(0)

for utr in utr_dict:
	utr_test_string += "\t".join(map(str,[utr,utr_dict[utr][0],utr_dict[utr][1]]))
