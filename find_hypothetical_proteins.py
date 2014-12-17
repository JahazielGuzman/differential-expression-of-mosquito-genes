import csv

# read in the reference genome, and read from it as a tab delimited file
noribo_ref = open("/media/jaxi/differential_expression/noriboaegypti.gff3", 'r')
hp_list = open("/media/jaxi/differential_expression/tests/hp_list.txt", "ab")

aegypti_reader = csv.reader(noribo_ref, delimiter='\t')

# we will add all the hypothetical proteins to a list
s = []
for gene in aegypti_reader:
	if len(gene) == 9:
		if "hypothetical protein" in gene[8]:			
			j = gene[8].split('; ')
			if "Parent" in j[1]:
				s.append([j[0], j[1]])
			elif "Parent" in j[2]:
				s.append([j[0], j[2]])

# count the number of unique genes which are hypothetical proteins
hp_count = 0
parent = ''
for i in s:
	if i[1] != parent:
		parent = i[1]
		hp_count += 1

# write the id of each hypothetical protein to the file hp_list.txt
for m in s:
	hp_list.write(str(m[0][8:]) + '\n')
