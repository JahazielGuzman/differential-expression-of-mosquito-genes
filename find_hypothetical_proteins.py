import csv

noribo_ref = open("/media/jaxi/differential_expression/noriboaegypti.gff3", 'r')
aegypti_reader = csv.reader(noribo_ref, delimiter='\t')

s = []
for gene in aegypti_reader:
	if len(gene) == 9:
		if "hypothetical protein" in gene[8]:
			j = gene[8].split('; ')
			if "Parent" in j[1]:
				s.append([j[0], j[1]])
			elif "Parent" in j[2]:
				s.append([j[0], j[2]])
