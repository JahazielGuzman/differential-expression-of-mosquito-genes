Differential Expression
=======================
These Python, R and BASH scripts have been used to process BAM files containing Aedes Aegypti reads aligned to a reference genome obtained from VectorBase. Many of the scripts in this repository cannot be properly used on another computer, since many of BAM files are too large to upload.

These BAM files each contain a set of "reads" which are digitized RNA sequences output from an RNA sequencer. These srquences are converted to cDNA sequences which are strings in the alphabet {A,T,G,C}*. each sequence corresponds to a fragment of an mRNA molecule extracted from an organism under a particular experimental condition. The sequences in the BAM files have been mapped to an "reference genome" and thus contain the corresponding numeric positions (also called "chromosomal coordinates") to which they map in the reference.

copy_bam.sh: moves the BAM files stored on the local computer to an external hard drive labeled JAXI. The files are then processed with a program called featureCounts. featureCounts takes every read in each BAM file supplied to it and checks whether that read overlaps with any annotated sequences in an "annotated genome" which contains several metadata for each known sequence in the organisms genome including the starting and ending chromosomal coordinates as well as its length.
The main output of this script is a set of read summarization files each containing a table, where rows correspond to sequence identifiers,
and the values of interest are in a column corresponding to the name of the supplied bam file. The values in this column are integer count values.

get_counts.R:
------------
Inputs: Opens the read summarization files output from copy_bam.sh.

Outputs: Merges the counts for each sequence across bam files into a matrix of count values and 
outputs the matrix to the files exon_counts.txt (for exon sequences) and hp_counts.txt (for hypothetical protein sequences).

filter_rRNA_genome.py:
----------------------
Inputs: opens the annotated genome file aegypti.gff3 				
Outputs: An edited annotated genome file noriboaegypti.gff3 which has rRNA sequence metadata removed. 

find_hypothetical_proteins.py:									
------------------------------
Inputs: opens the annotated genome file noriboaegypti.gff3 				
Outputs: An edited annotated genome file hypothetical_proteins.gff3 which has rRNA sequence metadata removed. Also outputs a list of hypothetical protein identifiers and writes them to the file hp_list.txt.

get_3prime_UTRs.py:											
-------------------
Inputs: opens the annotated genome file noriboaegypti.gff3
Outputs: prints the identifiers and metadata for 3’-UTR sequences to the file 3_prime_utr_list.txt

aegypti_de.R:												
-------------
Inputs: Opens the files noriboaegypti.gff3, aegypti_metadata.txt, exon_counts.txt, and hp_counts.txt.											Outputs: Produces several heat-maps and scatter plots displaying the log¬2 transformed count data for certain sequences of interest. Also produces the file comparison_table.txt which displays the rank of highly expressed sequences across several samples.

get_bam.sh:													
-----------
Inputs: accesses several read files in the bam_files folder. 					Outputs: Outputs the contents of each read file to a plain text file and compresses it using pigz (a multi-threaded implementation of gzip).

utr_annotation_errors.py: 											
-------------------------
Inputs: Opens the file 3_prime_utr_list.txt and the read file the users supplies as an argument to the program. 						Outputs: Determines the number of reads which are prefixes for each 3’-UTR and number of reads containing prefixes that are also suffixes of the 3’-UTR.
