#### aegypti_de.R by Jahaziel Guzman ##############
#### This R script will read the text files exon_counts.txt, 
#### aegypti_metadata.txt and hp_counts.txt
#### both exon_counts.txt and hp_counts.txt contain a
#### table of integer values, where each row corresponds to an exon
#### which is represented as a string, and each column corresponds to a set
#### of reads obtained from a particular sequencing experiment, also known as a sample. 
#### The count value of the ith row and jth column of the table thus represents the number of reads in
#### sample j that match the sequence of gene i.
#### these count values will then be plotted in heatmaps and scatterplots, showing the level of expression

library(dplyr)
library(DESeq)
library(RColorBrewer)
library(gplots)
library(ggplot2)
 
### Read in the table of exon counts into aegypti_table, then read
### the metadata for all the samples into aegypti_filter
aegypti_table <- read.table("tests/exon_counts.txt", row.names = 1, header = TRUE)
aegypti_filter <- read.table("tests/aegypti_metadata.txt")

# use counts for all exons, from sample 01-05 and 1-5
aegypti_samples <- aegypti_table[,1:10]

# now start differential expression, use metadata corresponding to the
# samples 01-05 and 1-5, in particular the time of mRNA extraction and the
# body part from which the mRNA was extracted
aegypti_cds <- newCountDataSet(aegypti_samples, aegypti_filter[1:10,1:2])
aegypti_cds <- estimateSizeFactors(aegypti_cds)

# store the count data in aegypti_counts
aegypti_counts <- counts(aegypti_cds)

# hmcol will hold the range of colors to be displayed in the heatmap
hmcol = colorRampPalette(brewer.pal(9,"OrRd"))(100)

# most_expressed will contain the 100 most expressed genes, this statement
# initializes most_expressed to have 100 rows, this is done so the cbind() function
# will work
most_expressed <- data.frame(row.names(aegypti_counts)[1:100])

# generate 5 pdf files, each pdf file i (where 0 < i <= 5)
# will contain the i*100 most expressed genes across samples 1-5 and 01-05

for (i in 1:5) {
    pdf(file = paste0(100*i,"tests/most expressed genes for each sample"), onefile = TRUE, width=8, height=(15*j))
    for (j in 1:10) {

     	  #sample_hm will refer to the current sample under consideration
     	  sample_hm <- cbind(
       	            ifelse(aegypti_counts[,j] != 0,log2(aegypti_counts[,j]),0), replicate(
                      nrow(aegypti_counts), 0))
	  
     	  colname = colnames(aegypti_counts)[j]
   
	  # order the genes in sample_hm by gene expression from highest to lowest
     	  sample_hm <- sample_hm[order(sample_hm[,1], decreasing=TRUE),]

	  # store the names of the 100 most expressed genes for the current
	  # sample as a column in the most_expressed data.frame
 	  most_expressed[,j] <- row.names(sample_hm)[1:100]
     	  
	  # add a dummy column because heatmap.2 needs at least 2 columns to map
     	  colnames(sample_hm) <- c(colname, "dummy")
	  
	  # generate the heatmap
       	  heatmap.2(sample_hm[1:(100*i),], col=hmcol, trace="none", margin=c(10,6), main="most expressed genes each sample", srtCol=0)  
     
     }
     dev.off()
}

# set the column names for most_expressed
colnames(most_expressed) = colnames(aegypti_counts)

# create an empty vector
m_exp_names <- c()

n = ncol(most_expressed)

# assign all the names belonging to each sample
for (i in 1:n) {
	m_exp_names <- c(m_exp_names,most_expressed[,i])
}

# then get the unique ones
m_exp_names <- unique(m_exp_names)

# this is the number of unique gene id's we have
n = length(m_exp_names)

# this matrix will hold the comparison results
# where each row i represents a gene, and each column j represents
# a sample, the ijth entry is the rank of gene i in sample j
most_exp_table <- matrix(0,ncol=10,nrow=n)

for (i in 1:n) {
	m = most_exp_table[i,]
	m <- sapply(most_expressed,function(x) which(m_exp_names[i]==x))	
	m[is.na(as.numeric(m))] = 0
	most_exp_table[i,] = as.numeric(m)

}

row.names(most_exp_table) = m_exp_names
colnames(most_exp_table) = colnames(most_expressed)

write.table(most_exp_table, "tests/comparison_table.txt")

for (i in ncol(most_exp_table)) {
	m = most_exp_table[,i]
 	m = ifelse(m != 0, log2(m), 0)
	most_exp_table[,i] = m
}

split = n / 4

pdf(file = "tests/ranking of genes weighted by expression", onefile = TRUE, width=8, height=50)
heatmap.2(most_exp_table[1:mid,], col=hmcol, trace="none", margin=c(10,6), main="comparing rank of genes weighted by expression pt. 1")
heatmap.2(most_exp_table[(mid+1):n,], col=hmcol, trace="none", margin=c(10,6), main="comparing rank of genes weighted by expression pt. 2")
dev.off()

aegypti_counts_norm <- counts(aegypti_cds, normalized=TRUE)
aegypti_counts_norm <- aegypti_counts_norm[,c("accepted_1", "accepted_01")]
aegypti_rep_1 <- aegypti_counts_norm[order(rowMeans(aegypti_counts_norm)),]
aegypti_rep_1 <- aegypti_rep_1[rowMeans(aegypti_rep_1) > 0,]

aegypti_rep_1[,1] <- ifelse(aegypti_rep_1[,1] != 0, log2(aegypti_rep_1[,1]), 0)
aegypti_rep_1[,2] <- ifelse(aegypti_rep_1[,2] != 0, log2(aegypti_rep_1[,2]), 0)

pdf(file="tests/replicate_1_comparison", width=8, height=50)
heatmap.2(aegypti_rep_1[(aegypti_rep_1[,1] > 0 & aegypti_rep_1[,2] > 0),], col=hmcol, trace="none", margin=c(10,6))
dev.off()

# sample 4 vs sample 04

aegypti_counts_norm <- counts(aegypti_cds, normalized=TRUE)
aegypti_counts_norm <- aegypti_counts_norm[,c("accepted_4", "accepted_04")]
aegypti_rep_1 <- aegypti_counts_norm[order(rowMeans(aegypti_counts_norm)),]
aegypti_rep_1 <- aegypti_rep_1[rowMeans(aegypti_rep_1) > 0,]

aegypti_rep_1[,1] <- ifelse(aegypti_rep_1[,1] != 0, log2(aegypti_rep_1[,1]), 0)
aegypti_rep_1[,2] <- ifelse(aegypti_rep_1[,2] != 0, log2(aegypti_rep_1[,2]), 0)

pdf(file="tests/replicate_4_comparison", width=8, height=100)
heatmap.2(aegypti_rep_1[(aegypti_rep_1[,1] > 0 & aegypti_rep_1[,2] > 0),], col=hmcol, trace="none", margin=c(10,6))
dev.off()

# do differential expression for 24 hour replicates and 72 hour replicates

aegypti_24x72 <- aegypti_filter$time == 24 | aegypti_filter$time == 72
aegypti_24x72 <- aegypti_24x72 * (aegypti_filter$bodyPart == "midgut" & aegypti_filter$treated == "no")
aegypti_24x72 <- as.logical(aegypti_24x72)

aegypti_exp_24x72 <- aegypti_table[,aegypti_24x72]
aegypti_hours <- aegypti_filter[aegypti_24x72,]$time

aegypti_cds <- newCountDataSet(aegypti_exp_24x72, aegypti_hours)
aegypti_cds <- estimateSizeFactors(aegypti_cds)
aegypti_cds <- estimateDispersions(aegypti_cds, fitType="local")
aegypti_res <- nbinomTest(aegypti_cds, "24", "72")

aegypti_resSig <- aegypti_res[order(aegypti_res$pval),]
aegypti_resSig <- counts(aegypti_cds, normalized=TRUE)[aegypti_resSig$id,]
for (i in 1:ncol(aegypti_resSig)) {
	m = aegypti_resSig[,i]
 	m = ifelse(m != 0, log2(m), 0)
	aegypti_resSig[,i] = m
}

pdf(file="tests/most differentially expressed genes for samples 2/02 and 4/04", width =8, height = 15)
heatmap.2(aegypti_resSig[1:100,], col=hmcol, trace="none", margin=c(10,6))
dev.off()

######### show expression for hypothetical proteins #################################

hyp_protein <- read.table("tests/hp_list.txt")
hyp_protein <- as.character(hyp_protein$V1)
aegypti_resSig <- data.frame(counts(aegypti_cds, normalized=TRUE))[hyp_protein,]
aegypti_resSig <- aegypti_resSig[order(rowmeans(aegypti_resSig), decreasing=TRUE),]
for (i in 1:ncol(aegypti_resSig)) {
	m = aegypti_resSig[,i]
 	m = ifelse(m != 0, log2(m), 0)
	aegypti_resSig[,i] = m
}

pdf(file="100 hypothetical proteins", width =8, height = 15)
heatmap.2(as.matrix(aegypti_resSig), col=hmcol, trace="none", margin=c(10,6))
dev.off()

for (i in 1:5) {
	mRNA_01 <- read.table("stranded/diff_expf/mRNA.01.txt", header=TRUE, row.names=1, skip=1, stringsAsFactors=FALSE)

	names(mRNA_01) <- "accepted_hits01"
	
	qplot(data=mRNA_01, x = row.names(mRNA_01), 
	  y = ifelse(accepted_hits01 == 0, 0, log2(accepted_hits01)), 
            main = "hypothetical protein counts in sample 01") + 
              ylim(0,20) + xlab("hypothetical proteins") + ylab("log2 of counts")

	ggsave("tests/hypothetical_proteins_sample_01.pdf")
}

###### show a scatter plot for the counts of hypothetical proteins in sample 01 ##################################
mRNA_01 <- read.table("stranded/diff_expf/mRNA.01.txt", header=TRUE, row.names=1, skip=1, stringsAsFactors=FALSE)

names(mRNA_01) <- "accepted_hits01"

qplot(data=mRNA_01, x = row.names(mRNA_01), 
     y = ifelse(accepted_hits01 == 0, 0, log2(accepted_hits01)), 
          main = "hypothetical protein counts in sample 01") + 
               ylim(0,20) + xlab("hypothetical proteins") + ylab("log2 of counts")

ggsave("tests/hypothetical_proteins_sample_01.pdf")

#### TRYING TO ORDER THE CHROMOSOME NAMES #####
i[order(substring(i$Chr,13,nchar(i$Chr))),]

###### show a scatter plot for the counts of hypothetical proteins in sample 03 ##################################
mRNA_03 <- read.table("stranded/diff_expf/mRNA.03.txt", header=TRUE, row.names=1, skip=1, stringsAsFactors=FALSE)

names(mRNA_03)[6] <- "accepted_hits03"

qplot(data=mRNA_03, x = row.names(mRNA_03), 
     y = ifelse(accepted_hits03 == 0, 0, log2(accepted_hits03)), 
          main="hypothetical protein counts in sample 03") + 
               ylim(0,20) + xlab("hypothetical proteins") + ylab("log2 of counts")

ggsave("tests/hypothetical_proteins_sample_02.pdf")


mRNA_02 <- read.table("stranded/diff_expf/mRNA.02.txt", header=TRUE, row.names=1, skip=1, stringsAsFactors=FALSE)

names(mRNA_02)[6] <- "accepted_hits02"

qplot(data=mRNA_02, x = row.names(mRNA_02), 
     y = ifelse(accepted_hits02 == 0, 0, log2(accepted_hits02)), 
          main="hypothetical protein counts in sample 02") +
               ylim(0,20) + xlab("hypothetical proteins") + ylab("log2 of counts")

ggsave("tests/hypothetical_proteins_sample_02.pdf")

###### plot for 04 #################################################################
mRNA_04 <- read.table("stranded/diff_expf/mRNA.04.txt", header=TRUE, row.names=1, skip=1, stringsAsFactors=FALSE)

names(mRNA_04)[6] <- "accepted_hits04"

qplot(data=mRNA_04, x = row.names(mRNA_04),
     y = ifelse(accepted_hits04 == 0, 0, log2(accepted_hits04)),
          main="hypothetical protein counts in sample 04") +
               ylim(0,20) + xlab("hypothetical proteins") + ylab("log2 of counts")

ggsave("tests/hypothetical_proteins_sample_04.pdf")

mRNA_05 <- read.table("stranded/diff_expf/mRNA.05.txt", header=TRUE, row.names=1, skip=1, stringsAsFactors=FALSE)

names(mRNA_05)[6] <- "accepted_hits05"

qplot(data=mRNA_05, x = row.names(mRNA_05),
     y = ifelse(accepted_hits05 == 0, 0, log2(accepted_hits05)),
          main="hypothetical protein counts in sample 05") +
               ylim(0,20) + xlab("hypothetical proteins") + ylab("log2 of counts")

ggsave("tests/hypothetical_proteins_sample_05.pdf")

rm(mRNA_01, mRNA_02, mRNA_03, mRNA_04, mRNA_05)

aegyptiHpCounts <- read.table("tests/hp_counts.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)

aegypti_cds <- newCountDataSet(aegyptiHpCounts, aegypti_filter[c(1:4, 6:10),1:2])
aegypti_cds <- estimateSizeFactors(aegypti_cds)

# use the normalized counts to draw heatmaps
aegypti_counts <- counts(aegypti_cds, normalized=TRUE)
# sort in descending order by most expressed genes
aegypti_hp_map <- aegypti_counts[order(rowMeans(aegypti_counts), decreasing=TRUE),]
# obtain the log2 of counts for each sample
cols <- ncol(aegypti_hp_map)
for (i in 1:cols) {
aegypti_hp_map[,i] <- ifelse(aegypti_hp_map[,i] == 0, 0, log2(aegypti_hp_map[,i]))
}
# generate heatmap for 100 most expressed hypothetical proteins for sample 1-4
pdf(file="tests/expression_hypothetical_protein_samples_1thru4", width =8, height = 15)
heatmap.2(aegypti_hp_map[1:100,1:4], col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, 1-4")
dev.off()
# do the same for samples 01-05
pdf(file="tests/expression_hypothetical_protein_samples_01thru04", width =8, height = 15)
heatmap.2(aegypti_hp_map[1:100,5:9], col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, 01-04")
dev.off()

aegypti_counts <- counts(aegypti_cds)
# sort in descending order by most expressed genes
aegypti_hp_map <- aegypti_counts[order(rowMeans(aegypti_counts), decreasing=TRUE),]
# obtain the log2 of counts for each sample
cols <- ncol(aegypti_hp_map)
for (i in 1:cols) {
aegypti_hp_map[,i] <- ifelse(aegypti_hp_map[,i] == 0, 0, log2(aegypti_hp_map[,i]))
}

toPlot <- aegypti_hp_map[1:50,]
# generate heatmap for 100 most expressed hypothetical proteins for sample 1-4
pdf(file="tests/expression_hypothetical_protein_samples_strandedunstranded", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, stranded & unstranded", dendogram="both")
dev.off()

toPlot <- aegypti_hp_map[50:100,]
pdf(file="tests/next_50_most_expressed_hypothetical_protein_samples_strandedunstranded", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, stranded & unstranded", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- aegypti_hp_map[101:150,]
pdf(file="tests/101-150_most_expressed_hypothetical_protein_samples_strandedunstranded", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, stranded & unstranded", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- aegypti_hp_map[151:200,]
pdf(file="tests/151-200_most_expressed_hypothetical_protein_samples_strandedunstranded", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, stranded & unstranded", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

for (i in mRNA_files) {
	mRNA_09 <- read.table(paste0("stranded/diff_expf/",i), header=TRUE, row.names=1, skip=1, stringsAsFactors=FALSE)
	qplot(data=mRNA_09, x = row.names(mRNA_09),
     	y = ifelse(mRNA_09[,6] == 0, 0, log2(mRNA_09[,6])),
          main=paste0("hypothetical protein counts in",i)) +
               ylim(0,20) + xlab("hypothetical proteins") + ylab("log2 of counts")

	ggsave(paste0("tests/plots_for_", i, ".pdf"))
}

###############################################################################
### now show the expression for hypothetical proteins in the cuticle samples###
###############################################################################

aegyptiHpCounts <- read.table("tests/hp_counts.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)

aegypti_cds <- newCountDataSet(aegyptiHpCounts[,10:14], aegypti_filter[11:15,1:2])
aegypti_cds <- estimateSizeFactors(aegypti_cds)

aegypti_counts <- counts(aegypti_cds)
# sort in descending order by most expressed genes
aegypti_hp_map <- aegypti_counts[order(rowMeans(aegypti_counts), decreasing=TRUE),]
# obtain the log2 of counts for each sample
cols <- ncol(aegypti_hp_map)
for (i in 1:cols) {
aegypti_hp_map[,i] <- ifelse(aegypti_hp_map[,i] == 0, 0, log2(aegypti_hp_map[,i]))
}

toPlot <- rbind(aegypti_hp_map[1:50,],c(20,20,20,20,20))
# generate heatmap for 100 most expressed hypothetical proteins for sample 06-010
pdf(file="tests/expression_hypothetical_protein_samples_cuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[50:100,],c(20,20,20,20,20))
pdf(file="tests/next_50_most_expressed_hypothetical_protein_samples_cuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[101:150,],c(20,20,20,20,20))
pdf(file="tests/101-150_most_expressed_hypothetical_protein_samples_cuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[151:200,],c(20,20,20,20,20))
pdf(file="tests/151-200_most_expressed_hypothetical_protein_samples_cuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hypothetical proteins, cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

###############################################################################
###############################################################################
#### Now compare midgut and cuticle samples at the same stages of larval development
###############################################################################

aegyptiCompare <- read.table("tests/exon_counts.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)

aegypti_cds <- newCountDataSet(aegyptiCompare[, 6:15], aegypti_filter[6:15, 1:2])
aegypti_cds <- estimateSizeFactors(aegypti_cds)

aegypti_counts <- counts(aegypti_cds)

# order rows by most significantly expressed
aegypti_hp_map <- aegypti_counts[order(rowMeans(aegypti_counts), decreasing=TRUE),]

# obtain the log2 of counts for each sample
cols <- ncol(aegypti_hp_map)

for (i in 1:cols) {
aegypti_hp_map[,i] <- ifelse(aegypti_hp_map[,i] == 0, 0, log2(aegypti_hp_map[,i]))
}

toPlot <- rbind(aegypti_hp_map[1:50,],c(20,20,20,20,20))
# generate heatmap for 100 most expressed exons for sample 06-010
pdf(file="tests/most_expressed_exon_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed exons, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[50:100,],c(20,20,20,20,20))
pdf(file="tests/next_50_most_expressed_exon_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed exons, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[101:150,],c(20,20,20,20,20))
pdf(file="tests/101-150_most_expressed_exon_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed exons, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[151:200,],c(20,20,20,20,20))
pdf(file="tests/151-200_most_expressed_exon_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed exons, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

###############################################################################
###############################################################################
##### hypothetical protein sequences, midgut vs cuticle #######################

aegyptiCompare <- read.table("tests/hp_counts.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)

aegypti_cds <- newCountDataSet(aegyptiCompare[, 5:14], aegypti_filter[5:14, 1:2])
aegypti_cds <- estimateSizeFactors(aegypti_cds)

aegypti_counts <- counts(aegypti_cds)

# order rows by most significantly expressed
aegypti_hp_map <- aegypti_counts[order(rowMeans(aegypti_counts), decreasing=TRUE),]

# obtain the log2 of counts for each sample
cols <- ncol(aegypti_hp_map)

for (i in 1:cols) {
aegypti_hp_map[,i] <- ifelse(aegypti_hp_map[,i] == 0, 0, log2(aegypti_hp_map[,i]))
}

toPlot <- rbind(aegypti_hp_map[1:50,],c(20,20,20,20,20))
# generate heatmaps for 100 most expressed hypothetical proteins for sample 06-010
pdf(file="tests/most_expressed_hyp.pro_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[50:100,],c(20,20,20,20,20))
pdf(file="tests/next_50_most_expressed_hyp.pro_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[101:150,],c(20,20,20,20,20))
pdf(file="tests/101-150_most_expressed_hyp.pro_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[151:200,],c(20,20,20,20,20))
pdf(file="tests/151-200_most_expressed_hyp.pro_midgutvscuticle", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

###### hypothetical protein sequences in the methoxyfenozide treated samples #####

aegyptiCompare <- read.table("tests/hp_counts.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)

aegypti_cds <- newCountDataSet(aegyptiCompare[, 15:18], aegypti_filter[16:19, 1:2])
aegypti_cds <- estimateSizeFactors(aegypti_cds)

aegypti_counts <- counts(aegypti_cds)

# order rows by most significantly expressed
aegypti_hp_map <- aegypti_counts[order(rowMeans(aegypti_counts), decreasing=TRUE),]

# obtain the log2 of counts for each sample
cols <- ncol(aegypti_hp_map)

for (i in 1:cols) {
aegypti_hp_map[,i] <- ifelse(aegypti_hp_map[,i] == 0, 0, log2(aegypti_hp_map[,i]))
}

toPlot <- rbind(aegypti_hp_map[1:50,],c(20,20,20,20))

# generate heatmaps for 100 most expressed hypothetical proteins for sample 06-010
pdf(file="tests/most_expressed_hyp.pro_midgutvscuticle_methoxy", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[50:100,],c(20,20,20,20))
pdf(file="tests/next_50_most_expressed_hyp.pro_midgutvscuticle_methoxy", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[101:150,],c(20,20,20,20))
pdf(file="tests/101-150_most_expressed_hyp.pro_midgutvscuticle_methoxy", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

toPlot <- rbind(aegypti_hp_map[151:200,],c(20,20,20,20))
pdf(file="tests/151-200_most_expressed_hyp.pro_midgutvscuticle_methoxy", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most expressed hp, midgut vs cuticle", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()

### test for differential expression between midgut samples of untreated and methoxyfenozide treated samples


aegypti_cds <- newCountDataSet(aegyptiCompare[,c(3,7,15,16)], aegypti_filter[c(3,8,16,17), 3])
aegypti_cds <- estimateSizeFactors(aegypti_cds)
aegypti_cds <- estimateDispersions(aegypti_cds, fitType="local")

aegypti_res <- nbinomTest(aegypti_cds, "no", "methoxyfenozide")
lowest_pval_exons <- aegypti_res[order(aegypti_res$pval),]

aegypti_hp_map <- counts(aegypti_cds)[order(aegypti_res$pval),]

# obtain the log2 of counts for each sample
cols <- ncol(aegypti_hp_map)

for (i in 1:cols) {
aegypti_hp_map[,i] <- ifelse(aegypti_hp_map[,i] == 0, 0, log2(aegypti_hp_map[,i]))
}

toPlot <- rbind(aegypti_hp_map[1:50,],c(20,20,20,20))

# generate heatmaps for 100 most expressed hypothetical proteins for sample 06-010
pdf(file="tests/most_significantly_diff_exp_midgut_48_methoxyfenozidevsuntreated", width =8, height = 15)
heatmap.2(toPlot, col=hmcol, trace="none", margin=c(10,6), main="most sig diff exp hp, treated vs untreated", Rowv=FALSE, Colv=FALSE, dendrogram="none")
dev.off()
