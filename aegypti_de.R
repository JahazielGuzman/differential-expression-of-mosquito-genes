library(dplyr)
library(DESeq)
library(RColorBrewer)
library(gplots)
library(ggplot2)


### these two commands will read in the count table and
### then read in the metadata for the count table respectively
aegypti_table <- read.table("tests/exon_counts.txt", row.names = 1, header = TRUE)
aegypti_filter <- read.table("tests/aegypti_metadata.txt")

# create logical vector to use only metadata for 24 hour and 48 hour
test_hours <- aegypti_filter$time == 24 | aegypti_filter$time == 48
test_hours <- as.logical(test_hours * (aegypti_filter$bodyPart == "midgut" & aegypti_filter$treated == "no"))
test_conditions <- aegypti_filter$time[test_hours]

# only use replicate samples from 24 hour and 48 hour midgut

# experimental conditions
aegypti_24x48 <- aegypti_table[,test_hours]

# now start differential expression
aegypti_cds <- newCountDataSet(aegypti_24x48, test_conditions)
aegypti_cds <- estimateSizeFactors(aegypti_cds)
sizeFactors(aegypti_cds)
aegypti_cds <- estimateDispersions(aegypti_cds, fitType="local")
aegypti_res <- nbinomTest(aegypti_cds, "24", "48")

# y axis is dispersion for each gene, genes are on the x axis, sorted
# by their average expression

plotDispEsts(aegypti_cds, main = "avg gene expression in 24 hour and 48 hour midgut samples")
plotMA(aegypti_res, main = "log2 fold change from 24 hours to 48 hours")
hist(aegypti_res$pval, breaks = 100, col="skyblue", border="slateblue", main="p-values for conditions 24 vs 48")
head(arrange(aegypti_res, foldChange, pval))

###### randomly select genes to display in a heat-map#############

significant_genes <- sample(row.names(aegypti_cds), 30, replace=FALSE)
aegypti_resSig <- aegypti_cds[significant_genes,]
hmcol = colorRampPalette(brewer.pal(9,"OrRd"))(200)

pdf(file="tests/most expressed for experiments 1 - 5", onefile=TRUE)

## print random list of genes
for (i in 1:10)
heatmap.2(counts(aegypti_resSig, normalized=TRUE)[,1:2], col=hmcol, trace="none", margin=c(10,6), main=paste0"most expressed genes for ")

dev.off()


######################################################################
######################################################################
######## map the most expressed genes for samples 01-05 and 1-5 ##########

library(dplyr)
library(DESeq)
library(RColorBrewer)
library(gplots)


### these two commands will read in the count table and
### then read in the metadata for the count table respectively
aegypti_table <- read.table("tests/exon_counts.txt", row.names = 1, header = TRUE)
aegypti_filter <- read.table("tests/aegypti_metadata.txt")

# only use replicate samples from 24 hour and 48 hour midgut

# experimental conditions
aegypti_24x48 <- aegypti_table[,1:10]

# now start differential expression
aegypti_cds <- newCountDataSet(aegypti_24x48, aegypti_filter[1:10,1:2])
aegypti_cds <- estimateSizeFactors(aegypti_cds)

# get the counts from the aegypti table and assign it to variable
# aegypti_counts, this is done for readability
aegypti_counts <- counts(aegypti_cds)

# heatmap will have varying values of red and orange
hmcol = colorRampPalette(brewer.pal(9,"OrRd"))(100)
most_expressed <- data.frame(row.names(aegypti_counts)[1:100])

# generate 5 pdf files, each pdf file i (where i > 0) 
# contains the i*100 most expressed genes across samples 1-5 and 01-05

for (j in 1:5) {
    pdf(file = paste0(100*j,"tests/most expressed genes for each sample"), onefile = TRUE, width=8, height=(15*j))
    for (i in 1:10) {

     	  #sample_hm will refer to the current sample under consideration
     	  sample_hm <- cbind(
       	            ifelse(aegypti_counts[,i] != 0,log2(aegypti_counts[,i]),0), replicate(
                      nrow(aegypti_counts), 0))
	  
     	  colname = colnames(aegypti_counts)[i]
   
	  # order the genes in sample_hm by gene expression from highest to lowest
     	  sample_hm <- sample_hm[order(sample_hm[,1], decreasing=TRUE),]

	  # store the names of the 100 most expressed genes for the current
	  # sample as a column in the most_expressed data.frame
 	  most_expressed[,i] <- row.names(sample_hm)[1:100]
     	  
	  # add a dummy column because heatmap.2 needs at least 2 columns to map
     	  colnames(sample_hm) <- c(colname, "dummy")
	  
	  # generate the heatmap
       	  heatmap.2(sample_hm[1:(100*j),], col=hmcol, trace="none", margin=c(10,6), main="most expressed genes each sample", srtCol=0)  
     
     }
     dev.off() 
}

# set the colum names for most_expressed
colnames(most_expressed)[i] = colnames(aegypti_counts)

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
