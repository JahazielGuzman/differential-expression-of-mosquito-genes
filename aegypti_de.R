
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

pdf(file="most expressed for experiments 1 - 5", onefile=TRUE)

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
aegypti_table <- read.table("exon_counts.txt", row.names = 1, header = TRUE)
aegypti_filter <- read.table("aegypti_metadata.txt")

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
    pdf(file = paste0(100*j," most expressed genes for each sample"), onefile = TRUE, width=8, height=(15*j))
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

write.table(most_exp_table, "comparison_table.txt")

for (i in ncol(most_exp_table)) {
	m = most_exp_table[,i]
 	m = ifelse(m != 0, log2(m), 0)
	most_exp_table[,i] = m
}

split = n / 4

pdf(file = "ranking of genes weighted by expression", onefile = TRUE, width=8, height=50)
heatmap.2(most_exp_table[1:mid,], col=hmcol, trace="none", margin=c(10,6), main="comparing rank of genes weighted by expression pt. 1")
heatmap.2(most_exp_table[(mid+1):n,], col=hmcol, trace="none", margin=c(10,6), main="comparing rank of genes weighted by expression pt. 2")
dev.off()

aegypti_counts_norm <- counts(aegypti_cds, normalized=TRUE)
aegypti_counts_norm <- aegypti_counts_norm[,c("accepted_1", "accepted_01")]
aegypti_rep_1 <- aegypti_counts_norm[order(rowMeans(aegypti_counts_norm)),]
aegypti_rep_1 <- aegypti_rep_1[rowMeans(aegypti_rep_1) > 0,]

aegypti_rep_1[,1] <- ifelse(aegypti_rep_1[,1] != 0, log2(aegypti_rep_1[,1]), 0)
aegypti_rep_1[,2] <- ifelse(aegypti_rep_1[,2] != 0, log2(aegypti_rep_1[,2]), 0)

pdf(file="replicate_1_comparison", width=8, height=50)
heatmap.2(aegypti_rep_1[(aegypti_rep_1[,1] > 0 & aegypti_rep_1[,2] > 0),], col=hmcol, trace="none", margin=c(10,6))
dev.off()

# sample 4 vs sample 04

aegypti_counts_norm <- counts(aegypti_cds, normalized=TRUE)
aegypti_counts_norm <- aegypti_counts_norm[,c("accepted_4", "accepted_04")]
aegypti_rep_1 <- aegypti_counts_norm[order(rowMeans(aegypti_counts_norm)),]
aegypti_rep_1 <- aegypti_rep_1[rowMeans(aegypti_rep_1) > 0,]

aegypti_rep_1[,1] <- ifelse(aegypti_rep_1[,1] != 0, log2(aegypti_rep_1[,1]), 0)
aegypti_rep_1[,2] <- ifelse(aegypti_rep_1[,2] != 0, log2(aegypti_rep_1[,2]), 0)

pdf(file="replicate_4_comparison", width=8, height=100)
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
