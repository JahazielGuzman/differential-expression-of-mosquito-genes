library(dplyr)
library(DESeq)
library(RColorBrewer)
library(gplots)
### these two commands will read in the count table and
### then read in the metadata for the count table respectively
aegypti_table <- read.table("exon_counts.txt", row.names = 1, header = TRUE)
aegypti_filter <- read.table("aegypti_metadata.txt")
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
plotDispEsts(aegypti_res, main = "avg gene expression in 24 hour and 48 hour midgut samples")
plotMA(aegypti_cds, main = "log2 fold change from 24 hours to 48 hours")
hist(aegypti_res$pval, breaks = 100, col="skyblue", border="slateblue", main="p-values for conditions 24 vs 48")
head(arrange(aegypti_res, foldChange, pval))
###### randomly select genes to display in a heat-map#############

significant_genes <- sample(row.names(aegypti_cds), 30, replace=TRUE)
aegypti_resSig <- aegypti_cds[significant_genes,]

hmcol = colorRampPalette(brewer.pal(9,"OrRd"))(100)
heatmap.2(counts(aegypti_resSig, normalized=TRUE), col=hmcol, trace="none", margin=c(10,6))

significant_genes <- arrange(aegypti_res, pval)[1:30]
aegypti_resSig <- aegypti_cds[significant_genes,]

heatmap.2(counts(aegypti_cds[aegypti_resSig, normalized=TRUE), col=hmcol, trace="none", margin=c(10,6))

significant_genes <- arrange(aegypti_res, log2FoldChange)[1:30]
aegypti_resSig <- aegypti_cds[significant_genes,]

heatmap.2(counts(aegypti_, normalized=TRUE), col=hmcol, trace="none", margin=c(10,6))
