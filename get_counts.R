library(dplyr)

# this list will hold the data.table for each read count file produced
# from feature counts
count_files <- list()

# 10 stranded sequence files, 5 unstranded sequence files
strand_num = 1:10
unstrand_num = 1:5


# get the data from the stranded files first
for (i in strand_num) {	
     
     count_files[[i]] <- read.table(paste0("/media/jaxi/differential_expression/stranded/diff_expf/accepted_hits0",i,".txt"), header=T)
     # rename the column containing actual integer counts
     names(count_files[[i]])[7] <- "s_counts"
     # use column 1 containing the exon identifier and column 7 containing the read counts
     count_files[[i]] <- count_files[[i]][,c(1,7)]
     # each exon and its counts should only correspond to one row, make sure this is the case
     count_files[[i]] <- count_files[[i]] %>% group_by(Geneid) %>% summarize(s_counts = sum(s_counts))
}

for (i in unstrand_num) {	

     count_files[[i+10]] <- read.table(paste0("/media/jaxi/differential_expression/unstranded/diff_expf/accepted_hits",i,".txt"), header=T)
      # rename the column containing actual integer counts
     names(count_files[[i+10]])[7] <- "u_counts"
     # use column 1 containing the exon identifier and column 7 containing the read counts
     count_files[[i+10]] <- count_files[[i+10]][,c(1,7)]
     # each exon and its counts should only correspond to one row, make sure this is the case
     count_files[[i+10]] <- count_files[[i+10]] %>% group_by(Geneid) %>% summarize(u_counts = sum(u_counts))
}	

# create a data frame with the Geneid's, this is done just so
# the cbind works
count_matrix <- data.frame(Geneid=count_files[[1]]$Geneid)

# change the name of the column appropriately
# and only keep the column that contains the read counts
for (i in 1:15) {
     
     if (i <= 10)
          new_name = paste0("accepted_0",i)
     else
          new_name = paste0("accepted_",i-10)
     
     count_matrix <- cbind(count_matrix, count_files[[i]][,2])
     names(count_matrix)[i+1] <- new_name
}

# set the row names to the identifier for each exon
row.names(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL

# write this table to a file for later use
write.table(count_matrix,"/media/jaxi/differential_expression/tests/exon_counts.txt")
