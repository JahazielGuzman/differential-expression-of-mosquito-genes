library(dplyr)
count_files <- list()

# 10 stranded files, 5 unstranded files
strand_num = 1:10
unstrand_num = 1:5

for (i in unstrand_num) {	

     count_files[[i]] <- read.table(paste0("/media/jaxi/differential_expression/unstranded/diff_expf/accepted_hits",i,".txt"), header=T)
     names(count_files[[i]])[7] <- "u_counts"
     count_files[[i]] <- count_files[[i]][,c(1,7)]
     count_files[[i]] <- count_files[[i]] %>% group_by(Geneid) %>% summarize(u_counts = sum(u_counts))
}	

for (i in strand_num) {	
     
     count_files[[i + 5]] <- read.table(paste0("/media/jaxi/differential_expression/stranded/diff_expf/accepted_hits0",i,".txt"), header=T)
     names(count_files[[i + 5]])[7] <- "s_counts"
     count_files[[i + 5]] <- count_files[[i + 5]][,c(1,7)]
     count_files[[i + 5]] <- count_files[[i + 5]] %>% group_by(Geneid) %>% summarize(s_counts = sum(s_counts))
}



count_matrix <- data.frame(Geneid=count_files[[1]]$Geneid)

for (i in 1:15) {
     
     if (i > 5)
          new_name = paste0("accepted_0",i-5)
     else
          new_name = paste0("accepted_",i)
     
     count_matrix <- cbind(count_matrix, count_files[[i]][,2])
     names(count_matrix)[i+1] <- new_name
}

row.names(count_matrix) <- count_matrix$Geneid
count_matrix$Geneid <- NULL

write.table(count_matrix,"/media/jaxi/differential_expression/tests/exon_counts.txt")
