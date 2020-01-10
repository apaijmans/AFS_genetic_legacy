# Clear all previous info
#rm(list=ls())

library(data.table)
library(pophelper)

# Load one of the structure files and save as slist
slist <- readQ("results/run_files/results_job_T1_1_f", filetype = "structure", indlabfromfile=T)

# Extract sample ID from slist (which got them from the structure file using indlabfromfile=T)
indnames <- attr(slist[[1]],"row.names")

# Load K=2 merged file, containing assignment of each individual to cluster 1 or 2
K2data <- read.table("results/run_files/pop-merged/pop_K2-combined-merged.txt")

# Add name list to the file and change column names
K2data$sampleID <- indnames
K2data$cluster1 <- K2data$V2
K2data$cluster2 <- K2data$V3

# Remove redundant columns
K2data[, 1:4] <- NULL

# Extract individuals with more than 0.1 assignment to cluster 2 (ie tropicalis) so more than 10% assigned to A tropicalis
Hyb <- K2data[which(K2data$cluster1 > 0.1),] # 154 individuals

# Save the datasets
system("mkdir results/run_files/hybs")

write.table(K2data, "results/run_files/hybs/pop_K2-combined-mergedINCLsampleID.txt", sep="\t", row.names = FALSE, quote = FALSE) 

write.table(Hyb, "results/run_files/hybs/list_of_hybrids.txt", sep="\t", row.names = FALSE, quote = FALSE) 

# We decide to use the 10% as a cut off, as we think it is better to be conservative. See also paper by Vaha 2006
