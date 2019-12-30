#------------------------------------------------------------
# Anneke Paijmans
# last edited: Mar 2019
#
# Seal microsat data - Making subsets from data for BOTTLENECK software
#                      Output: stru files, convert with PGD spider to genepop files
#                      1. per population 
#                      2. per cluster
#
#------------------------------------------------------------


# Set working directory and load libraries
library(here)

library(tidyverse)


##################################################
#   1. Make files per populations (locality)   ###
##################################################


#~~ Load seal data

seal <- read.table(here("Data", "Processed", "structure_in_gaz.stru"))

seal <- plyr::rename(seal, c("V1"="Sample", "V2"="pop"))
seal$Sample <- as.character(seal$Sample)

#~~ Due to a mistake in the data we needed to re-run the BOTTLENECK analysis for SG (pop2) and Mac (pop8)
#~~ So we repeat the code for these pops only

# First make subset of data with only these pops
seal2 <- rbind(seal[seal$pop==2,], seal[seal$pop==8,])

# Because each individual has 2 rows, it means there are duplicate indivuals
# therefore make data with unique sample IDs first
ind <- unique(seal2[,c("Sample", "pop")])


#~~ Make 1000 subset of 10 individuals for each pop each by looping

for (i in 1:1000) {
  set.seed(i)
  
  rand_samples <- tapply(ind$Sample, ind$pop, sample, replace=F, size=10) # this creates a list
  
  #duplicated(rand_samples[[5]]) # check to see if there are any duplicates (shouldnt be)
  
  rand_samples <- c(do.call("cbind",rand_samples)) # converts list into vector, so I can select rows matching this vector
  
  all_rand <- seal2[seal2$Sample %in% rand_samples, ]
  
  readr::write_delim(x =  all_rand, path = paste0("Data/Processed/bottleneck_files/raw/SG_MAC/SG_MACbottleneck_in", i, ".stru"), col_names = F)
  
}


#~~ convert resulting stru files with PGD spider to genepop files
#~~ Use runSpider.bat/runSpiderSG_MAC.bat file for fast converting

###########################################
#   2. Make files per genetic cluster   ###
###########################################


#~~ Load seal data

seal <- read.table(here("Data", "Processed", "structure_in_gaz.stru"))

seal <- plyr::rename(seal, c("V1"="Sample", "V2"="pop"))
seal$Sample <- as.character(seal$Sample)

# Relabel pops to clusters (ie 1=1, 2=2, 3=3, 4=4&5, 5= 6&7&8)
seal$cl <- gsub(5, 4, seal$pop, fixed = TRUE)
seal$cl <- gsub(6, 5, seal$cl, fixed = TRUE)
seal$cl <- gsub(7, 5, seal$cl, fixed = TRUE)
seal$cl <- gsub(8, 5, seal$cl, fixed = TRUE)

seal$pop <- seal$cl
seal$cl <- NULL


#~~ Due to a mistake in the data we needed to re-run the BOTTLENECK analysis for SG (pop2) and Mac (pop8)
#~~ So we repeat the code for these clusters only

# First make subset of data with only these clusters
seal2 <- rbind(seal[seal$pop==2,], seal[seal$pop==5,])

# Because each individual has 2 rows, it means there are duplicate indivuals
# therefore make data with unique sample IDs first
ind <- unique(seal2[,c("Sample", "pop")])


#~~ Make 1000 subset of 10 individuals for each pop each by looping

for (i in 1:1000) {
  set.seed(i)
  
  rand_samples <- tapply(ind$Sample, ind$pop, sample, replace=F, size=10) # this creates a list
  
  #duplicated(rand_samples[[5]]) # check to see if there are any duplicates (shouldnt be)
  
  rand_samples <- c(do.call("cbind",rand_samples)) # converts list into vector, so I can select rows matching this vector
  
  all_rand <- seal2[seal2$Sample %in% rand_samples, ]
  
  readr::write_delim(x =  all_rand, path = paste0("Data/Processed/bottleneck_files/raw/cluster/SG_MAC/SG_MAC_CLbottleneck_in", i, ".stru"), col_names = F)
  
}


#~~ convert resulting stru files with PGD spider to genepop files
#~~ Use runSpiderCl.bat/runSpiderSG_MAC_CL.bat file for fast converting
