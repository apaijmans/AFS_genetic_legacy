#------------------------------------------------------------
# Anneke Paijmans
# Created: Jun 2018
# Last edited: Dec 2019
#
# Seal microsat data - Prepare data without hybrids for
#                      STRUCTURE
#                      Make subsets: 
#                      1. all data per population (locality) 
#                      2. subset with equal sample sizes per population (locality)
#                      3. all data grouped per genetic cluster
#
#------------------------------------------------------------

# Clear all previous info
#rm(list=ls()) # better restart R

# Set working directory
library(here)


#~~ Load fur seal microsat data without hybrids

seal <- read.table(here("Data", "Processed", "Seal_no_hybs.txt"),
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   sep = "\t")

# Replace x with -9 for STRUCTURE
seal[seal == "x" ] <- -9


#~~ Split data and bind later, so that we get per individual two rows

seala <- seal[c(T,F)]
sealb <- seal[c(F,T)]

# Add pop column to seala
seala$pop <- seal$pop

# Reorder columns
seala <- seala[,c(1, 41, 2:40)]

# Add sample column to seala
sealb$Sample <- seal$Sample

# Reorder columns
sealb <- sealb[,c(41, 1:40)]


#~~ Rename the col names so that they bind under each other

names(sealb) <- sub(".b", "", names(seala), fixed=TRUE)#fixed true otherwise it overrides any b

# check whether col names are the same in both files
# names(seala) == names(sealb)


#~~ Bind data together

seal.tot <- rbind(seala, sealb)

#Order by sample and pop
seal.tot <- seal.tot[order(seal.tot$Sample),] 
seal.tot <- seal.tot[order(seal.tot$pop), ]


###########################################################
#   Make dataset 1: all data per population (locality)  ###
###########################################################

# Prep for running STRUCTURE on server
library(data.table)
library(dplyr)

# Check number of individuals/population
seal.tot %>%
  group_by(pop) %>%
  summarize(n = length(Sample)/2)

# Convert population names to numbers
struc <- seal.tot %>%
  mutate(pop = ifelse(grepl("SShetlands", pop), 1,
                      ifelse(grepl("SGeorgia", pop), 2,
                             ifelse(grepl("Bouvetoya", pop), 3,
                                    ifelse(grepl("Marion", pop), 4, 
                                           ifelse(grepl("Crozet", pop), 5,
                                                  ifelse(grepl("Kerguelen", pop), 6,
                                                         ifelse(grepl("Heard", pop), 7, 
                                                                ifelse(grepl("Macquarie", pop), 8, 9)))))))))


# Save as STRUCTURE file
fwrite(struc, here("Data", "Processed", "structure_in_gaz.stru"), col.names = F,
       quote = F, row.names = F, sep = " ")


############################################################################
#   Make dataset 2: equal samples sizes for each populations (locality)  ###
############################################################################


# Check minimum pop size
seal.tot %>%
  group_by(pop) %>%
  summarize(n = length(Sample)/2)

# pop              n
# 1 1_SShetlands   197
# 2 2_SGeorgia    1042
# 3 3_Bouvetoya    396
# 4 4_Marion       166
# 5 5_Crozet        18 <- smallest pop size
# 6 6_Kerguelen     51
# 7 7_Heard         22
# 8 8_Macquarie    108

# Initiate random number generator engine. This important for users to reproduce the analysis.
set.seed(42) 

# Because each individual has 2 rows, it means there are duplicate indivuals 
# therefore make data with unique sample IDs first
ind <- unique(seal.tot[,c("Sample", "pop")])

rand_samples <- tapply(ind$Sample, ind$pop, sample, replace=F, size=18) # this creates a list
#duplicated(rand_samples[[5]]) # check to see if there are any duplicates (shouldnt be)
rand_samples <- c(do.call("cbind",rand_samples)) # converts list into vector, so I can select rows matching this vector

all_rand <- seal.tot[seal.tot$Sample %in% rand_samples, ]

all_rand %>%
  group_by(pop) %>%
  summarize(n = length(Sample)/2)

# Convert population names to numbers
struc <- all_rand %>%
  mutate(pop = ifelse(grepl("SShetlands", pop), 1,
                      ifelse(grepl("SGeorgia", pop), 2,
                             ifelse(grepl("Bouvetoya", pop), 3,
                                    ifelse(grepl("Marion", pop), 4, 
                                           ifelse(grepl("Crozet", pop), 5,
                                                  ifelse(grepl("Kerguelen", pop), 6,
                                                         ifelse(grepl("Heard", pop), 7, 
                                                                ifelse(grepl("Macquarie", pop), 8, 9)))))))))

# Save as STRUCTURE file
fwrite(struc, here("Data", "Processed", "structure_in_equal_size.stru"), col.names = F,
       quote = F, row.names = F, sep = " ")


#########################################################################################################
#   Make dataset 3: all data per genetic cluster (as identified by STRUCTURE in: 3c_parse_structure)  ###
#########################################################################################################

# Prep for running STRUCTURE on server
library(data.table)
library(dplyr)

# Check number of individuals/population
seal.tot %>%
  group_by(pop) %>%
  summarize(n = length(Sample)/2)

# Convert population names to clusters
struc <- seal.tot %>%
  mutate(pop = ifelse(grepl("SShetlands", pop), 1,
                      ifelse(grepl("SGeorgia", pop), 2,
                             ifelse(grepl("Bouvetoya", pop), 3,
                                    ifelse(grepl("Marion", pop), 4, 
                                           ifelse(grepl("Crozet", pop), 4,
                                                  ifelse(grepl("Kerguelen", pop), 5,
                                                         ifelse(grepl("Heard", pop), 5, 
                                                                ifelse(grepl("Macquarie", pop), 5, 9)))))))))

struc %>%
  group_by(pop) %>%
  summarize(n = length(Sample)/2)

# Save as STRUCTURE file
fwrite(struc, here("Data", "Processed", "structure_in_gaz_cl.stru"), col.names = F,
       quote = F, row.names = F, sep = " ")

