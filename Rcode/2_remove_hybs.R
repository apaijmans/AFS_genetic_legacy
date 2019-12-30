#------------------------------------------------------------
# Anneke Paijmans
# Created: Oct 2017
# Last edited: Mar 2019
#
# Seal microsat data - Remove hybrids identified using STRUCTURE
#
#------------------------------------------------------------

# Set working directory
library(here)


#~~ Load hybrid data (generated in: 1e_extract_hyb)

hybs <- read.table(here("Data", "Raw", "list_of_hybrids.txt"), header = TRUE)

# Add column with species info
hybs$sp <-ifelse(hybs$cluster1 > 0.90, "A.Tropicalis", "Hybrid")

library(tidyverse)
hybs %>% group_by(sp) %>% count(sp)
# sp               n
# <chr>        <int>
# 1 A.Tropicalis   135
# 2 Hybrid          20

hybs2 <- hybs

# Add population column 
hybs2$pop <- as.factor(ifelse(grepl("Bouvet", hybs2$sampleID), "Bouvetoya",
                              ifelse(grepl("Crozet", hybs2$sampleID), "Crozet",
                                     ifelse(grepl("Kergu", hybs2$sampleID), "Kerguelen",
                                            ifelse(grepl("Marion", hybs2$sampleID), "Marion",
                                                   ifelse(grepl("Mac", hybs2$sampleID), "Macquarie", "Tropicalis"))))))
hybs2 %>% group_by(sp) %>% count(pop)                                                           


# Load seal data
seal <- read.csv(here("Data", "A_gazella_microsatellite_dataset 20_03_19.csv"),
                 header = TRUE,
                 stringsAsFactors = FALSE,
                 sep = ",")

# Check data
# names(seal)
# head (seal)
# dim(seal)
# str(seal)
# length(seal)
# class(seal)

# Replace spaces in tissue name to compare hybrids list
seal$Sample <- gsub(" ", "_", seal$Sample, fixed = TRUE)


#~~ Remove all hybrids and A. tropicalis
# use list of hybrids to remove all hybrids and tropicalis animals
# because hybs sample ID only contains first 11 characters, use substring on seal dataset when comparing sample ID's

allgaz <- seal[!substring(seal$Sample, 1,11) %in% hybs$sampleID, ] #or subset(df, ID %in% keep)

# Add population column (-tropicalis because they are now removed)
allgaz$pop <- as.factor(ifelse(grepl("Macquarie", allgaz$Sample), "8_Macquarie",
                               ifelse(grepl("Heard", allgaz$Sample), "7_Heard",
                                      ifelse(grepl("Kergu", allgaz$Sample), "6_Kerguelen",
                                             ifelse(grepl("Crozet", allgaz$Sample), "5_Crozet",
                                                    ifelse(grepl("Marion", allgaz$Sample), "4_Marion",
                                                           ifelse(grepl("Bouvet", allgaz$Sample), "3_Bouvetoya",
                                                                  ifelse(grepl("Shet", allgaz$Sample), "1_SShetlands", "2_SGeorgia"))))))))

# Move pop col to front
allgaz <- allgaz[,c(1, 80, 2:79)]

# Order by sample and pop
allgaz <- allgaz[order(allgaz$Sample),] 
allgaz <- allgaz[order(allgaz$pop), ]


#~~ Save dataset 

write.table(allgaz, here("Data", "Processed", "Seal_no_hybs.txt"), sep="\t", row.names = FALSE, quote = FALSE) 

