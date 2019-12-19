#------------------------------------------------------------
# Anneke Paijmans
# Created: June 2017
# Last edited: Mar 2019
#
# Seal microsat data - Prepare input file for STRUCTURE
#
#------------------------------------------------------------

# Clear all previous info 
#rm(list=ls()) # better restart R

# Set working directory
library(here)

#~~ Read fur seal microsat data

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


#~~ Check for repeats

library(dplyr)
seal %>% filter(Sample %in% unique(.[["Sample"]][duplicated(.[["Sample"]])]))


#~~ Replace x with -9 for STRUCTURE (but not in population or ID colums)

seal[seal == "x" ] <- -9


#~~ Make for each individual two rows containing genotypes

# Split data and bind later, so that we get per individual two rows
seala <- seal[c(F,T)]
seala <- cbind(seal[c(1)], seala[c(1:39)])
      
sealb <- seal[c(T,F)]
#sealb <- cbind(seal[c(1)], sealb)

# Rename the col names so that they bind under each other
names(sealb) <- sub(".b", "", names(seala), fixed=TRUE)#fixed true otherwise it overrides any b

# check whether col names are the same in both files
# names(seala) == names(sealb)

# Bind data together
seal.tot <- rbind(seala, sealb)


#~~ Make column with population info

seal.tot$pop <- as.factor(ifelse(grepl(" Macquarie", seal.tot$Sample), "8_Macquarie",
                                 ifelse(grepl("Heard", seal.tot$Sample), "7_Heard",
                                        ifelse(grepl("Kergu", seal.tot$Sample), "6_Kerguelen",
                                               ifelse(grepl("Crozet", seal.tot$Sample), "5_Crozet",
                                                      ifelse(grepl("Marion", seal.tot$Sample), "4_Marion",
                                                             ifelse(grepl("Bouvet", seal.tot$Sample), "3_Bouvetoya",
                                                                    ifelse(grepl("Tropicalis", seal.tot$Sample), "9_Trop_Mac",
                                                                           ifelse(grepl("Shet", seal.tot$Sample), "1_SShetlands", "2_SGeorgia")))))))))

# Move pop col to front
seal.tot <- seal.tot[,c(1, 41, 2:40)]


#~~ Tidying up

# Replace spaces in tissue name and population otherwise STRUCTURE will think they are different columns
seal.tot$Sample <- gsub(" ", "_", seal.tot$Sample, fixed = TRUE)
#seal.tot$pop <- gsub(" ", "_", seal.tot$pop, fixed = TRUE)

# Order the file by sample ID, so that the two entries for each individual are below each other
# Then order by location so that locations are grouped together
seal.tot <- seal.tot[order(seal.tot$Sample),] 
seal.tot <- seal.tot[order(seal.tot$pop), ]

#~~ Save as STRUCTURE file for running on server

library(data.table)

struc <- seal.tot %>% #fread("1 Data/ForPaper/SealStructureInput2019.txt") %>%
  mutate(pop = ifelse(grepl("SShetlands", pop), 1,
                             ifelse(grepl("SGeorgia", pop), 2,
                                    ifelse(grepl("Bouvetoya", pop), 3,
                                           ifelse(grepl("Marion", pop), 4, 
                                                  ifelse(grepl("Crozet", pop), 5,
                                                         ifelse(grepl("Kerguelen", pop), 6,
                                                                ifelse(grepl("Heard", pop), 7, 
                                                                       ifelse(grepl("Macquarie", pop), 8, 9)))))))))


fwrite(struc, here("Data", "Processed", "structure_in_trop.stru"), col.names = F,
       quote = F, row.names = F, sep = " ")

