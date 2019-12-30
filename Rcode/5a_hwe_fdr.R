#------------------------------------------------------------
# Anneke Paijmans
# Created: June 2017
# Edited: Mar 2019
#
# Seal microsat data - Data summary and calculate HWE
#                      1. all data per population (locality)
#                      2. equal sample size per population
#
#------------------------------------------------------------


# Clear all previous info
#rm(list=ls()) # better restart R

# Set working directory and load libraries
library(here)

library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet") 


##################################################
#   HWE 1: all data per population (locality)  ###
##################################################


#~~ Load seal data

seal <- read.structure(here("Data", "Processed", "structure_in_gaz.stru"), n.ind = 2000, n.loc = 39, 
                       onerowperind = F, col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

# number of genotypes: 2000
# number of markers: 39
# which column labels for genotypes: 1
# which column population factor: 2
# Optional columns: enter
# Row marker names: 0 (absent)
# Genotypes coded by single row? n

summary(seal)
# Group sizes: 197 1042 396 166 18 51 22 108

# Split per population
sealpop <- seppop(seal)

#test <- summary(sealpop[[1]])


#~~ Get observed Heterozygosity (Hobs) and expected heterozygosity (Hexp) for 
#   each population seperatly using a loop

# Set up values/vectors for loop
reps = 10000 # number of replicates for the Monte Carlo procedure 

allHWE = NULL
blank <- rep(NA, 39)

for(i in 1:8) {
  hwt <- hw.test(sealpop[[i]], B=reps)
  smry <- summary(sealpop[[i]])
  
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  pexact <- hwt[,4] #hw.test does chi2 test and exact test. We use p-values of exact test
  
  allHWE <- cbind(allHWE, Hobs, Hexp, pexact, blank) 
  
}


#~~ Check p values. Not all uniform, so FDR might not be optimal

hist(allHWE[, 3])
hist(allHWE[, 7])
hist(allHWE[, 11])
hist(allHWE[, 15])
hist(allHWE[, 19])
hist(allHWE[, 23])
hist(allHWE[, 27])
hist(allHWE[, 31])

# Run multiple comparisons tests: FDR and Bonferroni
# library("qvalue")
# q1 <- qvalue(allHWE[, 3])
# q2 <- qvalue(allHWE[, 7]) #doesnt work as no p value above 0.95
# q3 <- qvalue(allHWE[, 11])
# q4 <- qvalue(allHWE[, 15])
# q5 <- qvalue(allHWE[, 19])
# q6 <- qvalue(allHWE[, 23])
# q7 <- qvalue(allHWE[, 26])#doesnt work as no p value above 0.95
# q8 <- qvalue(allHWE[, 31])

# Apply Bonferroni correction
b1 <- p.adjust(allHWE[, 3], method = "bonferroni")
b2 <- p.adjust(allHWE[, 7], method = "bonferroni")
b3 <- p.adjust(allHWE[, 11], method = "bonferroni")
b4 <- p.adjust(allHWE[, 15], method = "bonferroni")
b5 <- p.adjust(allHWE[, 19], method = "bonferroni")
b6 <- p.adjust(allHWE[, 23], method = "bonferroni")
b7 <- p.adjust(allHWE[, 27], method = "bonferroni")
b8 <- p.adjust(allHWE[, 31], method = "bonferroni")

# Compare FDR with Bonferroni
# qVSbon <- NULL

# qVSbon <- cbind(qVSbon,q1$qvalues, b1, blank, b2, q3$qvalues, b3, q4$qvalues, b4, q5$qvalues, b5, q6$qvalues, b6, blank, b7, q8$qvalues, b8)

# Put Bonferroni into overal file
allHWE[, 4] <- b1
allHWE[, 8] <- b2
allHWE[, 12] <- b3
allHWE[, 16] <- b4
allHWE[, 20] <- b5
allHWE[, 24] <- b6
allHWE[, 28] <- b7
allHWE[, 32] <- b8


#~~ Save file

write.csv(allHWE, here("Figs-Tables", "HWE", "AllLocsHWE_all.csv"), quote=F, row.names = F)

############################################################
#   HWE 2: equal sample sizes per population (locality)  ###
############################################################

# Clear all previous info
rm(list=ls())


#~~ Load seal data

seal <- read.structure(here("Data", "Processed", "structure_in_equal_size.stru"), n.ind = 144, n.loc = 39, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

# number of genotypes: 144 

summary(seal)
# Group sizes: 18 18 18 18 18 18 18 18

# Split per population
sealpop <- seppop(seal)

#test <- summary(sealpop[[1]])


#~~ Get observed Heterozygosity (Hobs) and expected heterozygosity (Hexp) for 
#   each population seperatly using a loop

# Set up values/vectors for loop
reps = 10000 # number of replicates for the Monte Carlo procedure 

allHWE = NULL
blank <- rep(NA, 39)

for(i in 1:8) {
  hwt <- hw.test(sealpop[[i]], B=reps)
  smry <- summary(sealpop[[i]])
  
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  pexact <- hwt[,4] #hw.test does chi2 test and exact test. We use pvalues of exact test
  
  allHWE <- cbind(allHWE, Hobs, Hexp, pexact, blank) 
  
}


#~~ Check p-values. Not all uniform, so FDR might not be optimal

hist(allHWE[, 3])
hist(allHWE[, 7])
hist(allHWE[, 11])
hist(allHWE[, 15])
hist(allHWE[, 19])
hist(allHWE[, 23])
hist(allHWE[, 27])
hist(allHWE[, 31])

# Run multiple comparisons test
# library("qvalue")
# q1 <- qvalue(allHWE[, 3])
# q2 <- qvalue(allHWE[, 7])
# q3 <- qvalue(allHWE[, 11])
# q4 <- qvalue(allHWE[, 15])
# q5 <- qvalue(allHWE[, 19])
# q6 <- qvalue(allHWE[, 23])
# q7 <- qvalue(allHWE[, 26])
# q8 <- qvalue(allHWE[, 31])

# Apply Bonferroni correction
b1 <- p.adjust(allHWE[, 3], method = "bonferroni")
b2 <- p.adjust(allHWE[, 7], method = "bonferroni")
b3 <- p.adjust(allHWE[, 11], method = "bonferroni")
b4 <- p.adjust(allHWE[, 15], method = "bonferroni")
b5 <- p.adjust(allHWE[, 19], method = "bonferroni")
b6 <- p.adjust(allHWE[, 23], method = "bonferroni")
b7 <- p.adjust(allHWE[, 27], method = "bonferroni")
b8 <- p.adjust(allHWE[, 31], method = "bonferroni")

# Compare FDR with Bonferroni
#qVSbon <- NULL

# qVSbon <- cbind(qVSbon,q1$qvalues, b1, blank, b2, q3$qvalues, b3, q4$qvalues, b4, q5$qvalues, b5, q6$qvalues, b6, blank, b7, q8$qvalues, b8)

# Put Bonferroni into overal file
allHWE[, 4] <- b1
allHWE[, 8] <- b2
allHWE[, 12] <- b3
allHWE[, 16] <- b4
allHWE[, 20] <- b5
allHWE[, 24] <- b6
allHWE[, 28] <- b7
allHWE[, 32] <- b8


#~~ Save file

write.csv(allHWE, here("Figs-Tables", "HWE", "AllLocsHWE_equal_size.csv"), quote=F, row.names = F)
