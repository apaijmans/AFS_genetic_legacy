#------------------------------------------------------------
# Anneke Paijmans
# Created: June 2017
# Edited: Mar 2019
#
# Seal microsat data - Calculate Nei's Fst (also Gst. Fst is only for 2 alleles, 
#                       but as we have msat data, we calculate Gst)
#                      1. per population 
#                      2. per cluster
#                      3. per population but with equal sample sizes for all pops
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
library("hierfstat")


##################################################
#   Fst 1: all data per population (locality)  ###
##################################################


# All populations
pops <- c("South Shetlands", 
          "South Georgia", 
          "Bouvetoya", 
          "Marion Island", 
          "Crozet Island", 
          "Kerguelen Island", 
          "Heard Island", 
          "Macquarie Island")


#~~ Load seal data

seal <- read.structure(here("Data", "Processed", "structure_in_gaz.stru"), n.ind = 2000, n.loc = 39, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
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


#~~ Convert genind file to hierfstat file for bootstrapping

Fboot <- genind2hierfstat(seal)
Fboot$pop <- as.numeric(Fboot$pop)


#~~ Calculate pairwise Fst

Fstats <- basic.stats(Fboot) # gives Neis Fst 1987 (among other data)
FstNei <- pairwise.neifst(Fboot) # FST are calculated using Nei (1987) equations for Fst


#~~ Calculate bootstrapped confidence intervals

ci_fst <- boot.ppfst(dat=Fboot, nboot=1000, diploid=TRUE)

# boot.ppfis and boot.ppfst provide bootstrap confidence intervals (bootstrapping over loci) 
# for population specific Fis and pairwise Fst respectively.


#~~ Make result matrix: 
#   Fst value below diagonal, 
#   CI above consisting of both lower limit and upper limit seperated with -

ll <- ci_fst[[2]][upper.tri(ci_fst[[2]])] # lower limit: ci_fst[[2]], saved in upper part of matrix
ul <- ci_fst[[3]][upper.tri(ci_fst[[3]])] # upper limit: ci_fst[[2]], also saved in upper part of matrix

ll <- round(ll, 3)
ul <- round(ul, 3)

ci <- paste(ll,"-",ul)

m <- round(FstNei, 3)
m[upper.tri(m)] <- ci

diag(m) <- ""

rownames(m) <- pops
colnames(m) <- pops


#~~ Save file

write.csv(m, here("Figs-Tables", "Fst", "FstMatrix_gaz.csv"), quote = FALSE)


####################################
#   Fst 2: all data per cluster  ###
####################################

# All clusters
clusters <- c("South Shetlands", "South Georgia", "Bouvetoya", "Marion Island-Crozet Island", "Kerguelen Island-Heard Island-Macquarie Island")


#~~ Read seal data

seal <- read.structure(here("Data", "Processed", "structure_in_gaz_cl.stru"), n.ind = 2000, n.loc = 39, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)


#~~ Convert genind file to hierfstat file for bootstrapping

Fboot <- genind2hierfstat(seal)
Fboot$pop <- as.numeric(Fboot$pop)


#~~ Calculate pairwise Fst

Fstats <- basic.stats(Fboot) # gives Neis Fst 87 (among other data)
FstNei <- pairwise.neifst(Fboot) # FST are calculated using Nei (87) equations for FST'. Seems to fit better with the bootstrapping function, as the CI were outside of range of FST as calculated before...?


#~~ Calculate bootstrapped confidence intervals

ci_fst <- boot.ppfst(dat=Fboot, nboot=1000, diploid=TRUE)


#~~ Make result matrix: 
#   Fst value below diagonal, 
#   CI above consisting of both lower limit and upper limit seperated with -

ll <- ci_fst[[2]][upper.tri(ci_fst[[2]])] # lower limit: ci_fst[[2]], saved in upper part of matrix
ul <- ci_fst[[3]][upper.tri(ci_fst[[3]])] # upper limit: ci_fst[[2]], also saved in upper part of matrix

ll <- round(ll, 3)
ul <- round(ul, 3)

ci <- paste(ll,"-",ul)

m <- round(FstNei, 3)
m[upper.tri(m)] <- ci

diag(m) <- ""

rownames(m) <- clusters
colnames(m) <- clusters

# lt <- m[lower.tri(m)] # going first down, then right, down, etc.
# ut <- m[upper.tri(m)] # same, but for upper half


#~~ Save file

write.csv(m, here("Figs-Tables", "Fst", "FstMatrix_gaz_cl.csv"), quote = FALSE)


############################################################
#   Fst 3: equal sample sizes per population (locality)  ###
############################################################


# Clear all previous info, except pops info
rm(list=setdiff(ls(), "pops"))


#~~ Load seal data

seal <- read.structure(here("Data", "Processed", "structure_in_equal_size.stru"), n.ind = 144, n.loc = 39, 
                       onerowperind = F, col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

# number of genotypes: 152

summary(seal)


#~~ Convert genind file to hierfstat file for bootstrapping

Fboot <- genind2hierfstat(seal)
Fboot$pop <- as.numeric(Fboot$pop)


#~~ Calculate pairwise Fst

Fstats <- basic.stats(Fboot) # gives Neis Fst 87 (among other data)
FstNei <- pairwise.neifst(Fboot) # FST are calculated using Nei (87) equations for FST'. Seems to fit better with the bootstrapping function, as the CI were outside of range of FST as calculated before...?


#~~ Calculate bootstrapped confidence intervals

ci_fst <- boot.ppfst(dat=Fboot,nboot=1000,diploid=TRUE)


#~~ Make result matrix: 
#   Fst value below diagonal, 
#   CI above consisting of both lower limit and upper limit seperated with -

ll <- ci_fst[[2]][upper.tri(ci_fst[[2]])] # lower limit: ci_fst[[2]], saved in upper part of matrix
ul <- ci_fst[[3]][upper.tri(ci_fst[[3]])] # upper limit: ci_fst[[2]], also saved in upper part of matrix

ll <- round(ll, 3)
ul <- round(ul, 3)

ci <- paste(ll,"-",ul)

m <- round(FstAll, 3)
m[upper.tri(m)] <- ci

diag(m) <- ""

rownames(m) <- pops
colnames(m) <- pops

# lt <- m[lower.tri(m)] # going first down, then right, down, etc.
# ut <- m[upper.tri(m)] # same, but for upper half


#~~ Save file

write.csv(m, here("Figs-Tables", "Fst", "FstMatrix_equal_size.csv"), quote = FALSE)
