#------------------------------------------------------------
# Anneke Paijmans
# Created: Jun 2018
# Last edited: Mar 2019
#
# Seal microsat data - Summary statistics
#                      1. per population 
#                      2. per cluster
#                      3. per population but with equal sample sizes for all pops
#
#------------------------------------------------------------

# Clear all previous info
#rm(list=ls()) # better restart R

# Set working directory and load libraries
library(here)

library("adegenet") 
library("hierfstat")
library("diveRsity")

# All populations
pops <- c("South Shetlands", 
          "South Georgia", 
          "Bouvetoya", 
          "Marion Island", 
          "Crozet Island", 
          "Kerguelen Island", 
          "Heard Island", 
          "Macquarie Island")


########################################################
#   Sum stats 1: all data per population (locality)  ###
########################################################


#~~ Load seal data

seal <- read.structure(here("Data", "Processed", "structure_in_gaz.stru"), n.ind = 2000, n.loc = 39, 
                       onerowperind = F, col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)


# number of genotypes: 2000 <- Use here number of individuals, not number of rows (per individual two rows)
# number of markers: 39
# which column labels for genotypes: 1
# which column population factor: 2
# Optional columns: enter
# Row marker names: 0 (absent)
# Genotypes coded by single row? n

# # seems number of alleles per locus for each pop are incorrect, more info on Github
# test <- read.structure("1 Data/ForPaper/structure_test_mac.stru", n.ind = 108, n.loc = 39, onerowperind = F,
#                        col.lab = 1, col.pop = 2, col.others = NULL,
#                        row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
#                        ask = F, quiet = FALSE)
# The immediate solution is to specify drop = TRUE in your seppop call

summary(seal)
min(as.data.frame(summary(seal)[3])) #2
max(as.data.frame(summary(seal)[3])) #26
mean(as.numeric(summary(seal)[[3]])) #11.44

write.csv(summary(seal)[3], here("Figs-Tables", "Basics", "Num_alleles_locus_gaz.csv"), quote=F)


#~~ Calculate allelic richness

ar <- allelic.richness(seal) #gives allelic richness per population per locus. Get mean and sd over loci

alri <- cbind(apply(as.data.frame(ar[2]), 2, mean), apply(as.data.frame(ar[2]), 2, sd), summary(seal)[[2]]) # mean, SD, and group size

alri <- cbind(alri[,3], alri[,1:2])

rownames(alri) <- pops
colnames(alri) <- c("Sample size", "Ar (mean)", "Ar (SD)")


#~~ Calculate number of private alleles per population per locus

priv_alleles <- poppr::private_alleles(seal, form = locus ~ ., count.alleles = F) # if count.alleles=T you get number of individuals with that allele within the pop

# Get percentages.
#sweep(priv_alleles, 2, nAll(seal)[colnames(priv_alleles)], FUN = "/")

# Get sum of private alleles over loci
pal <- as.data.frame(rowSums(priv_alleles))

rownames(pal) <- pops
colnames(pal) <- "PA"


#~~ Calculate N alleles/locus, Hobs and Hexp mean across 39 loci for each pop

# Split per population
sealpop <- seppop(seal, drop = TRUE)

npop = NULL
allN = NULL
H0 = NULL
He = NULL
blank <- rep(NA, 39)

for(i in 1:8) {
  
  smry <- summary(sealpop[[i]])
  
  popn <- smry[[1]]
  Nall <- smry[[3]]
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  
  npop <- rbind(npop, popn)
  allN <- rbind(allN, Nall)
  H0 <- rbind(H0, Hobs)
  He <- rbind(He, Hexp)
  
}


#~~ Calculate Fis per pop (mean and sd across loci) using diveRsity package

Fis <- basicStats(infile="C:/Uni/03-Demography/Data/Processed/genepop_gaz.txt", outfile=NULL)

# Bind everything together
y <- cbind(npop,
           apply(allN, 1, mean), 
           apply(allN, 1, sd), 
           apply(H0, 1, mean), 
           apply(H0, 1, sd), 
           apply(He, 1, mean), 
           apply(He, 1, sd),
           apply(Fis$fis[, 2:9], 2, mean),
           apply(Fis$fis[, 2:9], 2, sd)) # mean, SD, and group size

colnames(y) <- c("N per pop", "Alleles/locus (mean)", "Alleles/locus (SD)", "Hobs (mean)", "Hobs (SD)", "Hexp (mean)", "Hexp (SD)", "Fis (mean)", "Fis (SD)")
rownames(y) <- pops

sumstats <- cbind(alri, pal, y)


#~~ Save file

write.csv(sumstats, here("Figs-Tables", "Basics", "AllStats_gaz.csv"), quote=F)

##########################################
#   Sum stats 2: all data per cluster  ###
##########################################

# All clusters
clusters <- c("South Shetlands", 
              "South Georgia", 
              "Bouvetoya", 
              "Marion Island-Crozet Island", 
              "Kerguelen Island-Heard Island-Macquarie Island")


#~~ Load seal data

seal <- read.structure(here("Data", "Processed", "structure_in_gaz_cl.stru"), n.ind = 2000, n.loc = 39, 
                       onerowperind = F, col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)


summary(seal)
min(as.data.frame(summary(seal)[3])) #2
max(as.data.frame(summary(seal)[3])) #26
mean(as.numeric(summary(seal)[[3]])) #11.44

write.csv(summary(seal)[3], here("Figs-Tables", "Basics", "Num_alleles_locus_gaz_cl.csv"), quote=F)


#~~ Calculate allelic richness

ar <- allelic.richness(seal, min.n=22) #gives allelic richness per population per locus. Get mean and sd over loci
# specify min.n thesame as per pop analysis to make it comparable with results per pop

alri <- cbind(apply(as.data.frame(ar[2]), 2, mean), apply(as.data.frame(ar[2]), 2, sd), summary(seal)[[2]]) # mean, SD, and group size

alri <- cbind(alri[,3], alri[,1:2])

rownames(alri) <- clusters
colnames(alri) <- c("Sample size", "Ar (mean)", "Ar (SD)")


#~~ Calculate number of private alleles per population per locus

priv_alleles <- poppr::private_alleles(seal, form = locus ~ ., count.alleles = F) # if count.alleles=T you get number of individuals with that allele within the pop

# Get sum of private alleles over loci
pal <- as.data.frame(rowSums(priv_alleles))

rownames(pal) <- clusters
colnames(pal) <- "PA"


#~~ Calculate N alleles/locus, Hobs and Hexp mean across 39 loci for each pop

# Split per population
sealpop <- seppop(seal, drop = TRUE)

npop = NULL
allN = NULL
H0 = NULL
He = NULL
blank <- rep(NA, 39)

for(i in 1:5) {
  
  smry <- summary(sealpop[[i]])
  
  popn <- smry[[1]]
  Nall <- smry[[3]]
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  
  npop <- rbind(npop, popn)
  allN <- rbind(allN, Nall)
  H0 <- rbind(H0, Hobs)
  He <- rbind(He, Hexp)
  
}


#~~ Calculate Fis per pop (mean and sd across loci) using diveRsity package

Fis <- basicStats(infile="C:/Uni/03-Demography/Data/Processed/genepop_gaz_cl.txt", outfile=NULL)

# Bind everything together
y <- cbind(npop,
           apply(allN, 1, mean), 
           apply(allN, 1, sd), 
           apply(H0, 1, mean), 
           apply(H0, 1, sd), 
           apply(He, 1, mean), 
           apply(He, 1, sd),
           apply(Fis$fis[, 2:6], 2, mean),
           apply(Fis$fis[, 2:6], 2, sd)) # mean, SD, and group size

colnames(y) <- c("N per cluster", "Alleles/locus (mean)", "Alleles/locus (SD)", "Hobs (mean)", "Hobs (SD)", "Hexp (mean)", "Hexp (SD)", "Fis (mean)", "Fis (SD)")
rownames(y) <- clusters

sumstats <- cbind(alri, pal, y)


#~~ Save file

write.csv(sumstats, here("Figs-Tables", "Basics", "AllStats_gaz_cl.csv"), quote=F)


##################################################################
#   Sum stats 3: equal sample sizes per population (locality)  ###
##################################################################


# Clear all previous info, except pops info
rm(list=setdiff(ls(), "pops"))

#~~ Read seal data
seal <- read.structure(here("Data", "Processed", "structure_in_equal_size.stru"), n.ind = 144, n.loc = 39, 
                       onerowperind = F, col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

# number of genotypes: 152 <- Use here number of individuals, not number of rows (per individual two rows)

summary(seal)


#~~ Calculate allelic richness

ar <- allelic.richness(seal) #gives allelic richness per population per locus. Get mean and sd over loci

alri <- cbind(apply(as.data.frame(ar[2]), 2, mean), apply(as.data.frame(ar[2]), 2, sd), summary(seal)[[2]]) # mean, SD, and group size

alri <- cbind(alri[,3], alri[,1:2])

rownames(alri) <- pops
colnames(alri) <- c("Sample size", "Ar (mean)", "Ar (SD)")


#~~ Calculate N alleles/locus, Hobs and Hexp mean across 39 loci for each pop

# Split per population
sealpop <- seppop(seal, drop = TRUE)

allN = NULL
H0 = NULL
He = NULL
blank <- rep(NA, 39)

for(i in 1:8) {
  
  smry <- summary(sealpop[[i]])
  
  Nall <- smry[[3]]
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  
  allN <- rbind(allN, Nall)
  H0 <- rbind(H0, Hobs)
  He <- rbind(He, Hexp)
}


#~~ Bind everything together

y <- cbind(apply(allN, 1, mean), 
           apply(allN, 1, sd), 
           apply(H0, 1, mean), 
           apply(H0, 1, sd), 
           apply(He, 1, mean), 
           apply(He, 1, sd)) # mean, SD, and group size

colnames(y) <- c("Alleles/locus (mean)", "Alleles/locus (SD)", "Hobs (mean)", "Hobs (SD)", "Hexp (mean)", "Hexp (SD)")
rownames(y) <- pops

sumstats <- cbind(alri, y)


#~~ Sace file

write.csv(sumstats, here("Figs-Tables", "Basics", "AllStats_equal_size.csv"), quote=F)
