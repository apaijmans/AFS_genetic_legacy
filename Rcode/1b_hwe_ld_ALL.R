#------------------------------------------------------------
# Anneke Paijmans
# Created: June 2017
# Last edited: Mar 2019
#
# Seal microsat data - Test HWE and LD in Antarctic 
#                      and Subantarctic fur seal samples
#
#------------------------------------------------------------

# Clear all previous info 
#rm(list=ls()) # better restart R

# Set working directory
library(here)

# Load data
library("adegenet") 

seal <- read.structure(here("Data", "Processed", "structure_in_trop.stru"), n.ind = 2155, n.loc = 39, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

# Check for NAs
poppr::info_table(seal, plot = TRUE)

# Split per population
sealpop <- seppop(seal)

#test <- summary(sealpop[[1]])


###########################################################
#   Test for HWE in Antarctic and subantarctic samples  ###
###########################################################

#~~ Get observed Heterozygosity (Hobs) and expected heterozygosity (Hexp) for 
#   each population seperatly using a loop

# Set up values/vectors for loop
reps = 10000 # number of replicates for the Monte Carlo procedure 

allHWE = NULL
blank <- rep(NA, 39)

# Run loop
for(i in 1:9) {
  hwt <- pegas::hw.test(sealpop[[i]], B=reps)
  smry <- summary(sealpop[[i]])
  
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  pexact <- hwt[,4] #hw.test does chi2 test and exact test. We use p-values of exact test which are given in 4th col
  
  allHWE <- cbind(allHWE, Hobs, Hexp, pexact, blank) 
  
}

#~~ Check p values for all 9 populations (8 Antarctic fur seal pops, 1 Subantarctic fur seal pop)

hist(allHWE[, 3])
hist(allHWE[, 7])
hist(allHWE[, 11])
hist(allHWE[, 15])
hist(allHWE[, 19])
hist(allHWE[, 23])
hist(allHWE[, 27])
hist(allHWE[, 31])
hist(allHWE[, 35])
#They are not all uniform, so FDR might not be optimal

# We therefore choose the bonferroni correction
b1 <- p.adjust(allHWE[, 3], method = "bonferroni")
b2 <- p.adjust(allHWE[, 7], method = "bonferroni")
b3 <- p.adjust(allHWE[, 11], method = "bonferroni")
b4 <- p.adjust(allHWE[, 15], method = "bonferroni")
b5 <- p.adjust(allHWE[, 19], method = "bonferroni")
b6 <- p.adjust(allHWE[, 23], method = "bonferroni")
b7 <- p.adjust(allHWE[, 27], method = "bonferroni")
b8 <- p.adjust(allHWE[, 31], method = "bonferroni")
b9 <- p.adjust(allHWE[, 35], method = "bonferroni")


# Put bonferroni results into overal file
allHWE[, 4] <- b1
allHWE[, 8] <- b2
allHWE[, 12] <- b3
allHWE[, 16] <- b4
allHWE[, 20] <- b5
allHWE[, 24] <- b6
allHWE[, 28] <- b7
allHWE[, 32] <- b8
allHWE[, 36] <- b9

#~~ Save the results

write.csv(allHWE, here("Figs-Tables", "HWE", "AllLocsHWE_inclTrop_2019.csv"), quote=F, row.names = F)

#~~ The software STRUCTURE assumes loci are in HWE, so before running STRUCTURE we need to remove the loci that are not in HWE
#~~ Remove loci that are not in HWE for pop9 (A. tropicalis): OrrFCB7 (L11), Lw10(L14), SSL301(L18) and ZcwD02 (L27)
#~~ Also remove locus Mang44(L38), as it has lots of NAs for A. tropicalis
#~~ We accept loci that are not in HWE in max 3 out of 9 pops, so we keep locus Agt10 (L20)

# Load structure file
struc <- fread(here("Data", "Processed", "structure_in_trop.stru"))

struc1 <- struc

# Remove loci that are not in HWE or have a lot of NAs

# struc1$OrrFCB7 <- NULL
# struc1$Lw10 <- NULL
# struc1$SSL301 <- NULL
# struc1$ZcwD02 <- NULL
# struc1$Mang44 <- NULL

struc1$V13 <- NULL
struc1$V16 <- NULL
struc1$V20 <- NULL
struc1$V29 <- NULL
struc1$V40 <- NULL

# Save new structure file
fwrite(struc1, here("Data", "Processed", "structure_trop_red.stru"), col.names = F,
       quote = F, row.names = F, sep = " ")



##########################################################
#   Test for LD in Antarctic and subantarctic samples  ###
##########################################################

# First convert structure file to genepop in PGDSpider!!!
# Use SPID_struct_genepop_39msats spid file for conversion settings

genepop::test_LD(here("Data", "Processed", "genepop_trop.txt"), here("Figs-Tables", "LD", "LD_trop.txt"), settingsFile = "",
                 dememorization = 1000, batches = 100, iterations = 1000,
                 verbose = interactive()) #NB doesnt seem to overwrite output file if it already exists...?

ld.out <- read.table(here("Figs-Tables", "LD", "LD_trop.txt"), skip=14, nrows=6669, fill=T, na.strings = c("No", "information"))

colnames(ld.out) <- c("pop", "locus_n", "locus_m", "p-value", "se", "switches")


pop.list <- split(ld.out, ld.out$pop)
rownames(pop.list[[1]]) <- paste0(pop.list[[1]][[2]], "_", pop.list[[1]][[3]])
rownames(pop.list[[2]]) <- paste0(pop.list[[2]][[2]], "_", pop.list[[2]][[3]])
rownames(pop.list[[3]]) <- paste0(pop.list[[3]][[2]], "_", pop.list[[3]][[3]])
rownames(pop.list[[4]]) <- paste0(pop.list[[4]][[2]], "_", pop.list[[4]][[3]])
rownames(pop.list[[5]]) <- paste0(pop.list[[5]][[2]], "_", pop.list[[5]][[3]])
rownames(pop.list[[6]]) <- paste0(pop.list[[6]][[2]], "_", pop.list[[6]][[3]])
rownames(pop.list[[7]]) <- paste0(pop.list[[7]][[2]], "_", pop.list[[7]][[3]])
rownames(pop.list[[8]]) <- paste0(pop.list[[8]][[2]], "_", pop.list[[8]][[3]])
rownames(pop.list[[9]]) <- paste0(pop.list[[9]][[2]], "_", pop.list[[9]][[3]])


#~~ Check p values for all 9 populations

hist(pop.list[[1]][[4]])                     
hist(pop.list[[2]][[4]]) 
hist(pop.list[[3]][[4]]) 
hist(pop.list[[4]][[4]]) 
hist(pop.list[[5]][[4]]) 
hist(pop.list[[6]][[4]]) 
hist(pop.list[[7]][[4]]) 
hist(pop.list[[8]][[4]]) 
hist(pop.list[[9]][[4]])

# Bonferroni corrections
adjust.p<-data.frame(matrix(NA, nrow = 741, ncol = 9))
adjust.p$X1 <- p.adjust(pop.list[[1]][[4]], method = "bonferroni")
adjust.p$X2 <- p.adjust(pop.list[[2]][[4]], method = "bonferroni")
adjust.p$X3 <- p.adjust(pop.list[[3]][[4]], method = "bonferroni")
adjust.p$X4 <- p.adjust(pop.list[[4]][[4]], method = "bonferroni")
adjust.p$X5 <- p.adjust(pop.list[[5]][[4]], method = "bonferroni")
adjust.p$X6 <- p.adjust(pop.list[[6]][[4]], method = "bonferroni")
adjust.p$X7 <- p.adjust(pop.list[[7]][[4]], method = "bonferroni")
adjust.p$X8 <- p.adjust(pop.list[[8]][[4]], method = "bonferroni")
adjust.p$X9 <- p.adjust(pop.list[[9]][[4]], method = "bonferroni")

rownames(adjust.p) <- rownames(pop.list[[1]])
colnames(adjust.p) <- c("SS", "BI", "CI", "HI", "KI", "MI", "trop", "MAR", "SG")

# Change all values over 0.05 to 1 for more clear plots
alpha  <- 0.05
newmat <- adjust.p
newmat[newmat > alpha] <- 1

# All non sig values are turned in to 1. Since each column contains a populations, the sum of the row will be 9 if all pops are non-significant
# So if we extract all rows that have a smaller sum than 9, we find all rows with at least one sig population
new_DF1 <- newmat[rowSums(newmat, na.rm=T) < 9,]
new_DF2 <- newmat[rowSums(newmat, na.rm=T) < 8,]
new_DF3 <- newmat[rowSums(newmat, na.rm=T) < 7,]

# Plot for quick view which pops have sig results
library("lattice")
levelplot(t(new_DF1))
levelplot(t(new_DF2))

pdf(here("Figs-Tables", "LD", "LD-trop.pdf"))
print(levelplot(t(new_DF1), aspect="fill"))
dev.off()

#~~ locus32-locus33 are in LD in 4 pops (CI, MI, MAR, SG), the rest in 3 pops or less  
