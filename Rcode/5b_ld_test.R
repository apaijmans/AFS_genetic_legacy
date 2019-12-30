#------------------------------------------------------------
# Anneke Paijmans
# Created: Oct 2018
# Last edited: Oct 2018
#
# Seal microsat data - LD
#
#------------------------------------------------------------

# Set working directory
library(here)

# First convert structure file to genepop in PGDSpider!!!
# Use SPID_struct_genepop_39msats spid file for conversion settings

genepop::test_LD(here("Data", "Processed", "genepop_gaz.txt"), here("Figs-Tables", "LD", "LD_gaz.txt"), settingsFile = "",
                 dememorization = 1000, batches = 100, iterations = 1000,
                 verbose = interactive()) #NB doesnt seem to overwrite output file if it already exists...

ld.out <- read.table(here("Figs-Tables", "LD", "LD_gaz.txt"), skip=14, nrows=5928, fill=T, na.strings = c("No", "information"))

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


#~~ Check p values for all 9 populations

hist(pop.list[[1]][[4]])                     
hist(pop.list[[2]][[4]]) 
hist(pop.list[[3]][[4]]) 
hist(pop.list[[4]][[4]]) 
hist(pop.list[[5]][[4]]) 
hist(pop.list[[6]][[4]]) 
hist(pop.list[[7]][[4]]) 
hist(pop.list[[8]][[4]]) 

# Bonferroni corrections
adjust.p<-data.frame(matrix(NA, nrow = 741, ncol = 8))
adjust.p$X1 <- p.adjust(pop.list[[1]][[4]], method = "bonferroni")
adjust.p$X2 <- p.adjust(pop.list[[2]][[4]], method = "bonferroni")
adjust.p$X3 <- p.adjust(pop.list[[3]][[4]], method = "bonferroni")
adjust.p$X4 <- p.adjust(pop.list[[4]][[4]], method = "bonferroni")
adjust.p$X5 <- p.adjust(pop.list[[5]][[4]], method = "bonferroni")
adjust.p$X6 <- p.adjust(pop.list[[6]][[4]], method = "bonferroni")
adjust.p$X7 <- p.adjust(pop.list[[7]][[4]], method = "bonferroni")
adjust.p$X8 <- p.adjust(pop.list[[8]][[4]], method = "bonferroni")

rownames(adjust.p) <- rownames(pop.list[[1]])
colnames(adjust.p) <- c("SS", "BI", "CI", "HI", "KI", "MI", "MAR", "SG")

# Change all values over 0.05 to 1 for more clear plots
alpha  <- 0.05
newmat <- adjust.p
newmat[newmat > alpha] <- 1

# All non sig values are turned in to 1. Since each column contains a populations, the sum of the row will be 9 if all pops are non-significant
# So if we extract all rows that have a smaller sum than 9, we find all rows with at least one sig population
new_DF1 <- newmat[rowSums(newmat, na.rm=T) < 8,] #to find any that are not 1
new_DF2 <- newmat[rowSums(newmat, na.rm=T) < 7,] #to find two that are not 1

# Plot
library("lattice")
levelplot(t(new_DF2))


#~~ Locus17-34 in LD for MI, SG and HI, all others 2 pops or less
