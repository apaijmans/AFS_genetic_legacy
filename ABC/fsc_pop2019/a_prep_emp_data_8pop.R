#------------------------------------------------------------
# AJ Paijmans
# Last edited: June 2019
#
# Seal microsat data - Prep empirical msat data for ABC
#                      
#------------------------------------------------------------

# Prep empirical msat data used in STRUCTURE for analysis

library("pegas")
library("adegenet")

# Read seal data
seal <- read.structure("structure_in_gaz.stru", n.ind = 2000, n.loc = 39, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

# number of genotypes: 2000 <- Use here number of individuals, not number of rows (per individual two rows)
# number of markers: 39
# which column labels for genotypes: 1
# which column population factor: 2
# Optional columns: enter
# Row marker names: 0 (absent)
# Genotypes coded by single row? n

namecol <- c("Pv9_a", "Pv9_b", "Hg63_a", "Hg63_b", "Hg810_a", "Hg810_b", "Hg13_a", "Hg13_b", "M11a_a", "M11a_b", "PvcA_a", "PvcA_b", "Zcwb07_a", "Zcwb07_b", "Agaz2_a", "Agaz2_b",
             "Ag3_a", "Ag3_b", "Agaz6_a", "Agaz6_b", "OrrFCB7_a", "OrrFCB7_b", "Ag2_a", "Ag2_b", "OrrFCB2_a", "OrrFCB2_b", "Lw10_a", "Lw10_b", "Zcwc01_a", "Zcwc01_b", "Agaz5_a", "Agaz5_b",
             "ZcwCgDhB14_a", "ZcwCgDhB14_b", "SSL301_a", "SSL301_b", "Ag7_a", "Ag7_b", "Agt10_a", "Agt10_b", "ZcwCgDh47_a", "ZcwCgDh47_b", "Zcwe05_a", "Zcwe05_b", "Ag1_a", "Ag1_b",
             "OrrFCB8_a", "OrrFCB8_b", "Agt47_a", "Agt47_b", "Zcwf07_a", "Zcwf07_b", "ZcwD02_a", "ZcwD02_b", "ZcwCgDh18_a", "ZcwCgDh18_b", "Aa4_a", "Aa4_b", "ZcCgDh58_a", "ZcCgDh58_b",
             "Agaz3_a", "Agaz3_b", "X9621_a", "X9621_b", "X5546_a", "X5546_b", "Zcwa12_a", "Zcwa12_b", "PvcE_a", "PvcE_b", "Zcwb09_a", "Zcwb09_b", "agaz10_a", "agaz10_b",
             "Mang44_a", "Mang44_b", "Mang36_a", "Mang36_b" )

# Convert genind object into loci file
seal1 <- genind2loci(seal)

# Split columns
seal2 <- splitstackshape::cSplit(seal1[-1], names(seal1[-1]), "/")

# Rename col names
names(seal2) <- namecol

# Merge pop column back in
seal3 <- cbind(seal1[1], seal2)

# Make rownames into column
seal3 <- cbind(rownames(seal3), seal3)
colnames(seal3)[1] <- "id"
colnames(seal3)[2] <- "pop"
row.names(seal3) <- NULL

# Clear all seal files that are not longer needed
rm("seal", "seal1", "seal2")

# Save data
write.table(seal3, "all_gaz.txt", quote = F, row.names = T)
