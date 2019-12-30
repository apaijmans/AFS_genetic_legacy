#------------------------------------------------------------
# Anneke Paijmans
# Created: June 2017
# Last edited: Mar 2019
#
# Seal microsat data - PCA Antarctic, Subantarctic and hybrid samples
#
#------------------------------------------------------------

# Clear all previous info 
#rm(list=ls()) # better restart R

# Set working directory
library(here)

#~~ Load seal data (dataset incl tropicalis but without loci that were not in HWE/were in LD)

library("adegenet")

seal <- read.structure(here("Data", "Processed", "structure_trop_red.stru"), n.ind = 2155, n.loc = 34, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

seal_df <- genind2df(seal, sep = "/")


#~~ Load hybrid data (generated in: 1e_extract_hyb)

hybs <- read.table(here("Data", "Raw", "list_of_hybrids.txt"), header = TRUE)
hybs$sampleID <- as.character(hybs$sampleID)

#Add column with species info
hybs$sp <-ifelse(hybs$cluster1 > 0.90, "A.Tropicalis", "Hybrid")

# Replace spaces in tissue name to compare hybrids list with full data set
seal_df$sampleID <- substring(row.names(seal_df), 1,11)


#~~ Merge hybrid info into seal data

seal1 <- merge(seal_df, hybs[, c(1,4)], by = "sampleID", all.x=T) # test <- seal1$sp[which(is.na(seal1$sp)==F)]
seal1$sp[which(is.na(seal1$sp)==T)] <- "A.gazella"

seal1 <- seal1[order(seal1$sampleID),] 
seal1 <- seal1[order(seal1$pop),] 

library(tidyverse)
seal1 %>% group_by(pop) %>% count(sp)
#    pop   sp               n
# <fct> <chr>        <int>
#  1 1     A.gazella      197
#  2 2     A.gazella     1042
#  3 3     A.gazella      396
#  4 3     Hybrid           1
#  5 4     A.gazella      166
#  6 4     A.Tropicalis    25
#  7 4     Hybrid           1
#  8 5     A.gazella       18
#  9 5     A.Tropicalis    23
# 10 5     Hybrid           1
# 11 6     A.gazella       51
# 12 6     Hybrid           1
# 13 7     A.gazella       22
# 14 8     A.gazella      108
# 15 8     A.Tropicalis     1
# 16 8     Hybrid          11
# 17 9     A.Tropicalis    86
# 18 9     Hybrid           5

# Save "pop ID" (actually species in this case) in seperate file
spID <- seal1$sp

seal1$sampleID <- NULL
seal1$pop <- NULL
seal1$sp <- NULL

# Save as genind object
obj.all <- df2genind(seal1, ploidy=2, sep = "/")
pop(obj.all) <- spID
obj.all

# Replace NA with mean for PCA
X <- scaleGen(obj.all, NA.method="mean")


#~~ PCA

#pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE) 
#based on the sreeplot I decided to keep the first 3 axis, since those showed the biggest change
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)

# Proportion of variance explained per PC
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent[1:3]

# 4.44 1.38 1.17, sum: 6.99
# total number of eigenvalues: 454


#~~ Plots of different axes

library(extrafont)

# Register fonts for Windows bitmap output
loadfonts(device="win")

# Define colours
#col <- c("#9999FF","#767676", "#2121D9") # blues
col <- c("#969696", "#525252", "#cccccc") # greys

# Plot axes 1-2, save as tiff
tiff(file=here("Figs-Tables", "PCA", "PCA colorplot axes 1-2 all data_grey.tif"), 3.3, 2.5, 
     units = "in", res = 600, family="Arial")
s.class(pca1$li, pop(obj.all),xax=1,yax=2, col=transp(col,1), clabel=0, # for blues use transp 0.8
        axesell=FALSE, cstar = 0, cellipse = 0, cpoint=0.85, grid=FALSE)
#add.scatter.eig(pca1$eig[1:20], 3,1,2, posi = "bottomright", csub = 0.95)
par(xpd=TRUE) # to make legend overlaying the margin
text(-4.2, 13.8, "(a)", cex = 1, family="Arial")
text(-3.8, 0.9, "PC1", cex = 0.85)
text(1.6, 13.9, "PC2", cex = 0.85)
dev.off()

# Plot axes 1-3 and add legend (species), save as tiff
species <- c(expression(italic("A. gazella")),
             "Hybrids",
             expression(italic("A. tropicalis")))

tiff(file=here("Figs-Tables", "PCA", "PCA colorplot axes 1-3 all data_grey.tif"), 3.3, 2.5, 
     units = "in", res = 600, family="Arial")
s.class(pca1$li, pop(obj.all),xax=1,yax=3, col=transp(col,1), clabel=0, # for blues use transp 0.8
        axesell=FALSE, cstar = 0, cellipse = 0, cpoint=0.85, grid=FALSE)
#add.scatter.eig(pca1$eig[1:20], 3,1,3, posi = "bottomright", csub = 0.95)
par(xpd=TRUE) # to make legend overlaying the margin
legend(-34.3,-6.5, legend=species, pch=19, col=transp(col,.8), pt.cex = 0.6, cex = .85, y.intersp=0.25, x.intersp=0.5, bty="n")
text(-28.5, 32.5, "(b)", cex = 1, family="Arial")
text(-27.3, 2.2, "PC1", cex = 0.85)
text(4.6, 32.6, "PC3", cex = 0.85)
dev.off()
