#------------------------------------------------------------
# Anneke Paijmans
# Created: June 2017
# last edited: Mar 2019
#
# Seal microsat data - PCA
#                      1. all data per population (locality)
#                      2. equal sample size per population
#
#------------------------------------------------------------

# Clear all previous info
#rm(list=ls()) # better restart R

# Set working directory and load libraries
library(here)
library("adegenet")

# All populations
pops <- c("South Shetlands", 
          "South Georgia", 
          "Bouvetøya", 
          "Marion Island", 
          "Crozet Islands", 
          "Kerguelen Islands", 
          "Heard Island", 
          "Macquarie Island")


##################################################
#   PCA 1: all data per population (locality)  ###
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


#~~ Simple PCA analysis

# First deal with NAs (replaced in this case by the mean allele frequency, scalegen also gives the option to centre data and scale data))
X <- scaleGen(seal, NA.method="mean")
#class(X)

#X[1:5,1:5]


#~~ PCA

#pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE) 
#based on the sreeplot I decided to keep the first 3 axis, since those showed the biggest change
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)


#~~ Eigenvalue plot absolute variance

barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

# Proportion of variance explained per PC
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent[1:3]
# 1.84 1.22 0.91, sum(eig_percent[1:3]): 3.97, total PCs: 406


#~~ Plots of different axes

library(extrafont)

# Define colours
col <- c( "#E7298A","#6A3D9A","#B2DF8A", "#33A02C", "#1F78B4", "#B15928", "#E31A1C", "#FFA500")

# Plot axes 1-2, save as tiff
tiff(file=here("Figs-Tables", "PCA", "PCA colorplot axes 1-2_gaz.tiff"), 3.3, 2.5, 
     units = "in", res = 600, family="Arial")
s.class(pca1$li, pop(seal), xax=1, yax=2, col=transp(col,.8), clabel=0, #label=c(1,2,3,4,5,6,7,8,9),  
        axesell=FALSE, cstar = 0, cellipse = 0, cpoint=0.85, grid=FALSE)
par(xpd=TRUE) # to make legend overlaying the margin
#add.scatter.eig(pca1$eig[1:20], nf=3,xax=1,yax=2, posi = "bottomright", ratio = 0.20, csub = 0.95)
text(-12.95, 7.5, "(a)", cex = 1, family="Arial")
text(-12.6, 0.7, "PC1", cex = 0.85)
text(1.1, 7.5, "PC2", cex = 0.85)
dev.off()

# Plot axes 1-3, add legend, save as tiff
tiff(file=here("Figs-Tables", "PCA", "PCA colorplot axes 1-3_gaz.tif"), 3.3, 2.5, 
     units = "in", res = 600, family="Arial")
s.class(pca1$li, pop(seal), xax=1, yax=3, col=transp(col,.8), clabel=0, #label=c(1,2,3,4,5,6,7,8,9),  
        axesell=FALSE, cstar = 0, cellipse = 0, cpoint=0.85, grid=FALSE)
par(xpd=TRUE) # to make legend overlaying the margin
legend(-17.4,15.1, pops, pch=19, col=transp(col,.8), pt.cex = 0.6, cex = .85, y.intersp=0.25, x.intersp=0.5, bty="n")
#add.scatter.eig(pca1$eig[1:20], nf=3,xax=1,yax=3, posi = "bottomright", ratio = 0.20, csub = 0.95)
text(-15.4, 13.3, "(b)", cex = 1, family="Arial")
text(-15, 0.8, "PC1", cex = 0.85)
text(1.4, 13.4, "PC3", cex = 0.85)
dev.off()


############################################################
#   PCA 2: equal sample sizes per population (locality)  ###
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


#~~ Simple PCA analysis

# First deal with NAs (replaced in this case by the mean allele frequency)
X <- scaleGen(seal, NA.method="mean")
#class(X)

#X[1:5,1:5]


#~~ PCA

#pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE) # Keep first 3
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)


#~~ Eigenvalue plot absolute variance

barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

# Proportion of variance explained per PC
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent[1:3]


#~~ Plots of different axes

# Define colours
col <- c( "#E7298A","#6A3D9A","#B2DF8A", "#33A02C", "#1F78B4", "#B15928", "#E31A1C", "#FFA500")

# Plot axes 1-2, save as tiff
tiff(file=here("Figs-Tables", "PCA", "PCA colorplot axes 1-2_equal_size.tif"), 3.3, 2.5, 
     units = "in", res = 600, family="Arial")
s.class(pca1$li, pop(seal), xax=1, yax=2, col=transp(col,.8), clabel=0,  
        axesell=FALSE, cstar = 0, cellipse = 0, cpoint=0.85, grid=FALSE)
par(xpd=TRUE)
#add.scatter.eig(pca1$eig[1:20], 3,1,2, posi = "bottomright", ratio = 0.20, csub = 0.95)
text(-18.35, 19.6, "(a)", cex = 1, family="Arial")
text(-17.8, 1.2, "PC1", cex = 0.85)
text(2.1, 19.6, "PC2", cex = 0.85)
dev.off()

# Plot axes 1-3, save as tiff
tiff(here("Figs-Tables", "PCA", "PCA colorplot axes 1-3_equal_size.tif"), 3.3, 2.5, units = "in", res = 600, family="Arial")
s.class(pca1$li, pop(seal), xax=1, yax=3, col=transp(col,.8), clabel=0,  
        axesell=FALSE, cstar = 0, cellipse = 0, cpoint=0.85, grid=FALSE)
par(xpd=TRUE)
#add.scatter.eig(pca1$eig[1:20], 3,1,3, posi = "bottomright", ratio = 0.20, csub = 0.95)
text(-12.55, 10, "(b)", cex = 1, family="Arial")
text(-12.2, 1, "PC1", cex = 0.85)
text(1.5, 10, "PC3", cex = 0.85)
dev.off()

# Plot legend separatly, save as tiff
tiff(here("Figs-Tables", "PCA", "PCA colorplot equal_size legend.tif"), 5, 2, units = "in", res = 600, family="Arial")
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("left", c("SSI","SG","BI","MAR","CI","KI","HI","MAC"), 
       pch=19, col=transp(col,.8), pt.cex = 0.6, cex = 0.85,  x.intersp=0.3, y.intersp=0.25,
       xpd=TRUE, horiz=TRUE, bty="n",
       text.width=c(0, .12, .115, .095, .115, .104, .095, .088))#sets width between text chunks in legend
dev.off()

# Assemble plots together in inkscape or any other picture editor