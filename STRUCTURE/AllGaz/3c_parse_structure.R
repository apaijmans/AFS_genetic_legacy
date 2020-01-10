# Clear all previous info
# rm(list=ls())

# Parse and plot STRUCTURE output using mainly pophelper()

library(dplyr)
library(tidyr)
library(stringr)
# install.packages('devtools',dependencies=T)
# library(devtools)
# install_github("royfrancis/pophelper", force = T)
library(pophelper)
library(data.table)
options(scipen=999)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Collect output files                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Full

system("mkdir results/run_files")
system("mv results/*_f results/run_files/")

system("mkdir results/run_files/evanno")
system("mkdir results/run_files/struct-plots")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Load files and collect clumpp output            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load_and_clumpp <- function(path_to_structure_out){
  
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  # export for CLUMPP later on
  clumppExport(qlist=slist, parammode=3, useexe=T)# if CLUMPP takes to long set parammode to 3, otherwise leave out
  # collect CLUMPP output
  collectClumppOutput(filetype="merged") # aligned
  system("rm -r pop_K*")
  
  # move clump files to correct directory
  system(paste0("mv pop-* ", path_to_structure_out))
  
}


load_and_clumpp("results/run_files/")

# Function to get K summary stats from run files

load_and_K <- function(path_to_structure_out){
  
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  # Get summary stats
  
  em <- evannoMethodStructure(summariseQ(tabulateQ(slist)))
  plot(em$deltaK)
  
  evannoMethodStructure(data=em, exportplot=T, writetable = T)
  system(paste0("mv evannoMethodStructure* ", path_to_structure_out, "evanno/"))
  
  em <- em
  
}

ks <- load_and_K("results/run_files/")


#~~ optimal K 

# if mean estimated ln probability of data is highest at K=1, K=1, 
# if not, use evanno-method (delta K) to choose k

optimal_k <- function(x) {
  elpd <- which.max(x$elpdmean)
  if (elpd == 1) {
    return(1)
  } else {
    delta_k <- which.max(x$deltaK)
    delta_k
  }
}

optimal_k(ks)


# output best k's

both_k <- function(x) {
  elpd <- which.max(x$elpdmean)
  deltaK <- which.max(x$deltaK)
  ks <- c(lnk = x$k[elpd], deltak = x$k[deltaK])
  # write outfile
  #write.table(ks, paste0(path_to_structure_out, "Ks.txt"))
  ks <- ks
}

bothKs <- both_k(ks)


#~~ Make deltak and elpd plot for Supp

require(gridExtra)
library(ggthemr)
ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.5, text_size = 16, type = "outer", spacing = 2)


#~~ Microsatellite Ks

dk <- ggplot(ks, aes(x=k, y = deltaK)) +
  geom_point(size = 1, col = "grey30") +
  geom_line(size = 1, col = "grey30") +
  theme(axis.text.x = element_text(face = "plain", size=20),
        axis.text.y = element_text(face = "plain", size=20),
        axis.title = element_text(size = 20)) +
  labs(x = "k", y = expression(paste(Delta,italic("k")))) +
  scale_x_continuous(breaks=c(1:10), labels=c(1:10),limits=c(1,10)) +
  ggtitle('(a)') + theme(plot.title=element_text(hjust=0, size = 20))

elpd <- ggplot(ks, aes(x=k, y = elpdmean)) +
  geom_point(size = 1, col = "grey30") + # 1/1.5
  geom_line(size = 1, col = "grey30") +
  geom_errorbar(aes(ymin = elpdmean - elpdsd, ymax= elpdmean + elpdsd), colour="grey30", width=0) +
  theme(axis.text.x = element_text(face = "plain", size=20),
        axis.text.y = element_text(face = "plain", size=20),
        axis.title = element_text(size = 20)) +
  labs(x = "k", y = expression(paste("Ln Pr(",italic("X"),"|",italic("k"),")"))) +
  scale_x_continuous(breaks=c(1:10), labels=c(1:10),limits=c(1,10)) +
  ggtitle('(b)') + theme(plot.title=element_text(hjust=0, size = 20))


png("deltak_elpd.png", units = "in", res = 300, width = 15, height = 9)
grid.arrange(dk, elpd, ncol=2, nrow=1)
dev.off()

system("mv deltak_elpd.png results/run_files/evanno")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Structure Plots               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#devtools::install_github('cttobin/ggthemr')

#~~ Structure runs not merged:

make_structure_plots <- function(path_to_structure_out, path_to_struc_file){
  
  # Get slist
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure")
  
  
  #~~ Customize plots
  
  # strip panel label showing k only
  fn1 <- function(x) attr(x,"k")
  spnames <- paste0("K=",sapply(slist,fn1))
  
  # custom colours
  library(ggthemr)
  ggthemr(palette = "solarized", layout = "clean",
          line_weight = 0.7, text_size = 20, type = "outer")
  swatch()
  
  # add group labels
  struc_file <- path_to_struc_file
  
  pops <- fread(struc_file) %>%
    select(V2) %>%
    mutate(V2 = ifelse(V2 == 1, "SShetlands",
                       ifelse(V2 == 2, "SGeorgia",
                              ifelse(V2 == 3, "Bouvetoya",
                                     ifelse(V2 == 4, "Marion", 
                                            ifelse(V2 == 5, "Crozet",
                                                   ifelse(V2 == 6, "Kerguelen",
                                                          ifelse(V2 == 7, "Heard", 
                                                                 ifelse(V2== 8, "Macquarie", "Trop_Mac")))))))))
  
  pops <- data.frame(pops[c(T,F),])
  colnames(pops) <- "location"
  pops$location <- as.character(pops$location)
  
  
  #~~ Write out plots
  
  #correctly align plot panel and label panel
  
  plotQ(qlist=(slist)[1:10], splab=spnames[1:10], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K1_reps"))
  
  plotQ(qlist=(slist)[11:20], splab=spnames[11:20], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K10_reps"))
  
  plotQ(qlist=(slist)[21:30], splab=spnames[21:30], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K2_reps"))
  
  plotQ(qlist=(slist)[31:40], splab=spnames[31:40], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K3_reps"))
  
  plotQ(qlist=(slist)[41:50], splab=spnames[41:50], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K4_reps"))
  
  plotQ(qlist=(slist)[51:60], splab=spnames[51:60], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K5_reps"))
  
  plotQ(qlist=(slist)[61:70], splab=spnames[61:70], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K6_reps"))
  
  plotQ(qlist=(slist)[71:80], splab=spnames[71:80], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K7_reps"))
  
  plotQ(qlist=(slist)[81:90], splab=spnames[81:90], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K8_reps"))
  
  plotQ(qlist=(slist)[91:100], splab=spnames[91:100], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K9_reps"))
  
}

make_structure_plots("results/run_files/", "structure_in_gaz.stru")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. CLUMPP Structure Plots        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Plot with merged files

make_structure_clumpp_plots <- function(path_to_structure_out, path_to_struc_file){
  
  # Get slist
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure", indlabfromfile=T)
  # File paths
  clumpp_files <- list.files(paste0(path_to_structure_out,"pop-merged/"))
  clumpp_out_paths <- paste0(path_to_structure_out,"pop-merged/", clumpp_files)
  
  # Read in files
  clist <- readQ(files=clumpp_out_paths, filetype = "clumpp")
  
  # Plot aligned structure plots
  
  # strip panel label showing k only
  fn1 <- function(x) attr(x,"k")
  spnames <- paste0("K=",sapply(slist,fn1))
  
  struc_file <- path_to_struc_file
  
  pops <- fread(struc_file) %>%
    select(V2) %>%
    mutate(V2 = ifelse(V2 == 1, "SShetlands",
                       ifelse(V2 == 2, "SGeorgia",
                              ifelse(V2 == 3, "Bouvetoya",
                                     ifelse(V2 == 4, "Marion", 
                                            ifelse(V2 == 5, "Crozet",
                                                   ifelse(V2 == 6, "Kerguelen",
                                                          ifelse(V2 == 7, "Heard", 
                                                                 ifelse(V2== 8, "Macquarie", "Trop_Mac")))))))))
  pops <- data.frame(pops[c(T,F),])
  colnames(pops) <- "location"
  pops$location <- as.character(pops$location)
  
  cols <- c( "#B2DF8A", "#6A3D9A", "#FFA500", "#E7298A", "#33A02C", "#1F78B4", "#B15928", "#E31A1C")
  
  plotQ(qlist=(clist)[3:5], splab=spnames[c(31,41,51)], splabsize=10, imgoutput = "join", grplab=pops, #showindlab = T,
        clustercol = cols, # c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K3-K5_merged"))
  
}

make_structure_clumpp_plots("results/run_files/", "structure_in_gaz.stru")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6b. CLUMPP Structure Plot for specific K        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Plot with merged files

make_structure_clumpp_plot_k <- function(path_to_structure_out, path_to_struc_file){
  
  library("adegenet")
  
  # Get slist
  all_files <- list.files(path_to_structure_out, pattern = "^results")
  
  # creates character vector containing path to relevant files
  struc_out_paths <- paste0(path_to_structure_out, all_files)
  
  slist <- readQ(files=struc_out_paths, filetype = "structure", indlabfromfile=T)
  # File paths
  clumpp_files <- list.files(paste0(path_to_structure_out,"pop-merged/"))
  clumpp_out_paths <- paste0(path_to_structure_out,"pop-merged/", clumpp_files)
  
  # Read in files
  clist <- readQ(files=clumpp_out_paths, filetype = "clumpp")
  
  # Plot aligned structure plots
  
  # strip panel label showing k only
  fn1 <- function(x) attr(x,"k")
  spnames <- paste0("K=",sapply(slist,fn1))
  
  struc_file <- path_to_struc_file
  
  pops <- fread(struc_file) %>%
    select(V2) %>%
    mutate(V2 = ifelse(V2 == 1, "SSI",
                       ifelse(V2 == 2, "SG",
                              ifelse(V2 == 3, "BI",
                                     ifelse(V2 == 4, "MAR", 
                                            ifelse(V2 == 5, "CI",
                                                   ifelse(V2 == 6, "KI",
                                                          ifelse(V2 == 7, "HI", 
                                                                 ifelse(V2== 8, "MAC", "Trop Macq")))))))))
  pops <- data.frame(pops[c(T,F),])
  colnames(pops) <- "location"
  pops$location <- as.character(pops$location)
  
  #Colours:
  #B2DF8A light green = BI
  #6A3D9A purple = SG
  #FFA500 orange = MAC
  #E7298A pink = SS
  #33A02C dark green = MAR
  #1F78B4 blue = CI
  #B15928 dark orange = KI
  #E31A1C red = HI
  
  
  plotQ(qlist=(clist)[2], splab=spnames[c(21)], splabsize=15, grplab=pops, grplabsize=3,
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#FFA500", "#6A3D9A"),
        height = 4, width=30,
        outputfilename = paste0(path_to_structure_out,"struct-plots/K2_merged"))
  
  plotQ(qlist=(clist)[3], splab=spnames[c(31)], splabsize=15, grplab=pops, grplabsize=3,
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#B2DF8A", "#6A3D9A", "#FFA500"),
        height = 4, width=30,
        outputfilename = paste0(path_to_structure_out,"struct-plots/K3_merged"))
  
  plotQ(qlist=(clist)[4], splab=spnames[c(41)], splabsize=15, grplab=pops, grplabsize=3,
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#B2DF8A", "#E7298A", "#FFA500", "#6A3D9A"), # cols,
        height = 4, width=30,
        outputfilename = paste0(path_to_structure_out,"struct-plots/K4_merged"))
  
  plotQ(qlist=(clist)[5], splab=spnames[c(51)], splabsize=15, grplab=pops, grplabsize=3,
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#33A02C", "#E7298A", "#FFA500", "#6A3D9A", "#B2DF8A"),
        height = 4, width=30,
        outputfilename = paste0(path_to_structure_out,"struct-plots/K5_merged"))
  
  plotQ(qlist=(clist)[6], splab=spnames[c(61)], splabsize=15, grplab=pops, grplabsize=3,
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#1F78B4", "#E7298A", "#33A02C", "#FFA500", "#B2DF8A",  "#6A3D9A"),
        height = 4, width=30,
        outputfilename = paste0(path_to_structure_out,"struct-plots/K6_merged"))
  
  plotQ(qlist=(clist)[7], splab=spnames[c(71)], splabsize=15, grplab=pops, grplabsize=3,
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#B15928", "#B2DF8A", "#33A02C", "#1F78B4", "#E7298A", "#6A3D9A", "#FFA500"),
        height = 4, width=30,
        outputfilename = paste0(path_to_structure_out,"struct-plots/K7_merged"))
  
  plotQ(qlist=(clist)[8], splab=spnames[c(81)], splabsize=15, grplab=pops, grplabsize=4.5, grplabheight=1,
        linesize=0.7, pointsize=4, linepos=.7,
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#FFA500", "#E31A1C", "#1F78B4", "#33A02C", "#B2DF8A", "#B15928", "#E7298A", "#6A3D9A"),
        height = 4, width=30,
        outputfilename = paste0(path_to_structure_out,"struct-plots/K8_merged"))
  
}

make_structure_clumpp_plot_k("results/run_files/", "structure_in_gaz.stru")
