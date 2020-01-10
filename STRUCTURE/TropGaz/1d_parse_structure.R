# Clear all previous info
#rm(list=ls())

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
  labs(x = "K", y = expression(paste(Delta,italic("K")))) +
  scale_x_continuous(breaks=c(1:10), labels=c(1:10),limits=c(1,10)) +
  ggtitle('A') + theme(plot.title=element_text(hjust=0, size = 20))

elpd <- ggplot(ks, aes(x=k, y = elpdmean)) +
  geom_point(size = 1, col = "grey30") + # 1/1.5
  geom_line(size = 1, col = "grey30") +
  geom_errorbar(aes(ymin = elpdmean - elpdsd, ymax= elpdmean + elpdsd), colour="grey30", width=0) +
  theme(axis.text.x = element_text(face = "plain", size=20),
        axis.text.y = element_text(face = "plain", size=20),
        axis.title = element_text(size = 20)) +
  labs(x = "K", y = expression(paste("Ln Pr(",italic("X"),"|",italic("K"),")"))) +
  scale_x_continuous(breaks=c(1:10), labels=c(1:10),limits=c(1,10)) +
  ggtitle('B') + theme(plot.title=element_text(hjust=0, size = 20))


png("deltaK_elpd.png", units = "in", res = 300, width = 15, height = 9)
grid.arrange(dk, elpd, ncol=2, nrow=1)
dev.off()

system("mv deltaK_elpd.png results/run_files/evanno")


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
        outputfilename = paste0(path_to_structure_out,"struct-plots/K2_reps"))
  
  plotQ(qlist=(slist)[21:30], splab=spnames[21:30], imgoutput = "join", grplab=pops,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K3_reps"))
  
}

make_structure_plots("results/run_files/", "structure_trop_red.stru")

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
  
  
  plotQ(qlist=(clist)[1:2], splab=spnames[c(11,21)], imgoutput = "join", grplab=pops, #showindlab = T,
        #clustercol=c("#1B9E77", "#D95F02", "#E6AB02", "#7570B3", "#E7298A", "#66A61E"),
        outputfilename = paste0(path_to_structure_out,"struct-plots/K2-K3_merged"))
  
}

make_structure_clumpp_plots("results/run_files/", "structure_trop_red.stru")

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
                                                                 ifelse(V2== 8, "MAC", "TM")))))))))
  pops <- data.frame(pops[c(T,F),])
  colnames(pops) <- "location"
  pops$location <- as.character(pops$location)
  
  plotQ(qlist=(clist)[1], splab=spnames[c(11)], splabsize=10, grplab=pops, grplabsize=4.5, grplabheight=1, 
        showlegend = F, #legendkeysize = 10, legendtextsize = 6,
        clustercol=c("#cccccc", "#969696"),
        height = 4, width=40, # if white lines appear, increase width
        outputfilename = paste0(path_to_structure_out,"struct-plots/K2_merged1_grey"))
  
}

make_structure_clumpp_plot_k("results/run_files/", "structure_trop_red.stru")
