#------------------------------------------------------------
# Anneke Paijmans
# last edited: Mar 2019
#
# Seal microsat data - Loading and analysing result from BOTTLENECK software
#                      1. per population 
#                      2. per cluster
#
#------------------------------------------------------------


# Set working directory
library(here)


#######################################
#   1. Per populations (locality)   ###
#######################################


#~~ Get data for pops SSI, BI, MAR, CI, KI, HI
#~~ Data for SG and MAC are imported later in this file

# Files 428 & 602 are missing because due to random selection of samples 
# these files contained too many missing values for a locus in a pop, and BOTTLENECK
# crashed while running these files

# Import files
filelist = list.files(path = here("Analyses", "7 BOTTLENECK", "dec2018"), pattern="^results", full.names=T)
myfiles = lapply(filelist, read.table, fill=T, stringsAsFactors = F)
raw_file <- do.call(rbind.data.frame, myfiles)

n = 998 # Number of runs contained in the BOTTLENECK file, ie number of files analysed in a row by BOTTLENECK

# rename columns, rename populations
pop <- rep(c("SSI", "SG", "BI", "MAR", "CI", "KI", "HI", "MAC"), n)
run <- rep(1:n, each = 8)

df <- data.frame(pop, run, stringsAsFactors=FALSE)

# Extract data of interest
exp.loci.het.ex <- raw_file$V8[grep("Expected", raw_file$V1)]

loci.het.def <- raw_file$V1[grep("Expected", raw_file$V1)+1]

loci.het.ex <- raw_file$V7[grep("Expected", raw_file$V1)+1]

prop_het.ex <- as.numeric(loci.het.ex)/39

# Bind all data of interest
all_data <- cbind(df, exp.loci.het.ex, loci.het.def, loci.het.ex, prop_het.ex)

# Remove pops that we will update in the next step
all_data2 <- all_data[!grepl("SG|MAC",all_data$pop),]


#~~ Import files for SG, MAC (had to be re-run due to changes in dataset in these pops)

# Files 185 & 809 are missing because due to random selection of samples 
# these files contained too many missing values for a locus in a pop, and BOTTLENECK
# crashed while running these files

# Import files
filelist = list.files(path = here("Analyses", "7 BOTTLENECK", "march2019", "SG_MAC_pop"), pattern="^results", full.names=T)
myfiles = lapply(filelist, read.table, fill=T, stringsAsFactors = F)
raw_file_SG_MAC <- do.call(rbind.data.frame, myfiles)

n = 998# Number of runs contained in the BOTTLENECK file, ie number of files analysed in a row by BOTTLENECK

# rename columns, rename populations
pop <- rep(c("SG", "MAC"), n)
run <- rep(1:n, each = 2)

df <- data.frame(pop, run, stringsAsFactors=FALSE)

# Extract data of interest
exp.loci.het.ex <- raw_file_SG_MAC$V8[grep("Expected", raw_file_SG_MAC$V1)]

loci.het.def <- raw_file_SG_MAC$V1[grep("Expected", raw_file_SG_MAC$V1)+1]

loci.het.ex <- raw_file_SG_MAC$V7[grep("Expected", raw_file_SG_MAC$V1)+1]

prop_het.ex <- as.numeric(loci.het.ex)/39

# Bind all data of interest
all_data_SG_MAC <- cbind(df, exp.loci.het.ex, loci.het.def, loci.het.ex, prop_het.ex)


#~~ Combine 1st set of results with 2nd set of results (SG and MAC)

all_data3 <- rbind(all_data2, all_data_SG_MAC)

all_data3 %>% 
  group_by(pop) %>%
  summarise(mean_loci = mean(as.numeric(as.character(loci.het.ex))))


#~~ change to long format for making figure (fig made in: 6_one_plot)

library(tidyverse)

prophet_long <- all_data3[c("prop_het.ex", "pop")] %>%
  gather("prop_het.ex", key = "stat", value = "values") %>%
  arrange(pop)

prophet_long %>% 
  group_by(pop) %>%
  summarise(mean_alleles = mean(values))


#~~ Save data

readr::write_delim(x =  prophet_long, path = "Figs-Tables/BOTTLENECK/BOTTLENECK_1-1000_pop_long.txt", col_names = T)



################################
#   2. Per genetic cluster   ###
################################


#~~ Get data for clusters SSI, BI, MAR_CI
#~~ Data for SG and HI_KI_MAC are imported later in this file

# Import files
filelist = list.files(path = here("Analyses", "7 BOTTLENECK", "dec2018", "clusters"), pattern="^results", full.names=T)
myfiles = lapply(filelist, read.table, fill=T, stringsAsFactors = F)
raw_file <- do.call(rbind.data.frame, myfiles)

n = 1000 # Number of runs contained in the BOTTLENECK file, ie number of files analysed in a row by BOTTLENECK
#n = 1

# rename columns, rename populations
pop <- rep(c("SSI", "SG", "BI", "MAR_CI", "HI_KI_MAC"), n)
run <- rep(1:n, each = 5)

df <- data.frame(pop, run, stringsAsFactors=FALSE)

# Extract data of interest
exp.loci.het.ex <- raw_file$V8[grep("Expected", raw_file$V1)]

loci.het.def <- raw_file$V1[grep("Expected", raw_file$V1)+1]

loci.het.ex <- raw_file$V7[grep("Expected", raw_file$V1)+1]

prop_het.ex <- as.numeric(loci.het.ex)/39

# Bind all data of interest
all_data <- cbind(df, exp.loci.het.ex, loci.het.def, loci.het.ex, prop_het.ex)

# Remove clusters that we will update in the next step
all_data2 <- all_data[!grepl("SG|HI_KI_MAC",all_data$pop),]


#~~ Import files for SG, HI_KI_MAC (had to be re-run due to changes in dataset in these pops)

# Files 185 & 809 are missing because due to random selection of samples 
# these files contained too many missing values for a locus in a pop, and BOTTLENECK
# crashed while running these files

# Import files
filelist = list.files(path = here("Analyses", "7 BOTTLENECK", "march2019", "SG_MAC_CL"), pattern="^results", full.names=T)
myfiles = lapply(filelist, read.table, fill=T, stringsAsFactors = F)
raw_file_SG_MAC <- do.call(rbind.data.frame, myfiles)

n = 998# Number of runs contained in the BOTTLENECK file, ie number of files analysed in a row by BOTTLENECK

# rename columns, rename populations
pop <- rep(c("SG", "HI_KI_MAC"), n)
run <- rep(1:n, each = 2)

df <- data.frame(pop, run, stringsAsFactors=FALSE)

# Extract data of interest
exp.loci.het.ex <- raw_file_SG_MAC$V8[grep("Expected", raw_file_SG_MAC$V1)]

loci.het.def <- raw_file_SG_MAC$V1[grep("Expected", raw_file_SG_MAC$V1)+1]

loci.het.ex <- raw_file_SG_MAC$V7[grep("Expected", raw_file_SG_MAC$V1)+1]

prop_het.ex <- as.numeric(loci.het.ex)/39

# Bind all data of interest
all_data_SG_MAC <- cbind(df, exp.loci.het.ex, loci.het.def, loci.het.ex, prop_het.ex)


#~~ Combine 1st set of results with 2nd set of results (SG and MAC)

all_data3 <- rbind(all_data2, all_data_SG_MAC)


#~~ change to long format for making figure (fig made in: 6_one_plot)

prophet_long <- all_data3[c("prop_het.ex", "pop")] %>%
  gather("prop_het.ex", key = "stat", value = "values") %>%
  arrange(pop)

prophet_long %>% 
  group_by(pop) %>%
  summarise(mean_alleles = mean(values))


#~~ Save data

readr::write_delim(x =  prophet_long, path = "Figs-Tables/BOTTLENECK/BOTTLENECK_1-1000_cluster_long.txt", col_names = T)
