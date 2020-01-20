#------------------------------------------------------------
# AJ Paijmans
# Last edited: Apr 2019
#
# Seal microsat data - Rarify and calculate summary statistics
#                      sample size: 181 (=min cluster size), n loci=39
#
#                      Below: Sum stats with sample size 10 for fig3
#                      So that it is visually comparable with the fig per locality (S4)
#
#------------------------------------------------------------

###########################################
# Sum stats for ABC (sample size = 181) ####

# Read in data
seal <- read.table("all_gaz.txt", header = T, row.names=1, stringsAsFactors = F)

# Make subsets for each cluster
library(dplyr)

ss <- seal %>%
  filter(pop == 1)
sg <- seal %>%
  filter(pop == 2)
bi <- seal %>%
  filter(pop == 3)
mar.ci <- seal %>%
  filter(pop == 4 | pop == 5)
ki.hi.mac <- seal %>%
  filter(pop == 6 | pop == 7 | pop == 8)


# use slightly adapted function from sealABC::mssumstats to calculate sum stats, and return not only mean and SD but the complete dataset (ie all nresamp)

library(sealABC)
source("mssumstatsAP.R") #function links to other functions within sealABC package, so sealABC has to be loaded

ss_stats <- mssumstatsAP(ss, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 181, nloc = 39)
sg_stats <- mssumstatsAP(sg, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 181, nloc = 39)
bi_stats <- mssumstatsAP(bi, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 181, nloc = 39)
mar.ci_stats <- mssumstatsAP(mar.ci, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 181, nloc = 39)
ki.hi.mac_stats <- mssumstatsAP(ki.hi.mac, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 181, nloc = 39)

# Put means/sd/ci together in one data frame
sumstats_full <- rbind(ss_stats[[1]], sg_stats[[1]], bi_stats[[1]], mar.ci_stats[[1]], ki.hi.mac_stats[[1]])
row.names(sumstats_full) <- c("ss", "sg", "bi", "mar.ci", "ki.hi.mac")

# Save file
write.table(sumstats_full, "all_sumstats_full_all_gaz_5cluster.txt", quote = F, row.names = T)

# # Put all 1000 calculations per population in 1 dataframe
# 
# ss_stats[[2]]$pop <- "SSI"
# sg_stats[[2]]$pop <- "SG"
# bi_stats[[2]]$pop <- "BI"
# mar.ci_stats[[2]]$pop <- "MAR_CI"
# ki.hi.mac_stats[[2]]$pop <- "KI_HI_MAC"
# 
# 
# all_sumstats_full <- rbind(ss_stats[[2]], sg_stats[[2]], bi_stats[[2]], mar.ci_stats[[2]], ki.hi.mac_stats[[2]])
# 
# # Save file
# write.table(all_sumstats_full, "all_sumstats_5cluster_full_1000.txt", quote = F, row.names = T)


############################################
# Sum stats for plot (sample size = 10) ####
# Used for sina plots (ie fig3), to be visually comparable with fig S4

ss_stats <- mssumstatsAP(ss, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
sg_stats <- mssumstatsAP(sg, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
bi_stats <- mssumstatsAP(bi, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
mar.ci_stats <- mssumstatsAP(mar.ci, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
ki.hi.mac_stats <- mssumstatsAP(ki.hi.mac, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)

# Put all 1000 calculations per population in 1 dataframe

ss_stats[[2]]$pop <- "SSI"
sg_stats[[2]]$pop <- "SG"
bi_stats[[2]]$pop <- "BI"
mar.ci_stats[[2]]$pop <- "MAR_CI"
ki.hi.mac_stats[[2]]$pop <- "KI_HI_MAC"

all_sumstats_full <- rbind(ss_stats[[2]], sg_stats[[2]], bi_stats[[2]], mar.ci_stats[[2]], ki.hi.mac_stats[[2]])

# Save file
write.table(all_sumstats_full, "all_sumstats_5cluster_full_1000_10indv.txt", quote = F, row.names = T)
