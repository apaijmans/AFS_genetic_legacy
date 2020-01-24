#------------------------------------------------------------
# AJ Paijmans
# Last edited: June 2019
#
# Seal microsat data - Rarify and calculate summary statistics
#                      sample size: 18 (=min pop size), n loci=39
#
#                      Below: Sum stats with sample size 10 for fig S4
#                      So that it is visually comparable with the fig per cluster (fig3)
#
#------------------------------------------------------------

###########################################
# Sum stats for ABC (sample size = 18) ####

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
mar <- seal %>%
  filter(pop == 4)
ci <- seal %>%
  filter(pop == 5)
ki <- seal %>%
  filter(pop == 6)
hi <- seal %>%
  filter(pop == 7)
mac <- seal %>%
  filter(pop == 8)


# use slightly adapted function from sealABC::mssumstats to calculate sum stats, and return not only mean and SD but the complete dataset (ie all nresamp)

library(sealABC)
source("mssumstatsAP.R") #function links to other functions within sealABC package, so sealABC has to be loaded

ss_stats <- mssumstatsAP(ss, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)
sg_stats <- mssumstatsAP(sg, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)
bi_stats <- mssumstatsAP(bi, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)
mar_stats <- mssumstatsAP(mar, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)
ci_stats <- mssumstatsAP(ci, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)
ki_stats <- mssumstatsAP(ki, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)
hi_stats <- mssumstatsAP(hi, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)
mac_stats <- mssumstatsAP(mac, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 18, nloc = 39)


# Put means/sd/ci together in one data frame
sumstats_full <- rbind(ss_stats[[1]], sg_stats[[1]], bi_stats[[1]], mar_stats[[1]], ci_stats[[1]], ki_stats[[1]], hi_stats[[1]], mac_stats[[1]])
row.names(sumstats_full) <- c("ss", "sg", "bi", "mar", "ci", "ki", "hi", "mac")

# Save file
write.table(sumstats_full, "all_sumstats_full_all_gaz_8pop.txt", quote = F, row.names = T)


############################################
# Sum stats for plot (sample size = 10) ####
# Used for sina plots (fig S4), to be visually comparable with fig3

ss_stats <- mssumstatsAP(ss, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
sg_stats <- mssumstatsAP(sg, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
bi_stats <- mssumstatsAP(bi, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
mar_stats <- mssumstatsAP(mar, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
ci_stats <- mssumstatsAP(ci, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
ki_stats <- mssumstatsAP(ki, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
hi_stats <- mssumstatsAP(hi, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)
mac_stats <- mssumstatsAP(mac, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)

# Put all 1000 calculations per population in 1 dataframe
ss_stats[[2]]$pop <- "SSI"
sg_stats[[2]]$pop <- "SG"
bi_stats[[2]]$pop <- "BI"
mar_stats[[2]]$pop <- "MAR"
ci_stats[[2]]$pop <- "CI"
ki_stats[[2]]$pop <- "KI"
hi_stats[[2]]$pop <- "HI"
mac_stats[[2]]$pop <- "MAQ"

all_sumstats_full <- rbind(ss_stats[[2]], sg_stats[[2]], bi_stats[[2]], mar_stats[[2]], ci_stats[[2]], ki_stats[[2]], hi_stats[[2]], mac_stats[[2]])

# Save file
write.table(all_sumstats_full, "all_sumstats_8pops_full_1000_10indv.txt", quote = F, row.names = T)
