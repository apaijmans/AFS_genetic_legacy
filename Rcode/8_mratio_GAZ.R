#------------------------------------------------------------
# Anneke Paijmans
# Last edited: Mar 2019
#
# Seal microsat data - Calculate Mratio (mean and SD)
#                      1. per population 
#                      2. per cluster
#
#------------------------------------------------------------

# Set working directory
library(here)

library(tidyverse)
library(strataG)

#####################################################
#   Mratio 1: all data per population (locality)  ###
#####################################################

#~~ Load data (created in: a_prep_emp_data_8pop, folder: fsc_pop2019)

seal <- read.table(here("Analyses", "8 ABC", "fsc_pop2019", "all_gaz.txt"), header = T, row.names=1, stringsAsFactors = F)

#~~ Convert to StrataG file

g_types_geno <- strataG::df2gtypes(as.data.frame(seal),
                                   ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)

#~~ Calculate Mratio

t <- mRatio(g_types_geno, by.strata = TRUE, rpt.size = 8:1)


#~~ Get means and SD over cols

colMeans(t)

mratio_mean_SD <- rbind(apply(as.data.frame(t), 2, mean), apply(as.data.frame(t), 2, sd))


#~~ Save file

readr::write_delim(x =  round(as.data.frame(mratio_mean_SD), 3), path = "Figs-Tables/Mratio/mratio_pop.txt", col_names = T)


###############################################
#   Mratio 2: all data per genetic cluster  ###
###############################################


#~~ Load data (created in: a_prep_emp_data_5cluster, folder: fsc_cluster2019)

seal <- read.table(here("Analyses", "8 ABC", "fsc_cluster2019", "all_gaz.txt"), header = T, row.names=1, stringsAsFactors = F)

# Relabel populations to clusters
seal <- seal %>%
  mutate(pop = gsub(5, 4, pop),
         pop = gsub(6, 5, pop),
         pop = gsub(7, 5, pop),
         pop = gsub(8, 5, pop))


#~~ Convert to StrataG file

g_types_geno <- strataG::df2gtypes(as.data.frame(seal),
                                   ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)

#~~ Calculate Mratio

t <- mRatio(g_types_geno, by.strata = TRUE, rpt.size = 8:1)


#~~ Get means and SD over cols

colMeans(t)

mratio_mean_SD <- rbind(apply(as.data.frame(t), 2, mean), apply(as.data.frame(t), 2, sd))


#~~ Save file

readr::write_delim(x =  round(as.data.frame(mratio_mean_SD), 3), path = "Figs-Tables/Mratio/mratio_cl.txt", col_names = T)
