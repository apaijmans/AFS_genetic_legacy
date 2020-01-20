
##########################################################################################
# Summary statistics for land breeding species for multispecies plot
# using Martin Stoffel's data (https://github.com/mastoffel/pinniped_bottlenecks)


library(sealABC)
library(readr)


# bottleneck_results_30 file is needed to select the right data from the seal_data file
all_seals <- read_excel_sheets("data/seal_data_largest_clust_and_pop_30.xlsx")
bottleneck <- read_delim("data/bottleneck_results_30.txt", col_names = TRUE, delim = " ")

all_seals <- all_seals[bottleneck$id]


# select only land breeding species, excluding AFS as we replace this data with our own
land_seals <- all_seals[c("australian_fur_seal",
                          "california_sea_lion",
                          "galapagos_fur_seal",
                          "galapagos_sea_lion",
                          "grey_seal_orkneys",
                          "guadalupe_fur_seal",
                          "harbour_seal_waddensee",
                          "hawaiian_monk_seal",
                          "mediterranean_monk_seal",
                          "new_zealand_fur_seal",
                          "new_zealand_sea_lion",
                          "nes",
                          "northern_fur_seal",
                          "south_american_fur_seal",
                          "south_american_sea_lion",
                          "ses",
                          "stellers_sea_lion",
                          "subantarctic_fur_seal")]


source("mssumstatsAP.R") #function links to other functions within sealABC package, so sealABC has to be loaded

set.seed(1122)
sumstats <- lapply(land_seals, function(x) mssumstatsAP(x, start_geno = 4, mratio = "loose", rarefaction = TRUE, nresamp = 1000, nind = 10))


library(purrr)
all_sumstats <- map_df(sumstats, ~as.data.frame(.x[[2]]), .id="species") 

# .x[[2]]: list is made of sublists for each species, 
# of which the second list (sumstats[[i]][[2]]) contains the 1000 values that we need
# so with this piece of code we select the 2nd sublist for each species


readr::write_delim(all_sumstats, path = "data/all_sumstats_landseals_1000_10indv.txt", col_names = TRUE)

##########################################################################################
# Sum stats for AFS for multispecies plot 
# use entire population, ie not split per location (sample size = 10, thesame as for other species above)

# Read in data
seal <- read.table("~/fsc_pop2019/all_gaz.txt", header = T, row.names=1, stringsAsFactors = F)

# use slightly adapted function from sealABC::mssumstats to calculate sum stats, and return not only mean and SD but the complete dataset (ie all nresamp)
library(sealABC)
source("mssumstatsAP.R") #function links to other functions within sealABC package, so sealABC has to be loaded

set.seed(1122)
seal_stats <- mssumstatsAP(seal, by_pop = NULL, start_geno = 3, mratio = "loose", rarefaction = T, nresamp = 1000, nind = 10, nloc = 39)

seal_stats[[2]]$species <- "antarctic_fur_seal"

# # Save file
write.table(seal_stats[[2]], "data/all_sumstats_afs_1000_10indv.txt", quote = F, row.names = T)

