
### Plot species Ar

land_species <- read.table("data/all_sumstats_landseals_1000_10indv.txt", header=T)

library(tidyr)
library(dplyr)

# Keep only otariids
otariids <- land_species %>%
  filter(species != "grey_seal_orkneys") %>%
  filter(species != "harbour_seal_waddensee") %>%
  filter(species != "hawaiian_monk_seal") %>%
  filter(species != "mediterranean_monk_seal") %>%
  filter(species != "nes") %>%
  filter(species != "ses")


# Add our own data for AFS to the landspecies df
afs <- read.table("data/all_sumstats_afs_1000_10indv.txt", header=T)
otariids <- rbind(otariids, afs)

otariid_names <- c(
  "antarctic_fur_seal" = "Antarctic Fur Seal",
  "australian_fur_seal" = "Australian Fur Seal",
  "california_sea_lion" = "California Sea Lion",
  "galapagos_fur_seal" = "Galapagos Fur Seal",
  "galapagos_sea_lion" = "Galapagos Sea Lion",
  # "grey_seal_orkneys" = "Grey Seal",
  "guadalupe_fur_seal"= "Guadalupe Fur Seal",
  # "harbour_seal_waddensee" = "Harbour Seal",
  # "hawaiian_monk_seal" = "Hawaiian Monk Seal",
  # "mediterranean_monk_seal" = "Mediterranean Monk Seal",
  "new_zealand_fur_seal" = "New Zealand Fur Seal",
  "new_zealand_sea_lion"= "New Zealand Sea Lion",
  # "nes" = "Northern Elephant Seal",
  "northern_fur_seal" = "Northern Fur Seal",
  "south_american_fur_seal" = "South American Fur Seal",
  "south_american_sea_lion" = "South American Sea Lion",
  # "ses" = "Southern Elephant Seal",
  "stellers_sea_lion" = "Steller Sea Lion",
  "subantarctic_fur_seal" = "Subantarctic Fur Seal")

library(ggbeeswarm)
library(grid)

estimate_mode <- function(s) {
  d <- density(s, adjust = 2.5)
  d$x[which.max(d$y)]
}

source("~/fsc_pop2018/martin.R")
#source("martin.R")

#library(dplyr)
#library(extrafont)
p <- otariids %>%
  ggplot(aes(num_alleles_mean, 
             x=reorder(species, num_alleles_mean, median),
             color=factor(ifelse(species=="antarctic_fur_seal","Highlighted","Normal")))) + #, colour=species)) +
  geom_quasirandom(alpha = 0.07, size = 2, width = 0.37, bandwidth = 2.5) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "white", alpha = 0.5, size = 0.2) +
  stat_summary(fun.y = "estimate_mode", colour = "black", geom = "point", size = 2, shape = 21, fill = "grey") +
  scale_x_discrete(labels = otariid_names) +
  scale_y_continuous(name = expression(italic(A)[r])) +
  guides(color=FALSE)+ #remove legend (here use color, not fill, as we used color to play with colours)
  xlab("") +
  scale_color_manual(name="species", values=c("#FAA460","grey50")) + #pink #E7298A
  #scale_x_discrete(name="", limits = rev(levels(allstats_long$pop)), labels=rev(full_pop)) + #flip order of y axis, so that on top is SS, then SG, ...etc.
  theme_martin(base_family = "Arial", highlight_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x=element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))
p 

ggsave(filename = "Ar_otariids_beeswarm_AFS_median.jpg", plot = p, width = 9, height = 5)
