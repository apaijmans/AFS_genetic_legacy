
# Fig 3

library(tidyr)
library(dplyr)

# Read in data
all_sumstats_full <- read.table("all_sumstats_5cluster_full_1000_10indv.txt", header = T, row.names=1, stringsAsFactors = F)
all_sumstats_full$pop <- gsub("KI_HI_MAC", "HI_KI_MAC", all_sumstats_full$pop, fixed = TRUE)

# Select summary statistics that we want to plot + cluster info ("pop")
sumstats <- c("num_alleles_mean", 
              "obs_het_mean",
              "mratio_mean",
              "pop")  

all_sumstats <- all_sumstats_full[sumstats]


# long format for plotting
sumstats_long <- all_sumstats %>%
  gather("num_alleles_mean", 
         "obs_het_mean",
         "mratio_mean", 
         key = "stat", value = "values")

# add Nbot to long data for plotting all together in one plot
load("abc_estimates/abc_sims10000k_cluster181_bot_ALL_3.RData")

nbot_long <- abc_bot[c("nbot", "pop")] %>%
  gather("nbot", key = "stat", value = "values")

# add proportion of loci in het excess to long data for plotting all together in one plot
prophet <- read.table("BOTTLENECK_1-1000_cluster_long.txt", header = T, stringsAsFactors = F)

allstats_long <- rbind(sumstats_long, nbot_long, prophet)

allstats_long$pop <- as.factor(allstats_long$pop)
allstats_long$stat <- as.factor(allstats_long$stat)

allstats_long$pop <- ordered(allstats_long$pop, levels=c("SSI", "SG", "BI", "MAR_CI", "HI_KI_MAC"))
allstats_long$stat <- ordered(allstats_long$stat, levels=c("num_alleles_mean", "obs_het_mean", "mratio_mean", "prop_het.ex", "nbot"))

# Define colors
col <- c( "#E7298A","#6A3D9A","#B2DF8A", "#CFC855", "#FFA500")

library(ggbeeswarm)
library(grid)

estimate_mode <- function(s) {
  d <- density(s, adjust = 2.5)
  d$x[which.max(d$y)]
}

source("martin.R")
full_pop <- c("South\nShetlands", 
              "South\nGeorgia", 
              "Bouvetøya", 
              "Marion Island &\nCrozet Islands", 
              "Kerguelen Islands,\nHeard Island &\nMacquarie Island")

levels(allstats_long$stat)

allstats_long$stat <- plyr::revalue(allstats_long$stat, c(
  "num_alleles_mean" = expression(paste("(a) ",italic(A)[r])),
  "obs_het_mean" = "Hobs", #expression(paste("(b) ", italic(H)[obs])),
  "prop_het.ex" = expression(paste("(b) ", italic(prop)[het-exc])),
  "mratio_mean" = "M-ratio",
  "nbot" = expression(paste("(c) ", italic(N)[e], "bot"))))

levels(allstats_long$stat) # Check if levels assigned correctly


p <- allstats_long %>%
  filter(stat!="M-ratio") %>%
  filter(stat!="Hobs") %>%
  ggplot(aes(values, x=pop, colour=pop)) +
  geom_quasirandom(alpha = 0.07, size = 2, width = 0.37, bandwidth = 2.5) +
  
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "white", alpha = 0.5, size = 0.2) +
  stat_summary(fun.y = "estimate_mode", colour = "black", geom = "point", size = 2, shape = 21, fill = "grey") +
  
  scale_colour_manual(values=col)+ # only works if you add colour= inside the aes
  
  ylab("Value") +
  #scale_x_continuous(limits = c(5, 8.5), breaks = c(5, 5.5, 6, 6.5, 7, 7.5,8,8.5), labels = c(5, 5.5, 6, 6.5, 7, 7.5,8,8.5)) +
  guides(colour=FALSE)+ #remove legend (here use color not fill, as we used color to play with colours)
  xlab("") +
  coord_flip() +
  scale_x_discrete(name="", limits = rev(levels(allstats_long$pop)), labels=rev(full_pop)) + #flip order of y axis, so that on top is SS, then SG, ...etc.
  facet_grid(. ~ stat, scales = "free", labeller = labeller(stat= label_parsed)) +
  theme_martin(base_family = "Arial", highlight_family = "Arial")+
  theme(axis.text.y= element_text(size = 9.5),
        axis.text.x= element_text(size = 7.5), #changes size xaxis and yaxis so pops, and values. But not top labels or "Value"
        panel.spacing.x=unit(1, "lines")) # increases spacing between panels

p

ggplot2::ggsave(filename = "figs/Ar_Hobs_prophetexc_and_nbot_cluster.jpg", p, width = 8, height = 6)

