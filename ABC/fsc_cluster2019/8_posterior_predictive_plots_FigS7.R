# evaluate post.predictive checks
# Fig S7

library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
source("martin.R")

# simulated based on posteriors
# created in 7_posterior_predictive_checks
all_checks <- read_delim("model_evaluation/check5_postpred/sims_10000k_cluster181_bot_post_pred_checks1.txt", delim = " ")

head(all_checks)

#remove accidentally double pop column
all_checks$pop_1 <- NULL
all_checks$pop <- gsub("HI_KI_MAC", "KI_HI_MAC" , all_checks$pop)

# compare these summary statistics
sumstats <- c("num_alleles_mean", 
              #"exp_het_mean", 
              "mratio_mean", 
              "prop_low_afs_mean")
              #"mean_allele_range")

# empirical summary stats
all_stats <- read.table("all_sumstats_full_all_gaz_5cluster.txt") 
all_stats$pop <- toupper(row.names(all_stats))
all_stats$pop <- gsub("\\.", "_" , all_stats$pop)
all_stats$pop <- gsub("SS", "SSI" , all_stats$pop)


# long format for plotting
all_checks_long <- all_checks %>% 
  dplyr::select(pop, sumstats) %>% 
  gather(sumstat, value, -pop) 

# observed sumstats long format
all_sumstats_full_long <- all_stats %>% 
  dplyr::select(pop, sumstats) %>%
  gather(sumstat, value, -pop)

# Setting order of genetic clusters, color and names for plots
all_sumstats_full_long$pop <- ordered(all_sumstats_full_long$pop, levels=c("SSI", "SG", "BI", "MAR_CI", "KI_HI_MAC"))
all_checks_long$pop <- ordered(all_checks_long$pop, levels=c("SSI", "SG", "BI", "MAR_CI", "KI_HI_MAC"))

# Define colors
col <- c( "#E7298A","#6A3D9A","#B2DF8A", "#CFC855", "#FFA500")

# Define sumstat names
sumstat_names <- c(
  exp_het_mean = "Expected\nhetetozygosity",
  mean_allele_range = "Allelic range",
  mratio_mean = "M-ratio",
  num_alleles_mean = "Allelic richness",
  prop_low_afs_mean  = "Prop. of low\nfrequency alleles"
)

# Define genetic cluster names
full_pop <- c(SSI = "South\nShetlands", 
              SG = "South\nGeorgia", 
              BI = "Bouvetøya", 
              MAR_CI = "Marion Island &\nCrozet Islands", 
              KI_HI_MAC = "Kerguelen Islands,\nHeard Island &\nMacquarie Island")

# Make plot
p <- ggplot(all_checks_long, aes(value, fill = pop)) +
  geom_histogram() +
  geom_vline(aes(xintercept = value), all_sumstats_full_long) +
  facet_grid(pop ~ sumstat, scales = "free", labeller = labeller(
    sumstat = sumstat_names, pop=full_pop
  )) +
  theme_martin(base_family = "Arial", highlight_family = "Arial") +
  scale_fill_manual(values = col) +
  guides(fill=FALSE)+
  theme(strip.text.y = element_text(angle = 0),
        panel.spacing = unit(1, "lines")) 

p

ggsave(filename = "figs/post_pred_checks.jpg", plot = p, width = 7, height = 7)
