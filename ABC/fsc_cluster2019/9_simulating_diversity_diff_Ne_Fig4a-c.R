
# Here we simulate diversity loss using different Ne (1000, 10000 and 50000)
# Then make Fig 4 a, b & c

library(strataG)

## Not run: 
# Set fastsimcoal parameters
# Population information: 1 pop with Ne = 10,000, drawing 181 samples.
Ne_cur <- 10000
Ne_bot <- 200

sim_diversity <- function(Ne_cur, Ne_bot, start_bot, end_bot) {
  
  pop_info <- fscPopInfo(pop.size = Ne_cur, sample.size = 181)
  
  if (is.null(Ne_bot)) Ne_bot <- Ne_cur
  
  # Migration rates among the 3 populations
  # mig.rates <- matrix(c(0, 0.5, 0.005, 0.5, 0, 0.0005, 0.005, 0.0005, 0), ncol = 3)
  
  # Define historical events in which populations diverged 2000 generations in past
  hist_ev <- fscHistEv(
    num.gen = c(start_bot, end_bot), source.deme = c(0, 0), prop.migrants = 1,
    sink.deme = c(0, 0), new.sink.size = c(Ne_bot / Ne_cur, Ne_cur / Ne_bot)
  )
  
  # Define 39 microsatellite loci, with random mutation rates
  msat_params <- fscLocusParams(
    locus.type = "msat", num.loci = 39, 
    mut.rate = 1e-4, gsm.param = 0.2, range.constraint = 30, ploidy = 2
  )
  
  # Run simulation and display locus summary
  sim_msats <- fastsimcoal(pop_info, msat_params, hist.ev = hist_ev, exec = "/home/anneke/bin/fsc26_linux64/fsc26")
  
  # num_alleles, allel_richness, prop_unique_alleles, expt_het, obs_het
  # mean and sd
  num_alleles <- strataG::numAlleles(sim_msats)
  num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
}


library(purrr)

Ne_currs <- rep(c(1000, 10000, 50000), each = 100000) #100.000 simulations of 1000, 10000 and 50000 each (overall: 300) 

with_bot <- purrr::map_dbl(Ne_currs, sim_diversity, Ne_bot = 200, start_bot = 20, end_bot = 30)
without_bot <- purrr::map_dbl(Ne_currs, sim_diversity, Ne_bot = NULL, start_bot = 20, end_bot = 30)
## End(Not run)


#library(tidyverse)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)
all_divs <- tibble(div = c(with_bot, without_bot), mod = rep(c("bot", "non_bot"), each = length(Ne_currs)), Ne = rep(Ne_currs, 2))

# Save all_divs so we can create figures without having to re-run simulations
save(all_divs, file = "simulating_divs/Nbot200_100000sims.RData")

# Loading data to create figures
load("~/fsc_cluster2019/simulating_divs/Nbot200_100000sims.RData")

all_divs2 <- all_divs %>%
  group_by(Ne, mod) %>%
  summarise(mean_div = mean(div))

mean_divs <- all_divs %>%
  group_by(mod, Ne) %>%
  summarize(div = mean(div)) %>%
  spread(mod, div) %>%
  mutate(div_lost_perc = 100* (1- (bot / non_bot)))


# Adding text for plotting to all_divs
all_divs$Ne <- as.factor(all_divs$Ne)

levels(all_divs$Ne)

all_divs$Ne <- plyr::revalue(all_divs$Ne, c(
  `1000` = expression(paste("(a) ",italic(N)[e], "hist = 1,000: 6.21% lost")), # for break line: expression(atop(paste("(a) ",italic(N)[e], "hist = 1,000:"), "6.21% lost"))
  `10000` = expression(paste("(b) ",italic(N)[e], "hist = 10,000: 13.93% lost")),
  `50000` = expression(paste("(c) ",italic(N)[e], "hist = 50,000: 17.75% lost"))))

levels(all_divs$Ne)

# Adding text for plotting to all_divs2
all_divs2$Ne <- as.factor(all_divs2$Ne)

levels(all_divs2$Ne) 

all_divs2$Ne <- plyr::revalue(all_divs2$Ne, c(
  `1000` = expression(paste("(a) ",italic(N)[e], "hist = 1,000: 6.21% lost")),
  `10000` = expression(paste("(b) ",italic(N)[e], "hist = 10,000: 13.93% lost")),
  `50000` = expression(paste("(c) ",italic(N)[e], "hist = 50,000: 17.75% lost"))))

levels(all_divs2$Ne) 

scientific_10 <- function(x) {
  parse(text=x/1000)
}

source("martin.R")

p4 <- ggplot(all_divs, aes(div, fill = mod)) +
  geom_histogram(position = "identity", alpha = 1, bins = 30) +
  scale_fill_manual("Model", labels = c("Bottleneck", "No Bottleneck"), values=c("#FAA460","grey45")) + # "#E876B1" pink
  ylab(expression(Number~of~simulated~datasets~x~10^3)) +
  xlab("Mean number of alleles per locus") +
  geom_vline(data = all_divs2, aes(xintercept = mean_div, col = mod, linetype = mod)) +
  scale_color_manual("Model", labels = c("Bottleneck", "No Bottleneck"), values = c("black", "black")) +
  scale_linetype_manual("Model", labels = c("Bottleneck", "No Bottleneck"), values = c("solid", "dashed")) +
  facet_wrap(Ne ~ ., scales = "free", labeller = labeller(Ne= label_parsed)) +
  #ggtitle(label = "Diversity lost in bottleneck (Ne of 200 for 10 generations) for three different
  #population sizes (Ne = current = historical)") +
  theme_martin(base_family = "Arial", highlight_family = "Arial")+
  theme(strip.text.x = element_text(size = 10, hjust = 0),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_y_continuous(label=scientific_10)

p4

ggsave( "figs/diversity_loss_200_overlap_3ne.jpg", p4, width = 10, height = 3)
