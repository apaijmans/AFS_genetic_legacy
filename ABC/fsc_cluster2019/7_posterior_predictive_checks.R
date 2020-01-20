# Posterior predictive checks
# Simulations
# input for Fig S7, fig itself is created in 8_posterior_predictive_plots_FigS7

# Start of the analysis: coalescent simulations. Uses strataG functions as interface.
# See strataG manual on how to install fastsimcoal26, which is needed for 
# the following script.

library(strataG)
library(dplyr)
library(parallel)

# posterior predictive checks
library(tidyr)
library(dplyr)
library(readr)

# load abc posterior data
load("abc_estimates/abc_sims10000k_cluster181_bot_ALL_3.RData") 
abc_bot <- abc_bot%>% 
  mutate(post_mod = "bot")

# model selection
mod_select <- read_delim("results/model_probs/sims10000k_cluster181_model_selection_3.txt",
                         delim = " ", col_names = c("pop", "bot", "neut"), skip = 1) %>%
  mutate(mod = ifelse(bot > 0.5, "bot", "neut")) %>%
  select(pop, mod)

# all pops  are supporting bot model

# number of simulations for posterior predictive check
num_sim <- 1000


abc_pars <- abc_bot %>% 
  group_by(pop) %>% 
  sample_n(num_sim)  #%>%


all_params <- abc_pars %>% 
  mutate(post_mod = ifelse(post_mod == "bot", 1, 0)) %>%
  mutate(param_num = 1:num_sim)  %>%
  ungroup()
str(all_params)

all_params <- as.data.frame(all_params)

# Change cluster/population info to numeric for the code to work
all_params$pop <- gsub("SSI", 1 , all_params$pop)
all_params$pop <- gsub("SG", 2 , all_params$pop)
all_params$pop <- gsub("BI", 3 , all_params$pop)
all_params$pop <- gsub("MAR_CI", 4 , all_params$pop)
all_params$pop <- gsub("HI_KI_MAC", 5 , all_params$pop)

all_params$pop <- as.numeric(all_params$pop)

#param_set <- all_params[1,]


# Run simulations

run_sims <- function(param_set){
  
  model <- as.character(param_set[["post_mod"]])
  lab <- as.character(param_set[["param_num"]])
  pop_info <- strataG::fscPopInfo(pop.size = param_set[["pop_size"]], sample.size = param_set[["sample_size"]])
  mig_rates <- matrix(0)
  
  hist_ev <- strataG::fscHistEv(
    num.gen = c(param_set[["tbotend"]], param_set[["tbotstart"]]), source.deme = c(0, 0),
    sink.deme = c(0, 0), new.sink.size = c(param_set[["nbot_prop"]], param_set[["nhist_bot_prop"]])
  )
  
  msat_params <- strataG::fscLocusParams(
    locus.type = "msat", num.loci = param_set[["num_loci"]], 
    mut.rate = param_set[["mut_rate"]], gsm.param = param_set[["gsm_param"]], 
    range.constraint = param_set[["range_constraint"]], ploidy = 2
  )
  
  sim_msats <- strataG::fastsimcoal(pop.info = pop_info, locus.params = msat_params, 
                                    hist.ev = hist_ev, exec = "/home/anneke/bin/fsc26_linux64/fsc26", label = lab) # , 
  
  
  # calculate summary statistics
  
  # num_alleles, allel_richness, prop_unique_alleles, expt_het, obs_het
  # mean and sd
  num_alleles <- strataG::numAlleles(sim_msats)
  num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
  num_alleles_sd <- sd(num_alleles, na.rm = TRUE)
  # exp_het
  exp_het <- strataG::exptdHet(sim_msats)
  exp_het_mean <- mean(exp_het, na.rm = TRUE)
  exp_het_sd <- sd(exp_het, na.rm = TRUE)
  # obs_het
  obs_het <- strataG::obsvdHet(sim_msats)
  obs_het_mean <- mean(obs_het, na.rm = TRUE)
  obs_het_sd <- sd(obs_het, na.rm = TRUE)
  # mratio mean and sd
  mratio <- strataG::mRatio(sim_msats, by.strata = FALSE, rpt.size = 1)
  mratio_mean <- mean(mratio, na.rm = TRUE)
  mratio_sd <- stats::sd(mratio, na.rm = TRUE)
  # allele frequencies
  afs <- strataG::alleleFreqs(sim_msats)
  # prop low frequency alleles
  prop_low_af <- function(afs){
    # low_afs <- (afs[, "freq"] / sum(afs[, "freq"])) < 0.05
    low_afs <- afs[, "prop"] <= 0.05
    prop_low <- sum(low_afs) / length(low_afs)
  }
  # and mean/sd for all
  prop_low_afs <- unlist(lapply(afs, prop_low_af))
  prop_low_afs_mean <- mean(prop_low_afs, na.rm = TRUE)
  prop_low_afs_sd <- stats::sd(prop_low_afs, na.rm = TRUE)
  # allele range
  allele_range <- unlist(lapply(afs, function(x) diff(range(as.numeric(row.names(x))))))
  mean_allele_range <- mean(allele_range, na.rm = TRUE)
  sd_allele_range <- sd(allele_range, na.rm = TRUE)
  
  # allele size variance and kurtosis
  # create vector of all alleles per locus
  all_alleles <- function(afs_element){
    alleles <- as.numeric(rep(row.names(afs_element), as.numeric(afs_element[, "freq"])))
    size_sd <- stats::sd(alleles)
    size_kurtosis <- moments::kurtosis(alleles, na.rm = TRUE)
    out <- data.frame(size_sd = size_sd, size_kurtosis = size_kurtosis)
  }
  all_allele_size_ss <- do.call(rbind, lapply(afs, all_alleles))
  
  mean_allele_size_sd <- mean(all_allele_size_ss$size_sd, na.rm = TRUE)
  sd_allele_size_sd <- sd(all_allele_size_ss$size_sd, na.rm = TRUE)
  
  mean_allele_size_kurtosis <- mean(all_allele_size_ss$size_kurtosis, na.rm = TRUE)
  sd_allele_size_kurtosis <- sd(all_allele_size_ss$size_kurtosis, na.rm = TRUE)
  
  out <- data.frame(
    num_alleles_mean, num_alleles_sd,
    exp_het_mean, exp_het_sd,
    obs_het_mean, obs_het_sd,
    mean_allele_size_sd, sd_allele_size_sd,
    mean_allele_size_kurtosis, sd_allele_size_kurtosis,
    mean_allele_range, sd_allele_range,
    mratio_mean, mratio_sd,
    prop_low_afs_mean, prop_low_afs_sd
  )
}


# Run function on cluster with 40 cores
# cl <- makeCluster(getOption("cl.cores", 10))
# clusterEvalQ(cl, c(library("strataG")))
# 
# sims_all <- parApply(cl, all_params, 1, run_sims)
# sims_df <- as.data.frame(data.table::rbindlist(sims_all))
# 
# stopCluster(cl)

sims_all <- apply(all_params, 1, run_sims)
sims_df <- as.data.frame(data.table::rbindlist(sims_all))

# reshape data to get a clean data.frame
pops <- abc_pars %>%
  select(pop)

sims <- cbind(as.data.frame(pops), sims_df, all_params)

# save simulations in a txt file
write.table(sims, file = "model_evaluation/check5_postpred/sims_10000k_cluster181_bot_post_pred_checks1.txt", row.names = FALSE)
