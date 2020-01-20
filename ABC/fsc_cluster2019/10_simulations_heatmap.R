
# Simulations to calculate diversity loss using different Nebot and duration
# Input for Fig 4d, figure itself is made in script 11_plotting_heatmap

library(strataG)
library(dplyr)
# library(truncnorm)
library(parallel)
options(scipen = 999) # disables scientific notation, ie 1e+10

# number of cores to use on system
num_cores <- 40

# number of coalescent simulations
num_sim <- 1000

# will be used later for filenames
sim_name <- paste0("sims_heatmap")

# create data.frame with all parameter values ---------------

# effective current pop size, bottleneck size and historical pop size
# pop_size = nhist = 10,000
# nbot from 25-600 in steps of 25, plus for the non-bottleneck model nbot=ncur

nbot <- c(rep(10000, each = num_sim), rep(seq(25, 600, by=25), each = num_sim))
pop_size <- rep(10000, times=length(nbot))
nhist <- rep(10000, times=length(nbot))

all_N <- as.data.frame(cbind(pop_size, nbot, nhist))

# calculate popsizes relative to current effective popsize
all_N <- mutate(all_N, nbot_prop = nbot / pop_size)
all_N <- mutate(all_N, nhist_bot_prop = nhist / nbot)
all_N <- mutate(all_N, nhist_neut_prop = nhist / pop_size)

# Duration bottleneck 1 gen, then going up
all_N <- mutate(all_N, tbotend = 10)
all_N2<- mutate(all_N, tbotstart = 11)

for (i in 2:20){
  
  all_N1<- mutate(all_N, tbotstart = 10+i)
  all_N2 <- rbind(all_N2, all_N1)
}


# mutation model
# mutation rate fixed
all_N2<- mutate(all_N2, mut_rate = 1e-4)
# parameter of the geometric distribution: decides about the proportion 
# of multistep mutations, fixed
all_N2<- mutate(all_N2, gsm_param = 0.2)

all_N2<- mutate(all_N2, range_constraint = 30)

# sample size
sample_size <- rep(181, nrow(all_N2)) # equal to smallest pop size, 181 for clusters
# number of loci
num_loci <- rep(39, nrow(all_N2))

all_params <- data.frame(sample_size, num_loci, all_N2, param_num = 1:nrow(all_N2))

param_set <- all_params


# Function to run simulations using the parameter set we just created
run_sims <- function(param_set){
  
  lab <- as.character(param_set[["param_num"]])
  pop_info <- strataG::fscPopInfo(pop.size = param_set[["pop_size"]], sample.size = param_set[["sample_size"]])
  mig_rates <- matrix(0)
  
  # Define historical events in which populations diverged x generations in past
  
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
  
  # num_alleles mean and sd
  num_alleles <- strataG::numAlleles(sim_msats)
  num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
  
  # mratio mean and sd
  mratio <- strataG::mRatio(sim_msats, by.strata = FALSE, rpt.size = 1)
  mratio_mean <- mean(mratio, na.rm = TRUE)
  
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
  
  out <- data.frame(
    num_alleles_mean, mratio_mean, prop_low_afs_mean
  )
  
  
}


#all_sims <- parApply(cl, all_params, 1, run_sims)
#sims_df <- as.data.frame(data.table::rbindlist(all_sims))

# Actually running the function on cluster with x cores
cl <- makeCluster(getOption("cl.cores", num_cores))
clusterEvalQ(cl, c(library("strataG")))

all_sims <- parApply(cl, all_params, 1, run_sims)
sims_df <- as.data.frame(data.table::rbindlist(all_sims))

stopCluster(cl)

# reshape data to get a clean data.frame
sims <- cbind(sims_df, all_params)

sims$model <- rep(c(rep("neut", num_sim), rep("bot", nrow(all_N)-num_sim)), times = 20) # check whether this works or goes wrong

# save simulations in a txt file
write.table(sims, file = paste0(sim_name, ".txt"), row.names = FALSE)
