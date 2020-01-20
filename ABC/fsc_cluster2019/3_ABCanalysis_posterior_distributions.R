# ABC analysis part 2: Parameter estimation

# (1) definition of parameter values and wether to run a cross-validation
# on each parameter via pseudo-observed datasets taken from the prior distribution
# (2) calculates summary statistics with the sealABC package and transforms
# them into a clean data.frame
# (3) user selects summary statistics
# (4) simulations are read from a text file in split up into parameters
# and summary statistics
# (5) user chooses one of the models to do abc on (will be generalized to all models soon)
# (6) cross-validation with cv4abc
# (7) abc for all clusters/populations where a given model had the highest probability as
# saved under results/model_probs. abc can be done with different methods
# on all clusters/populations and paramaters and is saved under abc_estimates/

# Clear all previous info
#rm(list=ls())

# load packages 
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library(abc)
library(ggplot2)
library(readxl)
library(dplyr)
library(magrittr)
library(parallel)
library(readr)
# leave cores free
cores_not_to_use <- 15 # was 45

# (1) parameter definition ---------------------------------------------------------

# path to simulation
sim_name <- "sims10000k_cluster181"

#number of populations
num_pop <- 5

# How many summary statistics do we use?
# 3: "num_alleles_mean", "prop_low_afs_mean", "mratio_mean"
# 5: "num_alleles_mean", "prop_low_afs_mean", "mean_allele_range",  "mratio_mean",  "exp_het_mean")

num_sum_stats <- 3

# do a cross validation for the abc parameters? ?cv4abc
#calc_cv_for_abc <- TRUE
# if yes, which method?
method_cv <- "rejection"      # "rejection"    # "loclinear"
# how many pseudo-observed datasets should be evaluated? (value given in multiples of 5)
nval_cv <- 1000
# which tolerance level(s)?
tols_cv <- c(0.0005)


# tolerance level for abc
tol_abc <- 0.0005
## abc method choice, all three possible
all_methods <- c("rejection") # "ridge", "loclinear", "neuralnet"


# (2) load emperical data ----------------------------------------------------------

all_sumstats_full <- read.table("all_sumstats_full_all_gaz_5cluster.txt", header = T, row.names=1, stringsAsFactors = F)

# (3) select summary statistics for posteriors. ------------------------------------

ifelse(
  num_sum_stats == 5,
  sumstats <- c(
    "num_alleles_mean",
    "prop_low_afs_mean",
    "mean_allele_range",
    "mratio_mean",
    "exp_het_mean"
  ),
  ifelse(
    num_sum_stats == 3,
    sumstats <- c("num_alleles_mean",
                  "prop_low_afs_mean",
                  "mratio_mean"),
    print("Define number of summary statistics: 3 or 5")
  )
)

sumstats # check, right number of sum stats?


all_sumstats_full <- all_sumstats_full[sumstats]
#rownames(all_sumstats_full) <-


# (4) load simulations -------------------------------------------------------------

path_to_sims <- paste0(sim_name,".txt")
sims <- fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# create dataframes with parameters and summary statistics ---------------------

# parameter columns in simulation data.frame
param_start <- which(names(sims) == "sample_size")
param_end <- which(names(sims) == "range_constraint")
params <- c(param_start:param_end)
# create a character vector with models
models <- sims$model
# extract names of all models
model_names <- names(table(models))
# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]

# (5) models
best_mod <- "bot" # or "neut", ... depending on which is the best model

mod <- best_mod

# subset sims_stats and sims_param for the given model
stat_mod <- subset(sims_stats, subset = models == mod)
par_mod <- subset(sims_param, subset = models == mod)

# (6) cv in parallel ------------------------------------------------------------------

# before inference, we see whether a parameter can be estimated at all

cv_nbot <- function(iter, par_mod, stat_mod, method, pars, tols = tol){ 
  cv_res_rej <- cv4abc(data.frame(par_mod)[pars], stat_mod, nval = 5,
                       tols = tols, method = method)
}
pars <- c("pop_size", "nbot", "nhist", "tbotend", "tbotstart", "mut_rate", "gsm_param")

#try with lower nval for testing
#nval_cv <- 100

cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - cores_not_to_use))
clusterEvalQ(cl, c(library("abc")))
all_cv_res <- parLapply(cl, 1:(nval_cv/5), cv_nbot, par_mod, stat_mod, method_cv, pars, tols_cv)
stopCluster(cl)

# # first one
# all_cv <- all_cv_res[[1]]

#cmlist1=lapply(seq(1,length(clist),by=2),function(x)do.call("rbind", clist[x:(x+1)]))
### Merging metalist of cv4abc runs
# relevant outputs: cvsamples, true, estim, seed

#first make final list:

finlist <- all_cv_res[[1]]

#unlist first "layer" of meta list:
bla <- unlist(all_cv_res, recursive=F)

# concatenate all cv sample vector of each list, then overwrite into final list
finlist[[2]] <- unlist(bla[seq(2,length(bla), 7)], use.names = F)

# concatenate all true dataframes of each list, then overwrite into final list
finlist[[4]] <- bind_rows(bla[seq(4,length(bla), 7)])

# concatenate all estim matrixes of each list, then overwrite into final list
finlist[[5]][[1]] <- do.call(rbind, unlist(bla[seq(5,length(bla), 7)], recursive = F))

# concatenate all seed vector of each list, then overwrite into final list
finlist[[7]] <- unlist(bla[seq(7,length(bla), 7)], use.names = F)

summary(finlist) # seems to work!

# ### Debugging ###
# pars <- c("pop_size", "nbot", "nhist", "tbotend", "tbotstart", "mut_rate", "gsm_param")
# all_cv <- cv4abc(data.frame(par_mod)[pars], stat_mod, nval = 100,
#                      tols = tols_cv, method = method_cv)

all_cv <- finlist

pred_error <- summary(all_cv) # table for prediction error
write.table(x =  round(as.data.frame(pred_error[1,]), 3), 
            file = paste0("model_evaluation/check4_params/", sim_name, "pred_error_bot_", num_sum_stats, ".txt"),
            quote = F)
# plot(all_cv)

# extract values for binding

# true values
true_vals <- data.frame(do.call(rbind, lapply(all_cv_res, function(x) x$true)))
cv_samples <- as.numeric(unlist(data.frame(do.call(rbind, lapply(all_cv_res, function(x) data.frame(x$cvsamples))))))

# estimated values
estim_vals <- list()
for (i in 1:length(all_cv$estim)) {
  estim_vals[[i]] <- data.frame(do.call(rbind, lapply(all_cv_res, function(x) x$estim[[i]])))
}
names(estim_vals) <- names(all_cv$estim)

# cv_samples
all_cv$cvsamples <- cv_samples
all_cv$true <- true_vals
all_cv$estim <- estim_vals

out <- paste0("model_evaluation/check4_params/cv_param_rej_", sim_name,"_",mod,"_", num_sum_stats, ".RData")
# write.table(cv_res, file = out, row.names = FALSE)
save(all_cv, file = out)



# (7) run the actual abc analysis --------------------------------------------------

dir.create("abc_estimates/")

# in code below abc is only done on best model (ie bottleneck or neutral)

for (i in 1:num_pop) {
  
  ## load model probabilities
  model_probs <- read.table(paste0("results/model_probs/", 
                                   sim_name, "_model_selection_3.txt"), header = TRUE, row.names = 1)
  
  model_bot <- model_probs$bot > 0.5
  #model_neut <- model_probs$bot <= 0.5
  
  # extract pop names
  all_pops <- rownames(all_sumstats_full)
  
  if (model_bot[i] == TRUE) {
    
    
    mod <- "bot"  
    
    # subset sims_stats and sims_param for the given model
    stat_mod <- subset(sims_stats, subset = models == mod)
    par_mod <- subset(sims_param, subset = models == mod)
    
    abc_est <- abc(target = all_sumstats_full[i,], param = par_mod, 
                   sumstat = stat_mod, tol = tol_abc, method = all_methods) #had to change tol_abc to 0.5 otherwise error?
    
    save(abc_est, file = paste0("abc_estimates/abc_", sim_name,"_",mod,"_",all_pops[i], "_", num_sum_stats, ".RData"))
    
    abc_est <- NULL
    
    
  } else if (mmodel_bot[i] == F) {
    mod <- "neut"  
    
    # subset sims_stats and sims_param for the given model
    stat_mod <- subset(sims_stats, subset = models == mod)
    par_mod <- subset(sims_param, subset = models == mod)
    
    abc_est <- abc(target = all_sumstats_full[i,], param = par_mod, 
                   sumstat = stat_mod, tol = tol_abc, method = all_methods) #had to change tol_abc to 0.5 otherwise error?
    
    save(abc_est, file = paste0("abc_estimates/abc_", sim_name,"_",mod,"_",all_pops[i], "_", num_sum_stats, ".RData"))
  }
}



