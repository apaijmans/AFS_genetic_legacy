### ABC analysis part 1: Model selection and evaluation

## This script will take the simulations from script 1_simulate_diversity as input and
# (1) plot the different summary statistics as boxplots for all models
# (2) calculate the probabilities of the models
# (3) calculate the fit of the empirical data to the model
# (4) output all plots / txt files under plots and results

# Note: Script creates some folders and puts results in.

# files needed:
# (1) sumstats for emperical data: all_sumstats_full_all_gaz_5cluster.txt (generated in script a+b)
# (2) simulation results are saved as text file sims10000k_cluster181.txt

# Clear all previous info
#rm(list=ls())

# packages
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


######  preparation ######

# how many cores should be left free?
cores_not_to_use <- 15 # was 45

# How many summary statistics do we use?
# 3: "num_alleles_mean", "prop_low_afs_mean", "mratio_mean"
# 5: "num_alleles_mean", "prop_low_afs_mean", "mean_allele_range",  "mratio_mean",  "exp_het_mean")

num_sum_stats <- 3

# load genetic data #

# this has to be the exact name of the simulation test file
# and will be used later for filenames
sim_name <- paste0("sims10000k_cluster181") 

###### Load emperical data ###### 

# Read in data
all_sumstats_full <- read.table("all_sumstats_full_all_gaz_5cluster.txt", header = T, row.names=1, stringsAsFactors = F)

#########################################


#### select summary statistics ######

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

sumstats # check; right number of sum stats?

all_sumstats <- all_sumstats_full[sumstats]


####### run abc step 1 ########

### load simulations, stored in main folder ###
path_to_sims <- paste0(sim_name, ".txt")
sims <- fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

### subsetting and definition of model selection parameters

# parameter columns in simulation data.frame
param_start <- which(names(sims) == "sample_size")
param_end <- which(names(sims) == "range_constraint")
params <- c(param_start:param_end)

# create a character vector with models
models <- sims$model

# tolerance rate
tol <- 0.0005

# cross-validation replicates / number of replicates used to estimate the null distribution of the goodness-of-fit statistic
cv_rep <- 100 #max usually 100, 50 is faster

# method for model selection with approximate bayesian computation, see ?postpr
method <- 'rejection'

# extract names of all models
model_names <- names(table(models))

# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]



##### (1) first visual checks  -------------------------------------------------
# just calculate the goodness of fit for each model?

dir_check1 <- "model_evaluation/check1_sumstats/"
#if (!dir.exists(dir_check1)) dir.create(dir_check1)

# check whether sumstats are different across models
pdf(paste0(dir_check1, sim_name, "_", num_sum_stats, ".pdf"))
par(mfcol = c(2, 3), mar = c(4,4,1,1))
for (i in sumstats) {
  boxplot(sims[[i]] ~ models, main = i)
}
dev.off() 

# Make nicer version of this plot for supplementary materials (Fig S6), 
# see script Supp_conf_mat_plot_FigS6

### (2) can abc at all distinguish between the 2 models ? ----------------------
# leave-one-out cv and confusion matrix

dir_check2 <- "model_evaluation/check2_models/"
#if (!dir.exists(dir_check2)) dir.create(dir_check2)

# model cv, rejection method
cv.modsel <- cv4postpr(models, sims_stats, nval = cv_rep, tol = tol, method = method)
# summary
post_probs_summary <- summary(cv.modsel)
# model_missclassification probabilites
write_delim(x =  round(as.data.frame(post_probs_summary$probs),3), 
            path = paste0(dir_check2, sim_name, "_model_missclass_probs_", num_sum_stats, ".txt"))
# model_missclassification frequencies
write_delim(x =  as.data.frame(post_probs_summary$conf.matrix), 
            path = paste0(dir_check2, sim_name, "_model_missclass_freq_", num_sum_stats, ".txt"))

# plot
png(paste0(dir_check2, sim_name, "_confusion_mat", num_sum_stats, ".png"), width = 4, height = 4, units = "in", res = 300)
plot(cv.modsel, names.arg = model_names)
dev.off() 



### (3) model selection --------------------------------------------------------

dir_modselection <- "results/model_probs/"
#if (!dir.exists(dir_modselection)) dir.create(dir_modselection)

# for each genetic cluster
all_probs_ss <- abc::postpr(target = all_sumstats["ss",], index = models, sumstat = sims_stats, tol = tol, method = method)
all_probs_sg <- abc::postpr(target = all_sumstats["sg",], index = models, sumstat = sims_stats, tol = tol, method = method)
all_probs_bi <- abc::postpr(target = all_sumstats["bi",], index = models, sumstat = sims_stats, tol = tol, method = method)
all_probs_marci <- abc::postpr(target = all_sumstats["mar.ci",], index = models, sumstat = sims_stats, tol = tol, method = method)
all_probs_kihimac <- abc::postpr(target = all_sumstats["ki.hi.mac",], index = models, sumstat = sims_stats, tol = tol, method = method)

# if method is regression
out1 <- round(summary(all_probs_ss)$Prob, 3)
out2 <- round(summary(all_probs_sg)$Prob, 3)
out3 <- round(summary(all_probs_bi)$Prob, 3)
out4 <- round(summary(all_probs_marci)$Prob, 3)
out5 <- round(summary(all_probs_kihimac)$Prob, 3)

# # if method is mnlogistic
# out1 <- round(all_probs_ss$pred, 3)
# out2 <- round(all_probs_sg$pred, 3)
# out3 <- round(all_probs_bi$pred, 3)
# out4 <- round(all_probs_mar$pred, 3)
# out5 <- round(all_probs_ci$pred, 3)
# out6 <- round(all_probs_ki$pred, 3)
# out7 <- round(all_probs_hi$pred, 3)
# out8 <- round(all_probs_mac$pred, 3)

out <- rbind(out1, out2, out3, out4, out5)
row.names(out) <- c("ss", "sg", "bi", "mar.ci", "ki.hi.mac")

#stopCluster(cl)
write.table(out, file = paste0(dir_modselection, sim_name, "_model_selection_", num_sum_stats, ".txt"), row.names = T)



### (4) Does the preferred model provide a good fit to the data?----------------

dir_modeval <- "model_evaluation/check3_modeval/"
#if (!dir.exists(dir_modeval)) dir.create(dir_modeval)

# calculate all fits
all_fits_neut1 <- abc::gfit(target = all_sumstats["ss",], sumstat = sims_stats,
                            nb.replicate = cv_rep, tol = tol, subset = models == "neut")
all_fits_neut2 <- abc::gfit(target = all_sumstats["sg",], sumstat = sims_stats,
                            nb.replicate = cv_rep, tol = tol, subset = models == "neut")
all_fits_neut3 <- abc::gfit(target = all_sumstats["bi",], sumstat = sims_stats,
                            nb.replicate = cv_rep, tol = tol, subset = models == "neut")
all_fits_neut4 <- abc::gfit(target = all_sumstats["mar.ci",], sumstat = sims_stats,
                            nb.replicate = cv_rep, tol = tol, subset = models == "neut")
all_fits_neut5 <- abc::gfit(target = all_sumstats["ki.hi.mac",], sumstat = sims_stats,
                            nb.replicate = cv_rep, tol = tol, subset = models == "neut")

#summary(all_fits_bot)
#plot(all_fits_bot)

# pvalue (if not sig different it means model is very similar to real data)
p1 <- summary(all_fits_neut1)$pvalue
p2 <- summary(all_fits_neut2)$pvalue
p3 <- summary(all_fits_neut3)$pvalue
p4 <- summary(all_fits_neut4)$pvalue
p5 <- summary(all_fits_neut5)$pvalue


p_neut <- rbind(p1, p2, p3, p4, p5)
#row.names(p_neut) <- c("ss", "sg", "bi", "mar.ci", "ki.hi.mac")
colnames(p_neut) <- "p_neut"

# save the p_values
#write.table(p_neut, file = paste0(dir_modeval, sim_name, "_p_vals_fit_neut.txt"), row.names = T)

# calculate all fits
all_fits_bot1 <- abc::gfit(target = all_sumstats["ss",], sumstat = sims_stats,
                           nb.replicate = cv_rep, tol = tol, subset = models == "bot")
all_fits_bot2 <- abc::gfit(target = all_sumstats["sg",], sumstat = sims_stats,
                           nb.replicate = cv_rep, tol = tol, subset = models == "bot")
all_fits_bot3 <- abc::gfit(target = all_sumstats["bi",], sumstat = sims_stats,
                           nb.replicate = cv_rep, tol = tol, subset = models == "bot")
all_fits_bot4 <- abc::gfit(target = all_sumstats["mar.ci",], sumstat = sims_stats,
                           nb.replicate = cv_rep, tol = tol, subset = models == "bot")
all_fits_bot5 <- abc::gfit(target = all_sumstats["ki.hi.mac",], sumstat = sims_stats,
                           nb.replicate = cv_rep, tol = tol, subset = models == "bot")

#summary(all_fits_bot)
#plot(all_fits_bot)

# pvalue (if not sig different it means model is very similar to real data)
pb1 <- summary(all_fits_bot1)$pvalue
pb2 <- summary(all_fits_bot2)$pvalue
pb3 <- summary(all_fits_bot3)$pvalue
pb4 <- summary(all_fits_bot4)$pvalue
pb5 <- summary(all_fits_bot5)$pvalue

# combine into table
p_bot <- rbind(pb1, pb2, pb3, pb4, pb5)
colnames(p_bot) <- "p_bot"
all_p <- cbind(p_bot, p_neut)
row.names(all_p) <- c("ss", "sg", "bi", "mar.ci", "ki.hi.mac")

# save the p_values
write.table(all_p, file = paste0(dir_modeval, sim_name, "_p_vals_fit_", num_sum_stats, ".txt"), row.names = T)

