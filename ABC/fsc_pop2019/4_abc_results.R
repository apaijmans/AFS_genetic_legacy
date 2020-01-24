# save abc estimates to RData file
# replace filename neut with bot when working on bottleneck model estimates.

# Clear all previous info
rm(list=ls())


load("abc_estimates/abc_sims10000k_pop18_bot_ss_3.RData")
abc_est.ss <- abc_est
abc_posteriors.ss <- as.data.frame(abc_est.ss$unadj.values)
abc_posteriors.ss$pop <- "SSI"

load("abc_estimates/abc_sims10000k_pop18_bot_sg_3.RData")
abc_est.sg <- abc_est
abc_posteriors.sg <- as.data.frame(abc_est.sg$unadj.values)
abc_posteriors.sg$pop <- "SG"

load("abc_estimates/abc_sims10000k_pop18_bot_bi_3.RData")
abc_est.bi <- abc_est
abc_posteriors.bi <- as.data.frame(abc_est.bi$unadj.values)
abc_posteriors.bi$pop <- "BI"

load("abc_estimates/abc_sims10000k_pop18_bot_mar_3.RData")
abc_est.mar <- abc_est
abc_posteriors.mar <- as.data.frame(abc_est.mar$unadj.values)
abc_posteriors.mar$pop <- "MAR"

load("abc_estimates/abc_sims10000k_pop18_bot_ci_3.RData")
abc_est.ci <- abc_est
abc_posteriors.ci <- as.data.frame(abc_est.ci$unadj.values)
abc_posteriors.ci$pop <- "CI"

load("abc_estimates/abc_sims10000k_pop18_bot_ki_3.RData")
abc_est.ki <- abc_est
abc_posteriors.ki <- as.data.frame(abc_est.ki$unadj.values)
abc_posteriors.ki$pop <- "KI"

load("abc_estimates/abc_sims10000k_pop18_bot_hi_3.RData")
abc_est.hi <- abc_est
abc_posteriors.hi <- as.data.frame(abc_est.hi$unadj.values)
abc_posteriors.hi$pop <- "HI"

load("abc_estimates/abc_sims10000k_pop18_bot_mac_3.RData")
abc_est.mac <- abc_est
abc_posteriors.mac <- as.data.frame(abc_est.mac$unadj.values)
abc_posteriors.mac$pop <- "MAC"

abc_bot <- rbind(abc_posteriors.ss, abc_posteriors.sg, abc_posteriors.bi, abc_posteriors.mar, 
                 abc_posteriors.ci, abc_posteriors.ki, abc_posteriors.hi, abc_posteriors.mac)

abc_bot$pop <- as.factor(abc_bot$pop)

# Sort pops in right order for visual boxplot
abc_bot$pop <- ordered(abc_bot$pop, levels=c("SSI", "SG", "BI", "MAR", "CI", "KI", "HI", "MAC"))

save(abc_bot, file = "abc_estimates/abc_sims10000k_pop18_bot_ALL_3.RData")
