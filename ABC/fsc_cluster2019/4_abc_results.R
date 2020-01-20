# save abc estimates to RData file
# replace filename neut with bot when working on bottleneck model estimates.

# Clear all previous info
rm(list=ls())


load("abc_estimates/abc_sims10000k_cluster181_bot_ss_3.RData")
abc_est.ss <- abc_est
abc_posteriors.ss <- as.data.frame(abc_est.ss$unadj.values)
abc_posteriors.ss$pop <- "SSI"

load("abc_estimates/abc_sims10000k_cluster181_bot_sg_3.RData")
abc_est.sg <- abc_est
abc_posteriors.sg <- as.data.frame(abc_est.sg$unadj.values)
abc_posteriors.sg$pop <- "SG"

load("abc_estimates/abc_sims10000k_cluster181_bot_bi_3.RData")
abc_est.bi <- abc_est
abc_posteriors.bi <- as.data.frame(abc_est.bi$unadj.values)
abc_posteriors.bi$pop <- "BI"

load("abc_estimates/abc_sims10000k_cluster181_bot_mar.ci_3.RData")
abc_est.mar.ci <- abc_est
abc_posteriors.mar.ci <- as.data.frame(abc_est.mar.ci$unadj.values)
abc_posteriors.mar.ci$pop <- "MAR_CI"

load("abc_estimates/abc_sims10000k_cluster181_bot_ki.hi.mac_3.RData")
abc_est.ki.hi.mac <- abc_est
abc_posteriors.ki.hi.mac <- as.data.frame(abc_est.ki.hi.mac$unadj.values)
abc_posteriors.ki.hi.mac$pop <- "HI_KI_MAC"

abc_bot <- rbind(abc_posteriors.ss, abc_posteriors.sg, abc_posteriors.bi, abc_posteriors.mar.ci,
                  abc_posteriors.ki.hi.mac)

abc_bot$pop <- as.factor(abc_bot$pop)

# Sort pops in right order for visual boxplot
abc_bot$pop <- ordered(abc_bot$pop, levels=c("SSI", "SG", "BI", "MAR_CI", "HI_KI_MAC"))

save(abc_bot, file = "abc_estimates/abc_sims10000k_cluster181_bot_ALL_3.RData")

