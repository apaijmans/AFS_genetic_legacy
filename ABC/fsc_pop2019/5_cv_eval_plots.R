
# cv_eval_plots

# Cross-validation plots for all parameters

library(ggplot2)
library(tidyr)
library(dplyr)
#source("martin.R")

load("model_evaluation/check4_params/cv_param_rej_sims10000k_pop18_bot_3.RData")

# long format for plotting
all_cv_est <- gather(data.frame(all_cv$estim[[1]]), key = "parameter", value = "estimate")
all_cv_true <- gather(data.frame(all_cv$true), key = "parameter", value = "true")

all_cv_est$parameter <- NULL # order of parameters same for both all_cv_est and all_cv_true so can remove this column and then cbind

cv_params <- cbind(all_cv_true, all_cv_est)

cv_params$parameter <- as.factor(cv_params$parameter)

levels(cv_params$parameter) 

cv_params$parameter <- plyr::revalue(cv_params$parameter, c(
  "gsm_param" = expression(GSM[par]),
  "mut_rate" = expression(mu),
  "nbot" = expression(paste(italic(N)[e], "bot")),
  "nhist" = expression(paste(italic(N)[e], "hist")),
  "pop_size" = expression(paste(italic(N)[e], "cur")),
  "tbotend" = expression(paste(italic(t)[bot], "end")),
  "tbotstart" = expression(paste(italic(t)[bot], "start"))))

levels(cv_params$parameter) # Check if levels assigned correctly


p_all <- ggplot(cv_params, aes(true, estimate)) +
  geom_point(size = 2.3, alpha = 0.3) +
  #theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
  geom_abline(intercept = 0, slope=1, color='blue') +
  #geom_line(data = data.frame(x = 1:500, y = 1:500), mapping = aes(x, y)) +
  #scale_y_continuous(limits = c(0, 400)) +
  xlab("True value") +
  ylab("Parameter estimate") + 
  facet_wrap(.~parameter, ncol=3, scales = "free", labeller = labeller(parameter=label_parsed)) 
#stat_smooth(method="lm", se = FALSE)
p_all

ggsave(p_all, file = "figs/cv_plots_all.jpg", width = 9, height = 9)


# Ne bot cross-validation plots (just to have a better look, not included in final paper)

p_nbot <- cv_params %>%
  filter(parameter=="paste(italic(N)[e], \"bot\")") %>%
  ggplot(aes(true, estimate)) +
  geom_point(size = 2.3, alpha = 0.3) +
  geom_abline(intercept = 0, slope=1, color='blue') +
  xlab("True value") +
  ylab("Parameter estimate")  

p_nbot

ggsave(p_nbot, file = "figs/cv_plot_nbot.jpg", width = 4, height = 4)

p_nhist <- cv_params %>%
  filter(parameter=="paste(italic(N)[e], \"hist\")") %>%
  ggplot(aes(true, estimate)) +
  geom_point(size = 2.3, alpha = 0.3) +
  geom_abline(intercept = 0, slope=1, color='blue') +
  scale_x_continuous(limits = c(0, 200000)) + # limited to have a better look at firt part of plot, excludes 6 data points
  scale_y_continuous(limits = c(0, 200000)) +
  xlab("True value") +
  ylab("Parameter estimate")  

p_nhist
ggsave(p_nhist, file = "figs/cv_plot_nhist_200k.jpg", width = 4, height = 4) #excludes 6 data points
