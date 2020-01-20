# Figure 4d

library(data.table)
library(dplyr)
library(ggplot2)


# path to simulation created in script 10_simulations_heatmap
sim_name <- "sims_heatmap"

# load simulations -------------------------------------------------------------

path_to_sims <- paste0(sim_name,".txt")
sims <- fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)


# calculate mean num of alleles-------------------------------------------------

sim_divs <- sims %>%
  mutate(bot_dur = tbotstart-tbotend) %>%
  mutate(nbot = as.numeric(gsub(10000, 0, nbot))) %>% # for plotting we 'correct' the nbot for the neutral model to 0
  group_by(nbot, model, bot_dur) %>%
  summarise(mean_alleles = mean(num_alleles_mean)) 


# calculate diversity loss -----------------------------------------------------

sim_divs <- sim_divs %>%
  group_by(bot_dur) %>% 
  mutate(alleles_lost = mean_alleles[nbot == 0] - mean_alleles) %>% 
  mutate(div_lost_perc = 100 * (1- (mean_alleles/mean_alleles[nbot == 0])))



# plot diversity loss ----------------------------------------------------------

# Fig 4d
source("martin.R")

p1 <- sim_divs %>%
  filter(nbot > 0 & nbot<525) %>% # use only nbot values between 25-500
  ggplot(aes(bot_dur, nbot)) + 
  geom_tile(aes(fill = div_lost_perc), colour = "white") + 
  scale_fill_gradientn(colors = c("#4575b4", "#a6d96a", "#ffffbf", "#fdae61", "#d73027")) +
  #scale_fill_distiller(name = "Diversity\nlost (%)", palette = "RdYlBu")+
  labs(x = "Bottleneck duration (generations)", y = expression(paste(italic(N)[e], "bot")), fill = "Diversity\nloss (%)") + 
  theme_martin(base_family = "Arial", highlight_family = "Arial")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.ticks=element_line(size=0.4))

p1

ggsave(p1, file = "figs/heatmap_perc_div_loss.jpg", width = 7, height = 6)

# min(sim_divs$div_lost_perc) #0
# max(sim_divs$div_lost_perc) #56
# 
# scale_fill_distiller(name = "Diversity\nlost (%)", palette = "RdYlBu")
# #4575b4 blue
# #66bd63 green
# #ffffbf yellow
# #d73027 red

# # plot mean num of alleles------------------------------------------------------
# # This plot was not used in the paper
# 
# source("martin.R")
# p2 <- ggplot(sim_divs, aes(bot_dur, nbot)) + 
#   geom_tile(aes(fill = mean_alleles), colour = "white") + 
#   scale_fill_gradient2(low = "#E7298A", mid = "white", high = "grey50", midpoint = 5.5) +
#   labs(x = "Bottleneck duration (generation)", y = expression(paste(italic(N)[e], "bot")), fill = "Mean number\nof alleles") + 
#   theme_martin(base_family = "Arial", highlight_family = "Arial")+
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(axis.ticks=element_line(size=0.4))
# 
# p2
# 
# ggsave(p2, file = "figs/heatmap_num_alleles.jpg", width = 9, height = 9)
# 
# min(sim_divs$mean_alleles) #3.28
# max(sim_divs$mean_alleles) #7.52
