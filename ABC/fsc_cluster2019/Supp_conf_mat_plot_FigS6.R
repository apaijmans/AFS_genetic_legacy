# Make pretty Confusion Matrix plot
# Fig S6

library(forcats)
library(dplyr)
library(ggplot2)

# Read in data

conf_mat <- read.table(file = "model_evaluation/check2_models/sims10000k_cluster181_model_missclass_freq_3.txt", skip=1, header = F)

names(conf_mat) <- c("Simulated under", "Classified as", "Frequency")

levels(conf_mat$`Classified as`)
conf_mat$`Classified as` <- factor(conf_mat$`Classified as`, levels = rev(levels(conf_mat$`Classified as`)))
levels(conf_mat$`Classified as`)

# conf_mat <- data.frame("Bottleneck" = c(15, 85), "Non-bottleneck" = c(89, 11)) %>% 
#   gather(key = `Classified into`, value = Frequency) %>% 
#   mutate(`Simulated under` = c("Non-bottleneck", "Bottleneck", "Non-bottleneck", "Bottleneck")) %>% 
#   mutate(`Simulated under` = fct_relevel(`Simulated under`, c("Non-bottleneck", "Bottleneck")))

source("martin.R")

p <- ggplot(conf_mat, aes(x=`Simulated under`, y=Frequency, fill = `Classified as`)) +
  geom_bar(stat = "identity") +
  xlab("Simulated under") +
  scale_x_discrete(labels = c("Bottleneck", "Non-bottleneck")) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, by = 20)) +
  scale_fill_manual(name="Classified as", values = c("grey", "grey45"), breaks=c("bot","neut"), labels=c("Bottleneck", "Non-bottleneck")) +
  theme_martin(base_family = "Arial", highlight_family = "Arial") 

p

ggsave(filename = "model_evaluation/check2_models/sims10000k_cluster181_confusion_mat3.jpg",
       height = 3, width = 4.5)
