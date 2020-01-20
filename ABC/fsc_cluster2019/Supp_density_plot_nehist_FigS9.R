
# Make density plots for Nehist (Fig S9) per cluster ("pop")

library(tidyr)
library(dplyr)


# Nhist to long data for plotting all together in one plot
load("abc_estimates/abc_sims10000k_cluster181_bot_ALL_3.RData")

nhist_long <- abc_bot[c("nhist", "pop")] %>%
  gather("nhist", key = "stat", value = "values")


# Order data for plot
nhist_long$pop <- ordered(nhist_long$pop, levels=c("SSI", "SG", "BI", "MAR_CI", "HI_KI_MAC"))

col <- c( "#E7298A","#6A3D9A","#B2DF8A", "#CFC855", "#FFA500")

full_pop <- c(SSI = "South\nShetlands", 
              SG = "South\nGeorgia", 
              BI = "Bouvetøya", 
              MAR_CI = "Marion Island &\nCrozet Islands", 
              HI_KI_MAC = "Kerguelen Islands,\nHeard Island &\nMacquarie Island")


source("martin.R")

# Function to calculate the mode
#library(ggbeeswarm)
library(ggplot2)
library(grid)
library(ggridges)

estimate_mode <- function(s) {
  d <- density(s, adjust = 2.5)
  d$x[which.max(d$y)]
}

# Function to calculate CI
calc_CI <- function(x, CI) {
  out <- stats::quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
  names(out) <- c("lower_ci", "upper_ci")
  out
}


# Add priors for Nehist

priors <- round(round(rlnorm(25000, 10.5, 1)))

nhist_all <- cbind(nhist_long, priors)


# Squareroot transformed plot for Nehist

# Function to use boxplot.stats to set the box-and-whisker locations
# This function will calculate the values for the whiskers such that they display the 'original' values 
# (ie not calculated over the transformed data, but over the 'normal' data) in the log transformed plot
mybxp = function(x) {
  bxp = sqrt(boxplot.stats(x^2)[["stats"]])
  names(bxp) = c("ymin","lower", "middle","upper","ymax")
  return(bxp)
}

p <- nhist_all %>%
  ggplot() +
  
  stat_density(geom="line", data= nhist_all, aes(x=priors), linetype="dashed")+
  stat_density(geom="line", data= nhist_all, aes(x=values, colour= pop, fill=pop))+
  
  scale_fill_manual(values=col)+ 
  scale_colour_manual(values=col)+
  
  facet_grid(pop ~ ., labeller = labeller(pop=full_pop)) +
  
  xlab("Nehist") +
  ylab("Density") +
  scale_x_continuous(limits = c(0, 100000)) +
  #scale_y_continuous(labels = scales::comma) +
  
  guides(fill=FALSE, colour=FALSE)+ #remove legend (here use color not fill, as we used color to play with colours)
  
  theme_martin(base_family = "Arial", highlight_family = "Arial")+
  theme(axis.text.y= element_text(size = 9.5),
        axis.text.x= element_text(size = 7.5), #changes size xaxis and yaxis so pops, and values. But not top labels or "Value"
        strip.text.y = element_text(angle = 0))
p

ggplot2::ggsave(filename = "figs/Nhist_posteriors_density.jpg", p, width = 7, height = 7)
