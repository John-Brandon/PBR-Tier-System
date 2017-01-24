#
# Plot N_min percentile that meets MMPA management objective as a function
#  of the CV of the abundance estimates
# John Brandon
rm(list = ls())    # Clear workspace
library(magrittr)  # For Pipes %>% in R
library(dplyr)     # Data wrangling
library(ggplot2)   # Plotting
library(devtools)  # To source code snippet (aka 'gist') with custom plot theme
#
# Read in data.frame with solutions to Nmin percentile as function of CV_N -----
#
snmin = read.csv(file = "solve_Nmin.csv") %>% tbl_df() %>% 
  mutate(Tier = as.factor(Tier))

#
# Plot N_min precentile solution as function of CV_N ---------------------------
#

# Source R code from: https://gist.github.com/John-Brandon/484d152675507dd145fe
# Use package devtools to load plot theme code with `source_gist`
source_gist(id = "484d152675507dd145fe", filename = "mytheme_bw.R")

# Create plot 
ggplot(data = snmin, 
        aes(x = CV, 
        y = Nmin, 
        group = Tier, 
        fill = Tier, 
        color = Tier)) + 
  geom_line(lwd = 1.5) + 
  geom_point(size = 5, aes(shape = Tier, color = Tier)) +
  mytheme_bw + 
  scale_color_manual(values = c("green", "blue", "black")) +
  labs(y = expression(italic(N[MIN])*"  Percentile"), x = expression(CV[N])) +
  scale_x_continuous(breaks = seq(0.2, 1.0, by = 0.2)) +
  geom_vline(xintercept = 0.41, lty = 2) + 
  theme(legend.position = c(0.19, 0.85)) +
  theme(legend.title=element_blank()) 

# Save a high-res copy of the plot for publication
fig_path = ~/YourDirectory/  # where you want ggsave to write output plot 
fig_width = 9.58
fig_height = 8.39
ggsave(filename = "Figure_2_color.pdf", path = fig_path, width = fig_width, height = fig_height, dpi = 500)
