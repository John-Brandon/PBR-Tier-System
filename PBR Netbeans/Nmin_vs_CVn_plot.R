#
# Plot N_min percentile that meets MMPA management objective as a function
#  of the CV of the abundance estimates
#

# Read in data.frame with solutions to Nmin percentile as function of CV_N -----
library(magrittr)
library(dplyr)
library(ggplot2)

main_dir = getwd(); main_dir
tmp_wdir = "~/Documents/2015 work/PBR Tier System/Code/pbrr/R"

setwd(tmp_wdir)

snmin = read.csv(file = "solve_Nmin.csv")
snmin %<>% tbl_df() %>% mutate(Tier = as.factor(Tier))

setwd(main_dir)

# Plot N_min precentile solution as function of CV_N ---------------------------
plt = ggplot(data = snmin, aes(x = CV, y = Nmin, group = Tier, fill = Tier, color = Tier)) # , color = Tier
plt = plt + geom_line(lwd = 1.5) + mytheme_bw 
plt = plt + geom_point(size = 5, aes(shape = Tier, color = Tier))
# plt = plt + scale_color_manual(values = c("gray40", "gray60", "gray80"))  # For print version (B & W)
# plt = plt + scale_color_manual(values = c("green", "blue", "orange"))
# plt = plt + scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))
# plt = plt + scale_color_manual(values = c("#b2df8a", "#33a02c", "#1f78b4"))
# plt = plt + scale_color_manual(values = c("#a6cee3", "#33a02c", "#1f78b4"))  # For web version (color)
plt = plt + scale_color_manual(values = c("green", "blue", "black"))
plt = plt + labs(y = expression(italic(N[MIN])*"  Percentile"), x = expression(CV[N]))
plt = plt + scale_x_continuous(breaks = seq(0.2, 1.0, by = 0.2))
plt = plt + geom_vline(xintercept = 0.41, lty = 2)
plt = plt + theme(legend.position = c(0.19, 0.85))
plt = plt + theme(legend.title=element_blank())
plt

# Save a high-res copy of the plot for publication
fig_path = "~/Documents/2015 Work/PBR Tier System/Manuscript Text/Draft_for_journal/Revised Manuscript After Review/Revised-High-Res-Figs"
fig_width = 9.58
fig_height = 8.39
ggsave(filename = "Figure_2_color.pdf", path = fig_path, width = fig_width, height = fig_height, dpi = 500)
