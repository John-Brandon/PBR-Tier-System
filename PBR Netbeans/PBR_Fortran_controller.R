#----------------------------------------------------------------------------- #
# Code to run batches of trials, compile depletion statistics, and plot ouput
# OS: Mac OS 10.10.4
# R version 3.2.0 (2015-04-16)
# Author: John R. Brandon
#
# Notes:
#  1. Sourcing PBR_fortran_output.R also sources:
#     - PBR_FileIO.R
#     - PBR_create_input_files.R, and; 
#     - PBR_batch.R
#
#  2. Tier 1 has been reserved for an index of abundance approach in the code.
#     Tier 2 is coded as the standard approach / single abundance estimate (Wade, 1998 Mar Mamm Sci)
#     Tier 3 is coded as an arithmetic weighted average over the previous 8 yrs (NMFS 2005 GAMMS)
#     Tier 4 is coded as the weighted by time and precision approach (Brandao and Butterworth, 2014)
#
#  3. n_sims is set to 100 as a place-holder for the code repository. 
#     Batching runs with larger numbers 
#     (e.g n_sims = 2000) of simulations can take up a moderate amount 
#     (e.g. a few Gigabytes) of disk space, and several minutes (or more) to run. 
#     So, proceed with that in mind before running batches with large "n_sims". 
#
#  4. cet_file_names is a list of input files by trial for cetaceans,
#      created in PBR_create_input_files.R (see also pin_ and all_file_names)
#
#  5. Assumes these packages are installed (might be some redundancies):
#       ggplot2
#       stringr
#       tidyr
#       magrittr
#       dplyr 
#
#----------------------------------------------------------------------------- #
rm(list = ls())  # clear workspace

# Check installed packages and install if not already --------------------------
list_of_packages <- c("ggplot2", "stringr", "tidyr", "magrittr", "dplyr") 
new_packages = list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) install.packages(new_packages)

# Copy existing input.par file to temp file ------------------------------------
file.copy(from = "input.par", to = "input_par_copy.txt", overwrite = TRUE) 

# Create input files for trials and load helper functions ----------------------
# May throw warnings: 'incomplete final lines'. Safe to ignore.
source(file = "PBR_fortran_output.R") 

# Batch run of all cetacean trials ---------------------------------------------
lapply(X = cet_file_names, FUN = run_batch, n_sims = 100) # run all trials
# lapply(X = pin_file_names, FUN = run_batch, n_sims = 100) # all pinniped trials

# Store results from Tier 2 (standard approach) --------------------------------
dat_cet_trial02_0A = read.table(file = "N_agg_Cet_0A.out", header = TRUE)
dat_cet_trial02_0B = read.table(file = "N_agg_Cet_0B.out", header = TRUE)

# Compile depletion statistics -------------------------------------------------
df_depl = compile_depl_stats(file_names = cet_out_names, 
                             stock_id = 1, 
                             lower_percentile = 0.05, 
                             upper_percentile = 0.95)  # compile depletion

# Create some additional columns which are trial factors -----------------------
df_depl = tbl_df(df_depl)
df_depl %<>% mutate(taxa = substr(trial_id, start = 1, stop = 3),
                    trial_n = substr(trial_id, start = 5, stop = 5),
                    trial_l = substr(trial_id, start = 6, stop = 6),
                    CV_N = ifelse(trial_l %in% c("A", "C"), "a_low", "b_high"))
df_depl_cet = df_depl %>% filter(taxa == "Cet")
df_depl_cet_ab = df_depl_cet %>% filter(trial_l == "A" | trial_l == "B")
df_depl_cet_cd = df_depl_cet %>% filter(trial_l == "C" | trial_l == "D")
df_depl_pin = df_depl %>% filter(taxa == "Pin")
df_depl_pin_ab = df_depl_pin %>% filter(trial_l == "A" | trial_l == "B")
df_depl_pin_cd = df_depl_pin %>% filter(trial_l == "C" | trial_l == "D")

# Expressions for facet labels in zeh plots
trial_names = list( 
  'a_low' = expression(CV[N]*" = 0.20"),
  'b_high' = expression(CV[N]*" = 0.80")
)

trial_labeller = function(variable, value){ 
  # wrapper func for labels in zeh plots  
  return(trial_names[value])
}

# -- Plot depletion across trials (Zeh plots), c.f. Wade (1998) Fig 7. ---------
tmp_cet_plt = ggplot(data = df_depl_cet_cd, aes(x = trial_n, y = median)) + geom_point() + coord_cartesian(ylim = c(0,1.1))
tmp_cet_plt +
  facet_grid(CV_N ~ ., labeller = trial_labeller) + geom_errorbar(aes(ymin = lower_p, ymax = upper_p), width = 0.25) +
  geom_point(data = df_depl_cet_ab, aes(x = trial_n, y = lower_p), shape = 6) +
  geom_segment(aes(x = 0, y = 0.5, xend = 7.5, yend = 0.5), size = 1, colour = "red", linetype = "dashed") +
  geom_segment(aes(x = 7.5, y = 0.45, xend = 8.5, yend = 0.45), size = 1, colour = "red", linetype = "dashed") +
  geom_segment(aes(x = 8.5, y = 0.70, xend = 9.5, yend = 0.70), size = 1, colour = "red", linetype = "dashed") +
  xlab("Trial Number") + ylab("Depletion") + ylim(0, 1.5) +
  ggtitle(expression("Cetacean Depletion After 100 Years: "*F[R]*" = 0.50")) + mytheme_bw

# -- Depletion vs. N_min percentile plots --------------------------------------
# Cet_0A.txt is the base case trial with CV_N = 0.2
# See the notes in PBR_create_input_files.R for file naming conventions for trials
depl_matrix_cet0a_tier2 = batch_nmin(base_case_file = "Cet_0A.txt")
plot_depl_nmin(depl_matrix_cet0a_tier2, expression("Cetacean: "*CV[N]*" = 0.20"))

# Check out some spaghetti trajectory plots ------------------------------------
ribbon_depletion_plot(n_traj = 30, default_data = FALSE, N_agg = "N_agg_Cet_0A.out", stock_id = 1)

# Change Tier (GAMMS approach) and run another batch of trials -----------------
write_inits(par_name = "tier", par_val = formatC(3, format = "d"), # GAMMS approach
            infile = "input.par", outfile = "input.par")
write_inits(par_name = "n_yrs_avg", par_val = formatC(8, format = "d"), 
            infile = "input.par", outfile = "input.par")
source("PBR_create_input_files.R") # Create input files for each trial with two stocks (n = 72)
lapply(X = cet_file_names, FUN = run_batch, n_sims = 100) # run cetacean trials

# Store some results from Tier 3 (GAMMS approach) ------------------------------
dat_cet_trial03_0A = read.table(file = "N_agg_Cet_0A.out", header = TRUE)
dat_cet_trial03_0B = read.table(file = "N_agg_Cet_0B.out", header = TRUE)

# Compare inter-survey variation in PBR between tiers --------------------------
#   to compare standard approach with another tier
plot_prb_variation(default_data = FALSE, N_agg = "N_agg_Cet_0A.out", stock_id = 1)

# Make some worm plots of PBR through time -------------------------------------
tier2 = dat_cet_trial02_0A %>% tbl_df(.) %>% mutate(tier = 2) %>% filter(., sim %in% sample(1:pars$n_sims, 25))
tier3 = dat_cet_trial03_0A %>% tbl_df(.) %>% mutate(tier = 3) %>% filter(., sim %in% sample(1:pars$n_sims, 25))
tier2and3 = rbind(tier2, tier3)
tier2and3 %<>% filter(., stock == 1) %>% mutate(stock = as.factor(stock))
tier2and3  %<>% filter(yr > 1) %>% mutate(pbr_yr_sim_rescaled = pbr_yr_sim / pars$KK[1])
tier2and3  %<>% mutate(tier = as.factor(tier), sim = as.factor(sim), tiersim = as.factor(paste(tier,sim,sep="")))
levels(tier2and3$tier)[levels(tier2and3$tier) == "2"] = "Tier 1"
levels(tier2and3$tier)[levels(tier2and3$tier) == "3"] = "Tier 2"
names(tier2and3)[names(tier2and3) == "tier"] = "Tier"

ggplot(data = tier2and3, aes(x = yr, y = pbr_yr_sim_rescaled, colour = Tier, group = tiersim, alpha = Tier)) +
  geom_line() + mytheme_bw +
  scale_alpha_manual(values = c(0.8, 0.8)) +
  scale_color_manual(values = c("green", "blue")) +
  # scale_color_manual(values = c("darkgray", "black")) +
  labs(x = "Year", y = "PBR / K") + coord_cartesian(ylim = c(0, max(tier2and3$pbr_yr_sim_rescaled)*1.05)) +
  theme(legend.position = c(0,1), legend.justification = c(0,1), legend.title=element_blank())

# Example plots from two stock simulation --------------------------------------
write_inits(par_name = "n_stocks", par_val = formatC(2, format = "d"), # Set to two stocks
            infile = "input.par", outfile = "input.par")
source("PBR_create_input_files.R") # Create input files for each trial with two stocks (n = 72)
lapply(X = cet_file_names, FUN = run_batch, n_sims = 100) # run cetacean trials
plot_2stock_depl(read_data = FALSE, N_agg = "N_agg_Cet_0B.out") # Base case CV_N = 0.8; F_R = 1.0

# Revert input.par file (and trials) to Tier 1 and one stock -------------------
file.copy(from = "input_par_copy.txt", to = "input.par", overwrite = TRUE) 
source("PBR_create_input_files.R") # Re-write input files for each trial with one stock

