#----------------------------------------------------------------------------- #
# Code to read Fortran output file(s) for PBR Tier System
# Author   : John R. Brandon
# eMail    : jbrandon@gmail.com
# Date     : Fall 2015
# OS       : Mac OS 10.9.5
# R version: 3.2.0 (2015-04-16)
#
# Notes: 
# - Could compile and run the main program from R
# > compile_cmd = "gfortran A_Random_module.f90 BRENT.f90 Declare_variables_module.f90 Eigen_module.f90 Generate_random_numbers_module.f90 Initialize_pop_module.f90 PBR_Errorcheck_module.f90 PBR_FileIO_Module.f90 PBR_calcs_module.f90 main.f90 -o main -framework accelerate"
# > system(compile_cmd) # Link and compile main program
# > system("./main")    # Run main program
#----------------------------------------------------------------------------- #
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(magrittr)
library(dplyr)

# Fetch parameters from input.par file
source(file = "PBR_Batch.R") # Reads the same input file being used by Fortran code

# Theme for plotting
mytheme_bw = theme_bw() + theme(axis.title.x = element_text(size = rel(1.5), vjust = 0.0), 
                                axis.title.y = element_text(size = rel(1.5), vjust = 1.0), 
                                axis.text = element_text(size = rel(1.25), colour = "black"),
                                plot.title = element_text(size = rel(1.75)),
                                legend.text = element_text(size = 20), 
                                legend.title = element_text(size = 20),
                                panel.grid.major = element_line(colour = "darkgray"),
                                panel.grid.minor.y = element_blank(),
                                strip.text.x = element_text(size = 15), 
                                strip.text.y = element_text(size = 15))

# ---------------------------------------------------------------------------- #
ribbon_depletion_plot = function(n_traj = 25, default_data = TRUE, N_agg = NULL, stock_id = 1, title = FALSE){
###
# Plot depletion levels from a random sample of size = n_traj from the set of simulations 
# The simulation with index zero is the reference simulation  
# Stock "zero" = abundance of stock 1 summed across areas
#  TODO: add a "mutate" command to sum stock 1 & 2 seperately across areas for multi-stock trials  
###  
  pars = NULL; sims_sample = NULL; ref_sim = NULL; trial_id = NULL; mnpl = NULL

  if (default_data){
    pars = read_inits() # read input parameters for this trial  
    N_agg = read.table(file = "N_aggregated.out", header = TRUE)
    N_agg = tbl_df(N_agg) 
    # trial_id = "Default N_aggregated.out"
  } else if(is.null(N_agg)){
    stop("Error: Need to supply 'N_agg' if argument 'read_data' is FALSE")
  } else {
    trial_id = substr(N_agg, start = 7, stop = 12)
    in_file = paste(trial_id, ".txt", sep = "")
    pars = read_inits(input_file = in_file) # read input parameters for this trial  
    N_agg = read.table(file = N_agg, header = TRUE)
    N_agg = tbl_df(N_agg)  
  }
  
  if(n_traj > pars$n_sims) stop(paste("n_traj greater than number of available simulations:", pars$n_sims))
  ref_sim = 0 # reference simulation with zero human caused mortality
  sims_sample = sample(1:pars$n_sims, n_traj) # random sample of n_traj from set
  sims_sample = c(ref_sim, sims_sample)

  N_agg$stock = as.factor(N_agg$stock) # Convert stock ID number into factor

  N_agg_sims = N_agg %>% filter(., stock == stock_id)  # TODO: Extend 'stock' -- Developing, starting with single stock scenario
  N_agg_sims = N_agg_sims %>% mutate(., ref_case = ifelse(sim == 0, TRUE, FALSE))
  N_agg_sims_sample = N_agg_sims %>% filter(., sim %in% sims_sample) # just a sample to plot
  
# Calculate annual point-wise depletion percentiles (5th and 95th) across simulations
  percentiles = filter(N_agg_sims, sim > 0) %>% group_by(., yr) %>%
    summarise(., median = median(depl_yr_stock),
                lower_5th = quantile(depl_yr_stock, 0.05),
                upper_95th = quantile(depl_yr_stock, 0.95))

  plt = ggplot() + 
    geom_line(data = N_agg_sims_sample, aes(x = yr, y = depl_yr_stock, group = sim, 
      linetype = ref_case, size = ref_case, col = ref_case, alpha = ref_case)) + 
    scale_alpha_manual(values = c(0.50, 1.0)) +
    scale_size_manual(values = c(0.75, 1.0)) +     
    scale_color_manual(values = c("black", "red")) + # Can set different colors here
    scale_linetype_manual(values = c(1, 1)) +          # Can set different line-types for ref_case (TRUE or FALSE)
    labs(x = "Year", y = "Depletion", title = paste(trial_id, "[", pars$ref, "]")) + 
    theme_bw() + theme(legend.position = "none") + 
    geom_ribbon(data = percentiles, aes(x = yr, ymin = lower_5th, ymax = upper_95th), alpha = 0.30)  +
      coord_cartesian(ylim = c(0, 1.1)) + mytheme_bw + theme(legend.position="none")
  
  if(pars$theta == 1.0){
    plt + geom_hline(aes(yintercept = 0.50), colour = "red", linetype = "dashed", lwd = 1.2)
  }else if(pars$theta == 0.53){
    plt + geom_hline(aes(yintercept = 0.45), colour = "red", linetype = "dashed", lwd = 1.2)
  }else if(pars$theta == 5.04){
    plt + geom_hline(aes(yintercept = 0.70), colour = "red", linetype = "dashed", lwd = 1.2)
  }else{
    plt + geom_hline(aes(yintercept = 0.50), colour = "blue", linetype = "dashed", lwd = 1.2)
  }
}

# ---------------------------------------------------------------------------- #
survey_plot_agg = function(sim_n, read_data = TRUE, N_agg = NULL, stock_id = 1){
#
  pars = read_inits()
#  
  if (read_data){
    trial_id = pars$ref
    N_agg = read.table(file = "N_aggregated.out", header = TRUE)
    N_agg = tbl_df(N_agg)  
  } else if(is.null(N_agg)){
    stop("Error: Need to supply 'N_agg' if argument 'read_data' is FALSE")
  } else{
    # trial_id = substr(N_agg, start = 7, stop = 12)
    trial_id = pars$ref
    N_agg = read.table(file = N_agg, header = TRUE)
    N_agg = tbl_df(N_agg)  
  }
  
  N_agg$stock = as.factor(N_agg$stock) # Convert stock ID number into factor
  
  N_agg_surveys = N_agg %>% filter(., stock == 0, sim == sim_n, n_hat_yr > 0) %>% 
    mutate(., lower2.5 = qlnorm(0.025, meanlog = log(n_hat_yr), sdlog = sqrt(log(1 + pars$cv_n[1]^2))),
              upper2.5  = qlnorm(0.975, meanlog = log(n_hat_yr), sdlog = sqrt(log(1 + pars$cv_n[1]^2))) )  
  
  # N_agg_trajectories
  N_agg_first_sim = N_agg %>% filter(., sim == sim_n)
  
  ggplot(data = N_agg_first_sim, aes(x = yr, y = N_tot_area123, group = stock, col = stock)) + 
    geom_line() +
    geom_point(data = N_agg_surveys, aes(x = yr, y = n_hat_yr)) +
    geom_errorbar(data = N_agg_surveys, aes(x = yr, ymin = lower2.5, ymax = upper2.5, width = 1.2)) +
    expand_limits(y = 0) +  
    coord_cartesian(xlim = range(N_agg_first_sim$yr), ylim = c(0, max(N_agg_surveys$upper2.5) * 1.1)) +
    labs(x = "Year", y = "Abundance", title = trial_id ) + # "Abundance in Areas 1, 2, and 3"
    theme_bw() +  theme(legend.position="none") #theme(legend.position = c(0.1, 0.8))
}

# ---------------------------------------------------------------------------- #
plot_depl_nmin = function(depl, title){
#  
# `depl` is a data.frame created by the `batch_nmin()` function in the PBR_batch.R script
#
  lower_tail = as.numeric(names(depl))
  
  plot(0, 0, xlab=expression(N[MIN]*" Percentile"), ylab="Final Depletion", type="n", xlim=c(0.0,0.55), ylim=c(0,1.1), yaxs="i", xaxs="i")
  
  for (ii in 1:length(lower_tail)){
    Box(depl[,ii], ii = ii, xx = lower_tail)
  } 
  for (ii in 2:length(lower_tail)){
    lines(x = c(lower_tail[ii-1], lower_tail[ii]), y = c(median(depl[,ii-1]), median(depl[,ii])))
  }
  title(title)
}

# ---------------------------------------------------------------------------- #
plot_depl_nmin2 = function(depl1, depl2, title){
  #  Modified version of function to plot two lines (i.e. compare 2 tiers) on same graph
  # `depl` is a data.frame created by the `batch_nmin()` function in the PBR_batch.R script
  #
  lower_tail1 = as.numeric(names(depl1))
  lower_tail2 = as.numeric(names(depl2))
  
  plot(0, 0, xlab=expression(N[MIN]*" Percentile"), ylab="Depletion", type="n", xlim=c(0.0,0.55), ylim=c(0,1.1), yaxs="i", xaxs="i")
  
  for (ii in 1:length(lower_tail1)){
    Box(depl1[,ii], ii = ii, xx = lower_tail1)
  } 
  for (ii in 2:length(lower_tail1)){
    lines(x = c(lower_tail1[ii-1], lower_tail1[ii]), y = c(median(depl1[,ii-1]), median(depl1[,ii])), col = "green")
  }
  for (ii in 1:length(lower_tail2)){
    Box(depl2[,ii], ii = ii, xx = lower_tail2)
  } 
  for (ii in 2:length(lower_tail2)){
    lines(x = c(lower_tail2[ii-1], lower_tail2[ii]), y = c(median(depl2[,ii-1]), median(depl2[,ii])), col = "blue")
  }
  title(title)
}

# ---------------------------------------------------------------------------- #
# "Zeh" style plots for depletion as function of N_min(lower_tail), cf Wade (1998) Fig 4.
Box <- function(Dets, ii, xx, pch=16){
  points(x = xx[ii], y = median(Dets), pch = 16)
  lines(c(xx[ii], xx[ii]), c(quantile(Dets, 0.05), quantile(Dets, 0.95)), lty = 1)
  abline(h = 0.5, lty = 2, col = "black")
}

# ---------------------------------------------------------------------------- #
plot_prb_variation = function(default_data = TRUE, N_agg = NULL, stock_id = 1){
  # Calculate Average inter-survey variation and compare between results from naive tier (1) and another tier
  # Assumes lower 20th percentile. PBR_naive is the standard approach using only last abundance estimate.
  
  if (default_data){
    pars = read_inits() # read input parameters for this trial  
    N_agg = read.table(file = "N_aggregated.out", header = TRUE)
    N_agg = tbl_df(N_agg)  
  } else if(is.null(N_agg)){
    stop("Error: Need to supply 'N_agg' if argument 'read_data' is FALSE")
  } else {
    trial_id = substr(N_agg, start = 7, stop = 12)
    in_file = paste(trial_id, ".txt", sep = "")
    pars = read_inits(input_file = in_file) # read input parameters for this trial  
    N_agg = read.table(file = N_agg, header = TRUE)
    N_agg = tbl_df(N_agg)  
  }
  
  tier_n = pars$tier - 1
  
  N_agg %<>% filter(., stock == 1, n_hat_yr > 0, sim > 0)  
  
  # Calculate standard PBR given single N_hat in a given survey year    
  pbr_df = N_agg %>% mutate(
    pbr_naive = 0.5 * pars$R_max * pars$F_r[1] * calc_n_min(n_best = n_hat_yr, cv_n = pars$cv_n[1], z_score = abs(qnorm(pars$lower_tail)))
  )
  
  pbr_df %<>% group_by(., sim) %>% mutate(pbr_naive_dif = pbr_naive - lag(pbr_naive), 
                                          pbr_yr_dif = pbr_yr_sim - lag(pbr_yr_sim))
  
  pbr_df %<>% select(sim, pbr_naive_dif, pbr_yr_dif) 
  
  #   aav_naive = sum(abs(tmp_df$pbr_naive_dif), na.rm = TRUE) / sum(tmp_df$pbr_naive, na.rm = TRUE); aav_naive
  #   aav_tier = sum(abs(tmp_df$pbr_yr_dif), na.rm = TRUE) / sum(tmp_df$pbr_yr_sim, na.rm = TRUE); aav_tier
  
  pbr_df %<>% gather(tier, dif, -sim) # gather reshapes data.frame in preperation for ggplot
  
  pbr_df %<>% filter(!is.na(dif)) %>% 
    mutate(tier = as.factor(tier)) %>% 
    mutate(pbr_over_k = dif / pars$KK[1],  #, dif = abs(dif)
           aav = abs(pbr_over_k))  
  
  pbr_df$tier = ifelse(pbr_df$tier == "pbr_naive_dif", "Tier 1", paste("Tier", tier_n))
  
  levels(pbr_df$tier)[levels(pbr_df$tier) == "pbr_naive_dif"] = "Tier 1"
  levels(pbr_df$tier)[levels(pbr_df$tier) == "pbr_yr_dif"] = paste("Tier ", tier_n)
  
  # Plot the inter-survey variation in PBR
  plt = ggplot(data = pbr_df, aes(x = pbr_over_k, fill = tier)) # 
  plt + geom_histogram(position = "identity", alpha = 0.7, binwidth = .0002) + # , aes(y=..count../sum(..count..))
    scale_alpha_manual(values = c(0.5, 0.5)) +  
    scale_fill_manual(values = c("green", "blue")) +
    labs(x = "Inter-survey difference in PBR / K", y = "Frequency", title = "") +
    mytheme_bw + 
    theme(legend.position = c(0,1), legend.justification = c(0,1), legend.title=element_blank()) 
  
}

# ---------------------------------------------------------------------------- #
compile_depl_stats = function(file_names, stock_id = 1, lower_percentile = 0.05, upper_percentile = 0.95){
  #
  # Read each output file in file_names[], and compile depletion statistics (5th, median and 95th percentiles)  
  #  
  ii = NULL; out_file_names = NULL; max_sim = NULL; n_trials = NULL
  
  n_trials = length(file_names)
  depl_df = data.frame(
    trial_id = character(n_trials),
    lower_p = double(n_trials),
    median = double(n_trials),
    upper_p = double(n_trials), 
    stringsAsFactors = FALSE) 
  
  file_check(file_names = file_names) # check all files in list exist
  trial_id = substr(file_names, start = 7, stop = 12) # extract trial_id
  
  for(ii in 1:length(file_names)){
    print(paste("Starting trial", ii, "of", length(file_names),":", trial_id[ii]))
    N_agg = read.table(file = file_names[ii], header = TRUE)
    N_agg = tbl_df(N_agg)
    N_agg$stock = as.factor(N_agg$stock) # Convert stock ID number into factor
    # Start extracting relevant data    
    N_agg_sims = N_agg %>% filter(., stock == stock_id)  # TODO: Extend 'stock' -- Developing, starting with single stock scenario
    N_agg_sims = N_agg_sims %>% filter(., sim > 0) # exclude the reference set 
    final_depl = N_agg_sims %>% filter(., yr == 100) %>% select(., depl_yr_stock) %>% arrange(., depl_yr_stock)
    max_sim = max(N_agg_sims$sim)
    lower_ii = max_sim * lower_percentile
    upper_ii = max_sim * upper_percentile
    upper_ii = upper_ii + 1
    # Assign percentiles to data.frame
    depl_df$trial_id[ii] = trial_id[ii]
    depl_df$lower_p[ii] = final_depl$depl_yr_stock[lower_ii] # Extract lower_percentile 
    depl_df$median[ii] = median(final_depl$depl_yr_stock)    # Extract median 
    depl_df$upper_p[ii] = final_depl$depl_yr_stock[upper_ii] # Extract upper_percentile 
  }
  return(depl_df)
}

# ---------------------------------------------------------------------------- #
worm_plot = function(n_traj = 30, default_data = TRUE, N_agg = "N_aggregated.out", stock_id = 1){

  if (default_data){
    pars = read_inits() # read input parameters for this trial  
    N_agg = read.table(file = "N_aggregated.out", header = TRUE)
    N_agg = tbl_df(N_agg)  
  } else if(is.null(N_agg)){
    stop("Error: Need to supply 'N_agg' if argument 'read_data' is FALSE")
  } else {
    trial_id = substr(N_agg, start = 7, stop = 12)
    in_file = paste(trial_id, ".txt", sep = "")
    pars = read_inits(input_file = in_file) # read input parameters for this trial  
    N_agg = read.table(file = N_agg, header = TRUE)
    N_agg = tbl_df(N_agg)  
  }

  if(n_traj > pars$n_sims) stop(paste("n_traj greater than number of available simulations:", pars$n_sims))  
  N_agg$stock = as.factor(N_agg$stock) # Convert stock ID number into factor
  
  N_agg_sims = N_agg %>% filter(., stock == stock_id)  # TODO: Extend 'stock' -- Developing, starting with single stock scenario
  sims_sample = sample(1:pars$n_sims, n_traj) # random sample of n_traj from set
  N_agg_sims_sample = N_agg_sims %>% filter(., sim %in% sims_sample) # just a sample to plot
  N_agg_sims_sample  %<>% filter(yr > 1) %>% mutate(pbr_yr_sim_rescaled = pbr_yr_sim / pars$KK[1])
  ggplot(data = N_agg_sims_sample, aes(x = yr, y = pbr_yr_sim_rescaled, group = sim)) + geom_line() + mytheme_bw +
    labs(x = "Year", y = "PBR / K") + coord_cartesian(ylim = c(0, 0.050))
  
}

# ---------------------------------------------------------------------------- #
calc_n_min = function(n_best, cv_n, z_score){
  ####### +++ ### ### +++ ### ### +++ ### 
  #====== +++ === === +++ === === +++ === 
  # Function to calculate N_min, given (i) CV of abundance estimate (ii) N_best and (iii) standard normal variate, z
  # Uses equation (4) of Wade(1998 Mar Mamm Sci)
  # Note that N_best is the expectation (mean) of a log-normal distribution, not the median.        
  # Positive values for z-score return lower tail the way the Eqn 4 of Wade (1998) is derived. 
  # For example, for the lower 2.5th percentile, z_score = 1.96 (not -1.96) 
  # Of course, you could just write this function in one line of R-code
  #  <- qlnorm(p = z_score, meanlog = log(n_best), sdlog = sqrt(log(1 + cv_n * cv_n)))
  #====== +++ === === +++ === === +++ === 
  if (z_score < 0) {
    print("Error: Z-score assumed to be absolute value")
    print("Changing to absolute value")
    print("You should check that resulting number is correct")
    z_score = abs(z_score)
  }
  calc_n_min = log(1 + cv_n * cv_n)   # Start by calculating denominator
  calc_n_min = sqrt(calc_n_min)       # Standard deviation in log-space
  calc_n_min = z_score * calc_n_min
  calc_n_min = exp(calc_n_min)
  calc_n_min = n_best / calc_n_min # divide N_best by denominator
  
  return(calc_n_min)
}

# ---------------------------------------------------------------------------- #
plot_2stock_depl = function(n_traj = 30, read_data = TRUE, N_agg = NULL){
  #
  pars = read_inits()
  #  
  if (read_data){
    trial_id = pars$ref
    N_agg = read.table(file = "N_aggregated.out", header = TRUE)
    N_agg = tbl_df(N_agg)  
  } else if(is.null(N_agg)){
    stop("Error: Need to supply 'N_agg' if argument 'read_data' is FALSE")
  } else{
    # trial_id = substr(N_agg, start = 7, stop = 12)
    trial_id = pars$ref
    N_agg = read.table(file = N_agg, header = TRUE)
    N_agg = tbl_df(N_agg)  
  }
  
  if(n_traj > pars$n_sims) stop(paste("n_traj greater than number of available simulations:", pars$n_sims))
  # ref_sim = 0 # reference simulation with zero human caused mortality
  sims_sample = sample(1:pars$n_sims, n_traj) # random sample of n_traj from set
  
  # Remove reference simulation with zero human caused mortality 
  N_agg_sims = N_agg %>% filter(., sim > 0 & stock > 0) %>% mutate(stock = as.factor(stock))
#   N_agg_sims = N_agg %>% filter(., stock == stock_id)  # TODO: Extend 'stock' -- Developing, starting with single stock scenario
#   N_agg_sims = N_agg_sims %>% mutate(., ref_case = ifelse(sim == 0, TRUE, FALSE))
  N_agg_sims_sample = N_agg_sims %>% filter(., sim %in% sims_sample) # just a sample to plot
  # Calculate annual point-wise depletion percentiles (5th and 95th) across simulations
  percentiles = filter(N_agg_sims, sim > 0) %>% group_by(., yr, stock) %>% 
    summarise(., median = median(depl_yr_stock),
              lower_5th = quantile(depl_yr_stock, 0.05),
              upper_95th = quantile(depl_yr_stock, 0.95))
  names(percentiles)[names(percentiles) == "stock"] = "Stock" # Prettier for plotting
  names(N_agg_sims_sample)[names(N_agg_sims_sample) == "stock"] = "Stock" # Prettier for plotting
  
  N_agg_sims_sample$depl_yr_stock
  
  ggplot() +
  geom_line(data = N_agg_sims_sample, aes(x = yr, y = depl_yr_stock, group = interaction(sim, Stock), col = Stock)) +        
  expand_limits(y = 0) + 
  scale_alpha_manual(values = c(0.50, 0.5)) +
  coord_cartesian(xlim = range(N_agg_sims$yr), ylim = c(0, 1.05)) +
  labs(x = "Year", y = "Depletion", title = "" ) + # "Abundance in Areas 1, 2, and 3"
  theme_bw() + mytheme_bw +
  geom_hline(aes(yintercept = 0.50), colour = "black", linetype = "dashed", lwd = 1.2) +
  theme(legend.position = c(0,1), legend.justification = c(0,1)) + # , legend.title=element_blank()    
  scale_color_manual(values = c("blue", "red")) +
  geom_ribbon(data = percentiles, aes(x = yr, ymin = lower_5th, ymax = upper_95th, fill = Stock), alpha = 0.40) +
  scale_fill_manual(values = c("blue", "red")) #+ # Can set different colors here

}

