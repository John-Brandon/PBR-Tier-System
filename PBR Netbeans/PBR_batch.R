#----------------------------------------------------------------------------- #
# Functions for running batches: PBR Tier System
# Author   : John R. Brandon
# eMail    : jbrandon@gmail.com
# Date     : Fall 2015
# OS       : Mac OS 10.9.5
# R version: 3.2.0 (2015-04-16)
#
# Notes:
# To compile the Fortran code on Mac OS:
# gfortran A_Random_module.f90 BRENT.f90 Declare_variables_module.f90 Eigen_module.f90 Generate_random_numbers_module.f90 Initialize_pop_module.f90 PBR_Errorcheck_module.f90 PBR_FileIO_Module.f90 PBR_calcs_module.f90 main.f90 -o main -fbounds-check -framework accelerate
#----------------------------------------------------------------------------- #
library(dplyr)
library(ggplot2)
library(tidyr)   # For function `seperate`: Split trial_id column into two columns with taxa and trial factors
library(stringr)

# Fetch parameters from input.par file
# code_dir = "~/Documents/2015 Work/PBR Tier System/Code/R Code" # You'll need to change this (apologies for hard-coding)
# setwd(code_dir)
source(file = "PBR_FileIO.R") # Reads the same input file being used by Fortran code
source(file = "PBR_create_input_files.R")

# copy existing input.par file to temp file
file.copy(from = "input.par", to = "input_par_copy.txt", overwrite = TRUE) 

# Set-up list of files to run as a batch -- these vectors created in PBR_create_input_files.R
all_file_names = c(a_file_names, b_file_names, c_file_names, d_file_names)
cet_file_names = subset(all_file_names, substr(all_file_names, 1, 3) == "Cet")
pin_file_names = subset(all_file_names, substr(all_file_names, 1, 3) == "Pin")

# Create list with output file names
all_out_names = sub(x = all_file_names, pattern = "txt", replacement = "out")
all_out_names = paste("N_agg_", all_out_names, sep = ""); all_out_names
cet_out_names = sub(x = cet_file_names, pattern = "txt", replacement = "out")
cet_out_names = paste("N_agg_", cet_out_names, sep = ""); cet_out_names
pin_out_names = sub(x = pin_file_names, pattern = "txt", replacement = "out")
pin_out_names = paste("N_agg_", pin_out_names, sep = ""); pin_out_names

create_out_names = function(in_files_names){
  out_names = sub(x = in_files_names, pattern = "txt", replacement = "out")
  out_names = paste("N_agg_", out_names, sep = "")
  return(out_names)
}

### Run a batch of trials ------------------------------------------------------
run_batch = function(file_names, n_sims = 100, ...){
  #
  # Run a batch of simulation trials
  #

  n_trials = length(file_names)
  trial_depl = data.frame(trial_id = character(n_trials),
                          lower = numeric(n_trials),
                          median = numeric(n_trials),
                          upper = numeric(n_trials),
                          stringsAsFactors = FALSE)

  trial_id = substr(file_names, start = 1, stop = 6)

  file_check(file_names) # will stop and throw error if input file(s) don't exist

  file.copy(from = "input.par", to = "input_par_copy.txt", overwrite = TRUE) # copy existing input.par file to temp file

  for(ii in 1:length(file_names)){ # Loop over input files (trials)
    if(file_names != "input.par"){
      file.copy(from = file_names[ii], to = "input.par", overwrite = TRUE) # re-establish original input.par file
    }
    # file.copy(from = file_names[ii], to = "input.par", overwrite = TRUE)
    write_inits(par_name = "n_sims", par_val = formatC(n_sims, format = "d", digits = 3),
                infile = "input.par", outfile = "input.par")
    system2("./main", ...) # Mac OS system command to run main program
    trial_depl[ii,] = read.table(file = "depl_by_trial.out", stringsAsFactors = FALSE, header = TRUE)
    file.copy(from = "N_aggregated.out", to = paste("N_agg_", trial_id[ii], ".out", sep = ""), overwrite = TRUE)
  }

  file.copy(from = "input_par_copy.txt", to = "input.par", overwrite = TRUE) # re-establish original input.par file

  trial_depl %<>% mutate(taxa = substr(trial_id, start = 1, stop = 3),
                         trial_n = substr(trial_id, start = 5, stop = 5),
                         trial_l = substr(trial_id, start = 6, stop = 6)
  )

  return(trial_depl)
}

### Run base case trial over range of percentiles for N_min --------------------
batch_nmin = function(base_case_file = "Cet_0A.txt"){
#
# Compile depletion percentiles across range of lower_tail percentiles used to calculate N_min c.f. Wade Fig 4.
#
  file.copy(from = "input.par", to = "input_par_copy_nmin.txt", overwrite = TRUE) # copy old input.par
  if(base_case_file != "input.par"){
    file.copy(from = base_case_file, to = "input.par", overwrite = TRUE) # re-establish original input.par file
  }

  lower_tail = c(0.025, 0.05)
  lower_tail = c(lower_tail, seq(from = 0.10, to = 0.50, by = 0.05))

  pars = read_inits()
  depl = as.data.frame(matrix(0, nrow = pars$n_sims, ncol = length(lower_tail)))
  names(depl) = as.character(lower_tail)

  # Loop over percentiles for N_min
  for(ii in 1:length(lower_tail)){
    write_inits(par_name = "lower_tail", par_val = formatC(lower_tail[ii], format = "f", digits = 3),
                infile = "input.par", outfile = "input.par")
    system2("./main") # Mac OS system command to run main program
    N_agg = NULL
    N_agg = read.table(file = "N_aggregated.out", header = TRUE)
    N_agg = tbl_df(N_agg)
    depl[,ii] = N_agg %>% filter(., sim > 0, stock == 1, yr == 100) %>% dplyr::select(., depl_yr_stock)
    depl[,ii] = sort(depl[,ii])
  }
  file.copy(from = "input_par_copy_nmin.txt", to = "input.par", overwrite = TRUE) # re-establish original input.par file
  return(depl)
}

### Run base case trial over range of initial depletion levels -----------------
batch_init_depl = function(base_case_file = "Cet_0A.txt"){
  #
  # Compile final depletion as a function of initial depletion percentiles
  #
  file.copy(from = "input.par", to = "input_par_copy_nmin.txt", overwrite = TRUE) # copy old input.par
  if(base_case_file != "input.par"){
    file.copy(from = base_case_file, to = "input.par", overwrite = TRUE) # re-establish original input.par file
  }

  init_depl = c(0.025, 0.05)
  init_depl = c(init_depl, seq(from = 0.10, to = 1.0, by = 0.1))

  pars = read_inits()
  depl = as.data.frame(matrix(0, nrow = pars$n_sims, ncol = length(init_depl)))
  names(depl) = as.character(init_depl)

#  Loop over percentiles for N_min
  for(ii in 1:length(init_depl)){
    print(paste("Initial Depletion:", init_depl[ii]))
    write_inits(par_name = "init_depl", par_val = formatC(init_depl[ii], format = "f", digits = 3),
                infile = "input.par", outfile = "input.par")
    system2("./main") # Mac OS system command to run main program
    N_agg = NULL
    N_agg = read.table(file = "N_aggregated.out", header = TRUE)
    file.copy(from = "N_aggregated.out", to = paste("N_agg_depl_init_", init_depl[ii], ".out", sep = ""), overwrite = TRUE)
    N_agg = tbl_df(N_agg)
    depl[,ii] = N_agg %>% filter(., sim > 0, stock == 1, yr == 100) %>% dplyr::select(., depl_yr_stock)
    depl[,ii] = sort(depl[,ii])
  }
  file.copy(from = "input_par_copy_nmin.txt", to = "input.par", overwrite = TRUE) # re-establish original input.par file
  return(depl)
}


