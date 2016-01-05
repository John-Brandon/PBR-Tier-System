#===#===#===#===#===#===#===#===#===#===#
# Do some file IO
#===#===#===#===#===#===#===#===#===#===#
# rm(list = ls()) # Clear workspace
# dat_dir = "~/Documents/2015 Work/PBR Tier System/Code/PBR Netbeans/PBR Netbeans"
# setwd(dat_dir)

# Read parameter values from input file ----------------------------------------
read_inits = function(input_file = "input.par"){
#
# Read parameter values from input file.
# Returns a list, with parameters as elements.
#
#=======================================# Initialize vectors with stock specific values
  surv_freq = rep(x = 0, times = 2)     # Interval between availability of new abundance estimate
  KK = rep(x = 0, times = 2)            # Carrying capacity
  cv_n = rep(x = 0, times = 2)          # Coefficient of variation (CV) for abundance estimates
  cv_mortality = rep(x = 0, times = 2)  # CV for mortality (mortality is stochastic)
  F_r = rep(x = 0, times = 2)           # Recovery factor
  init_depl = rep(x = 0, times = 2)     # Depletion in first year (initial abundance / carrying capacity)
#====================================== # These vectors not read from input file, but it's convienient to initialize them here.
  sum_NPR_age = rep(x = 0, times = 2)   # Total numbers (summed across ages) per female recruit
  NPR_oneplus = rep(x = 0, times = 2)   # Total numbers (summed across ages aged 1+) per female recruit (note R indexing starts at 1, not zero, so this is sum(2:age_x))
  sum_NPR = rep(x = 0, times = 2)       #
  NPR_age = matrix(0, nrow = 2, ncol = 2)    # Numbers at age per female recruit (row) by stock (column)
  Nage_imm_0 = matrix(0, nrow = 2, ncol = 2) #
  Nage_mat_0 = matrix(0, nrow = 2, ncol = 2) # Vector of numbers of mature at each age
  prop_NPR = matrix(0, nrow = 2, ncol = 2)   # Vector with proportions in each age class for NPR, where sum(prop_NPR) = 1.0
  b_1 = rep(x = 0, times = 2)                # Birth rate in first year (given initial depletion)
  names(pars)
#====================================== # Read input to parameter data.frame.
  par_df = read.table(file = input_file, skip = 4, stringsAsFactors = FALSE)
  names(par_df) = c("Par", "Val")
  ref = strsplit(readLines(input_file, n = 2)[2], ":")[[1]][1]
#====================================== # Process input parameters
  cseed = par_df[1,2]
  iseed = as.integer(par_df[2,2])
  n_sims = as.integer(par_df[3,2])
  n_stocks = as.integer(par_df[4,2])
  yr_max = as.integer(par_df[5,2])

  surv_freq[1] = as.integer(par_df[6,2])
  # surv_freq[2] = as.integer(par_df[7,2])

  KK[1] = as.numeric(par_df[7,2])
  KK[2] = as.numeric(par_df[8,2])

  cv_n[1] = as.numeric(par_df[9,2])
  # cv_n[2] = as.numeric(par_df[10,2])

  cv_mortality[1] = as.numeric(par_df[10,2])
  cv_mortality[2] = as.numeric(par_df[11,2])

  theta = as.numeric(par_df[12,2])
  R_max = as.numeric(par_df[13,2])

  F_r[1] = as.numeric(par_df[14,2])
  F_r[2] = as.numeric(par_df[15,2])

  init_depl[1] = as.numeric(par_df[16,2])
  init_depl[2] = as.numeric(par_df[17,2])

  lower_tail = as.numeric(par_df[18,2])
  B_max = as.numeric(par_df[19,2])
  B_sex_ratio = as.numeric(par_df[20,2])
  S_adult = as.numeric(par_df[21,2])
  S_juv = as.numeric(par_df[22,2])
  a_t = as.integer(par_df[23,2]) + 1   # Need to add one, as booking step, to be comparible with Fortran code,
  a_m = as.integer(par_df[24,2]) + 1   #  because indexing has been defined to start at zero in Fortran.
  age_x = as.integer(par_df[25,2]) + 1
  a_r = as.integer(par_df[26,2]) + 1

  p_a1_s1 = as.numeric(par_df[27,2])
  p_a2_s1 = as.numeric(par_df[28,2])
  p_a2_s2 = as.numeric(par_df[29,2])
  p_a3_s2 = as.numeric(par_df[30,2])
  p_a4_s2 = as.numeric(par_df[31,2])

# TODO : Update list with omega parameters
# TODO : Update list with bias parameters e.g., m_bias, n_bias, r_bias etc.

  par_ls = list(ref = ref, cseed = cseed, iseed = iseed, n_sims = n_sims, n_stocks = n_stocks,
                yr_max = yr_max, surv_freq = surv_freq, KK = KK,
                cv_n = cv_n, cv_mortality = cv_mortality, theta = theta,
                R_max = R_max, F_r = F_r, init_depl = init_depl,
                lower_tail = lower_tail, b_max = B_max, b_sex_ratio = B_sex_ratio,
                S_adult = S_adult, S_juv = S_juv, a_t = a_t, a_m = a_m,
                age_x = age_x, a_r = a_r, p_a1_s1 = p_a1_s1, p_a2_s1 = p_a2_s1,
                p_a2_s2 = p_a2_s2, p_a3_s2 = p_a3_s2, p_a4_s2 = p_a4_s2,
                sum_NPR_age = sum_NPR_age, NPR_oneplus = NPR_oneplus, sum_NPR = sum_NPR)

  return(par_ls)
}

# Write parameter value to input file ------------------------------------------
write_inits = function(par_name, par_val, infile, outfile){
#
# Update value of `par_name` to `par_val`
# This function can be used to batch simulation runs,
#  e.g., over a sequence of values for lower tail / alternative N_min
#
#
  input_text = NULL; ii = NULL; current_val = NULL
  input_text = readLines(infile) # Read input.par file into a vector
  ii = grep(pattern = paste("^", par_name, sep = ""), x = input_text) # line number for this parameter
  # current_val = substr(input_text[ii], 17, 20)
  print(infile)
  print(input_text[ii])
  substr(input_text[ii], start = 17, stop = 24) = "          " # clean slate
  substr(input_text[ii], start = 17, stop = 17 + nchar(par_val)) = par_val
  # substr(input_text[ii], 17, 21) = par_val
  print(input_text[ii])
  # formatC(par_val, format = "f", digits = 3) # formatting needed or doesn't overwrite well
  writeLines(text = input_text, con = outfile)
}

# Check for files --------------------------------------------------------------
file_check = function(file_names){
#
# Check that all files in file_names exist in this directory
#  Used as a check when calling functions to run a batch list of file names.
#
  status = NULL
  file_check = file.exists(file_names)
  if (FALSE %in% file_check){
    missing = which(file_check == FALSE)
    stop(paste(file_names[missing], ": Files were not found. Stopping"))
    status = FALSE
  } else {
    print("All files exist")
    status = TRUE
  }
  return(status)
}

# Create file list -------------------------------------------------------------
create_file_list = function(n_trials, prefix){
#
# Create a list of file names for batching runs
#
  file_names = paste(prefix, n_trials, "A.txt", sep = "")
  file_names = c(file_names, paste(prefix, n_trials, "B.txt", sep = ""))
  file_names = c(file_names, paste(prefix, n_trials, "C.txt", sep = ""))
  file_names = c(file_names, paste(prefix, n_trials, "D.txt", sep = ""))
  return(file_names)
}

# Return par list --------------------------------------------------------------
pars = NULL
pars = read_inits()

