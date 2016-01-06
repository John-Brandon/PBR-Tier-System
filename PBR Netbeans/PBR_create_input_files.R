#----------------------------------------------------------------------------- #
# Automate the creation of several input files for PBR Tier System simulations
# Author   : John R. Brandon
# eMail    : jbrandon@gmail.com
# Date     : Fall 2015
# OS       : Mac OS 10.9.5
# R version: 3.2.0 (2015-04-16)
#
# Notes:
#  1. `create_file_list` & `file_check` functions defined in PBR_FileIO.R
#  2. File name convention = 'Txa_XY';
#    where: Txa = "Cet" or "Pin" for cetacean or pinniped
#      and: X = trial number (Wade 1998; Table 2)
#           Y = A (CV_N = 0.20 & F_R = 1.0)
#             = B (CV_N = 0.80 & F_R = 1.0)
#             = C (CV_N = 0.20 & F_R = 0.5)
#             = D (CV_N = 0.80 & F_R = 0.5)
#  3. For example, Cet_2D = Cetacean trial #2. For that trial:
#      Biased abundance estimates ->
#      N{hat} = 2.0 * N{true}
#      with CV_N = 0.8 and F_R = 0.50
#----------------------------------------------------------------------------- #

# Fetch parameters from input.par file -----------------------------------------
# source(file = "PBR_FileIO.R") # Functions for FileIO; also reads initial input.par

# Create vectors of file names -------------------------------------------------
cet_file_names = create_file_list(prefix = "Cet_", n_trials = 0:8); cet_file_names
pin_file_names = create_file_list(prefix = "Pin_", n_trials = 0:8); pin_file_names
all_file_names = c(cet_file_names, pin_file_names)

# File names for trials with CV_N = 0.20 and F_R = 1.0 ----------------------- #
a_file_names = paste("Cet_", 0:8, "A.txt", sep = "")
a_file_names = c(a_file_names, paste("Pin_", 0:8, "A.txt", sep = ""))

# File names for trials with CV_N = 0.80 and F_R = 1.0
b_file_names = paste("Cet_", 0:8, "B.txt", sep = "")
b_file_names = c(b_file_names, paste("Pin_", 0:8, "B.txt", sep = ""))

# File names for trials with CV_N = 0.20 and F_R = 0.5
c_file_names = paste("Cet_", 0:8, "C.txt", sep = "")
c_file_names = c(c_file_names, paste("Pin_", 0:8, "C.txt", sep = ""))

# File names for trials with CV_N = 0.80 and F_R = 0.5
d_file_names = paste("Cet_", 0:8, "D.txt", sep = "")
d_file_names = c(d_file_names, paste("Pin_", 0:8, "D.txt", sep = ""))

# Create files -----------------------------------------------------------------
# Boiler-plate input.par file available via GitHub repository
file.copy(from = "input.par", to = all_file_names, overwrite = TRUE)

# Edit generic parameters ------------------------------------------------------
# Set N_min percentile = 0.20 (other routines may alter this parameter in input.par)
# for (ii in 1:length(all_file_names)){
#   write_inits(par_name = "lower_tail", par_val = formatC(0.20, format = "f", digits = 3),
#               infile = all_file_names[ii], outfile = all_file_names[ii])
# }
# Make default r_max for pinnipeds = 0.12
for(ii in 1:length(pin_file_names)){
  write_inits(par_name = "r_max", par_val = formatC(0.12, format = "f", digits = 3),
              infile = pin_file_names[ii], outfile = pin_file_names[ii])
}

# Assume adult survival = 0.95 and max birth rate 1 yr for pinnipeds
for(ii in 1:length(pin_file_names)){
  write_inits(par_name = "S_adult", par_val = formatC(0.95, format = "f", digits = 3),
              infile = pin_file_names[ii], outfile = pin_file_names[ii])
  write_inits(par_name = "b_max", par_val = formatC(1.00, format = "f", digits = 3),
              infile = pin_file_names[ii], outfile = pin_file_names[ii])
}

# Set F_R = 0.5 for stock 1 in "C" trials
for(ii in 1:length(c_file_names)){
  write_inits(par_name = "F_r\\(1\\)", par_val = formatC(0.5, format = "f", digits = 3),
              infile = c_file_names[ii], outfile = c_file_names[ii])
}

# Set F_R = 0.5 for stock 1 in "D" trials
for(ii in 1:length(d_file_names)){
  write_inits(par_name = "F_r\\(1\\)", par_val = formatC(0.5, format = "f", digits = 3),
              infile = d_file_names[ii], outfile = d_file_names[ii])
}

# Make CV_N and CV_N_TRUE = 0.80 for "B" trials
for(ii in 1:length(b_file_names)){
  write_inits(par_name = "cv_n ", par_val = formatC(0.8, format = "f", digits = 3),
              infile = b_file_names[ii], outfile = b_file_names[ii])
  write_inits(par_name = "cv_n_true", par_val = formatC(0.8, format = "f", digits = 3),
              infile = b_file_names[ii], outfile = b_file_names[ii])
}

# Make CV_N and CV_N_TRUE = 0.80 for "D" trials
for(ii in 1:length(d_file_names)){
  write_inits(par_name = "cv_n ", par_val = formatC(0.8, format = "f", digits = 3),
              infile = d_file_names[ii], outfile = d_file_names[ii])
  write_inits(par_name = "cv_n_true", par_val = formatC(0.8, format = "f", digits = 3),
              infile = d_file_names[ii], outfile = d_file_names[ii])
}

### Edit trial parameters ------------------------------------------------------
# Trial 0 is base case and should be set through code above
t_n = substr(all_file_names, start = 5, stop = 5)
t_n = as.numeric(t_n)
t_1 = which(t_n == 1)
t_2 = which(t_n == 2)
t_3 = which(t_n == 3)
t_4 = which(t_n == 4)
t_5 = which(t_n == 5)
t_6 = which(t_n == 6)
t_7 = which(t_n == 7)
t_8 = which(t_n == 8)

# Create Trial #1 --------------------------------------------------------------
# Estimated mortality equal to one-half the actual mortality
for (ii in 1:length(t_1)){
  write_inits(par_name = "m_bias", par_val = formatC(2.0, format = "f", digits = 3),
              infile = all_file_names[t_1[ii]], outfile = all_file_names[t_1[ii]])
}

# Create Trial #2 --------------------------------------------------------------
# Estimated N twice actual N
for (ii in 1:length(t_2)){
  write_inits(par_name = "n_bias", par_val = formatC(2.0, format = "f", digits = 3),
              infile = all_file_names[t_2[ii]], outfile = all_file_names[t_2[ii]])
}

# Create Trial #3 --------------------------------------------------------------
# Assummed R_max is twice true R_max (true R_max = 0.04, assumed = 0.08)
for (ii in 1:length(t_3)){
  write_inits(par_name = "r_bias", par_val = formatC(0.5, format = "f", digits = 3),
              infile = all_file_names[t_3[ii]], outfile = all_file_names[t_3[ii]])
}

# This commented code would would aim to recreate Wade's (1998) methods
# cr_ii = substr(all_file_names[t_3], 1, 3) %in% "Cet" # get index for cetacean input files for trial 3
# cr_ii = t_3[which(cr_ii)]
# pr_ii = substr(all_file_names[t_3], 1, 3) %in% "Pin" # get index for pinniped input files for trial 3
# pr_ii = t_3[which(pr_ii)]
#
# for(ii in 1:length(cr_ii)){
#   write_inits(par_name = "r_max", par_val = formatC(0.08, format = "f", digits = 3),
#               infile = all_file_names[cr_ii[ii]], outfile = all_file_names[cr_ii[ii]])
# }
#
# for(ii in 1:length(pr_ii)){
#   write_inits(par_name = "r_max", par_val = formatC(0.24, format = "f", digits = 3),
#               infile = all_file_names[pr_ii[ii]], outfile = all_file_names[pr_ii[ii]])
# }

# Create Trial #4 --------------------------------------------------------------
# CV_N < Actual CV_N (cv_n_true)
for (ii in 1:length(t_4)){
  t_l = substr(all_file_names[t_4[ii]], start = 6, stop = 6) # Trial letter
  if(t_l == "A" | t_l == "C"){
    write_inits(par_name = "cv_n_true", par_val = formatC(0.8, format = "f", digits = 3),
                infile = all_file_names[t_4[ii]], outfile = all_file_names[t_4[ii]])
  }else{
    write_inits(par_name = "cv_n_true", par_val = formatC(1.6, format = "f", digits = 3),
                infile = all_file_names[t_4[ii]], outfile = all_file_names[t_4[ii]])
  }
}

# Create Trial #5 --------------------------------------------------------------
# Estimated CV_M = 0.25 of actual CV_M (cv_m_true = 1.2)
for (ii in 1:length(t_5)){
  write_inits(par_name = "cv_mortality\\(1\\)", par_val = formatC(1.2, format = "f", digits = 3),
              infile = all_file_names[t_5[ii]], outfile = all_file_names[t_5[ii]])
  write_inits(par_name = "cv_mortality\\(2\\)", par_val = formatC(1.2, format = "f", digits = 3),
              infile = all_file_names[t_5[ii]], outfile = all_file_names[t_5[ii]])
}

# Create Trial #6 --------------------------------------------------------------
# Survey frequency every 8 yrs
for (ii in 1:length(t_6)){
  write_inits(par_name = "surv_freq", par_val = formatC(8, format = "d"),
              infile = all_file_names[t_6[ii]], outfile = all_file_names[t_6[ii]])
}

# Create Trial #7 --------------------------------------------------------------
# True MNPL equal to 0.45 * K (theta = 0.53)
for (ii in 1:length(t_7)){
  write_inits(par_name = "theta", par_val = formatC(0.53, format = "f", digits = 3),
              infile = all_file_names[t_7[ii]], outfile = all_file_names[t_7[ii]])
}

# Create Trial #8 --------------------------------------------------------------
# Mortality bias as in trial 1, with true MNP equal to 0.70 * K (theta = 5.04)
for (ii in 1:length(t_8)){
  write_inits(par_name = "theta", par_val = formatC(5.04, format = "f", digits = 3),
              infile = all_file_names[t_8[ii]], outfile = all_file_names[t_8[ii]])
  write_inits(par_name = "m_bias", par_val = formatC(2.0, format = "f", digits = 3),
              infile = all_file_names[t_8[ii]], outfile = all_file_names[t_8[ii]])
}

# Do headers -------------------------------------------------------------------
file_headers_cet = readLines(con = "file_headers_cet.txt") # Cetacean trials
file_headers_pin = readLines(con = "file_headers_pin.txt") # Pinniped trials

# Write headers to cetacean input files across trials
for (ii in seq_along(cet_file_names)){
  t_id = substr(cet_file_names[ii], start = 1, stop = 6)
  tmp_text = readLines(con = cet_file_names[ii])
  header_line = grep(pattern = paste(t_id, sep=""), x = file_headers_cet)
  tmp_text[1] = file_headers_cet[header_line]
  tmp_text[2] = file_headers_cet[header_line + 1]
  writeLines(tmp_text, con = cet_file_names[ii])
}

# Write headers to pinniped input files across trials
for (ii in seq_along(pin_file_names)){
  t_id = substr(pin_file_names[ii], start = 1, stop = 6)
  tmp_text = readLines(con = pin_file_names[ii])
  header_line = grep(pattern = paste(t_id, sep=""), x = file_headers_pin)
  tmp_text[1] = file_headers_pin[header_line]
  tmp_text[2] = file_headers_pin[header_line + 1]
  writeLines(tmp_text, con = pin_file_names[ii])
}

# End --------------------------------------------------------------------------


