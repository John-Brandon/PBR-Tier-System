File descriptions: PBR-Tier-System
===============

## Fortran Files

1. **main.f90** 
  * This is the main program file. At the top of the file are a list of "use" statements. The use statements include the following modules of code: 

2. **Declare_variables_module.f90**  
  * Module containing variable declarations, which are potentially global in scope.

3. **PBR_FileIO_Module.f90**  
  * Module containing file input / output routines.

4. **PBR_Errorcheck_module.f90**  
  * Module for error checking input file. 

5. **Initialize_pop_module.f90**  
  * Module for initializing age structure.

6. **BRENT.f90**  
  * Code for finding a root using Brent's method (Press et al. Num Recipes).

7. **PBR_calcs_module.f90**  
  * Module containing routines for calculating PBR.

8. **A_Random_module.f90**  
  * Module that seeds random number generator (given user input).

9. **Generate_random_numbers_module.f90**  
  * Module with random number generating procedures.

10. **Eigen_module.f90**  
  * Module with wrapper functions for eigen analysis (calculating R_max).

11. **input.par**  
  * This file contains the list of input variables (e.g., survey frequency, etc.) that differ between trials. 

## R Files 
R is a freely available programming environment for statistical computing and graphics. Instructions and downloads for installing R are provided for Linux, Windows and Mac OS at <a href="https://cran.r-project.org/" target="_blank">the CRAN website</a>.

The R code assumes that the Fortran code has been compiled into an executable, named "main", residing in the same directory. Instructions for compiling the Fortran code are <a href="https://github.com/John-Brandon/PBR-Tier-System/blob/master/README.md" target="_blank">available here.</a>

[//]: # (At present, some modification would be needed to run the batch scripts under Windows. They involve `system()` calls, which differ between operating systems.)

The code in the **PBR_Fortran_controller.R** script will `source()` the functions in the other R files, and contains several calls to batch runs of simulations. It is therefore recommended, at least as a first pass, that the 'controller' R script be run using <a href="https://www.rstudio.com/" target="_blank">RStudio</a>, or another IDE that allows individual lines or segments of code to be run (segments in the controller script are delineated by commented lines). This will also provide a better sense of what the R code is doing, and  finer control over the batching process at first. In addition to batching runs, the controller script compiles depletion statisics and creates plots of from the simulation output.  

If you're running a Linux OS (e.g. Ubuntu), and don't have the pre-requisite R packages installed ("ggplot2", "stringr", "tidyr", "magrittr", and "dplyr"), you may need to start your R session using the `sudo` terminal command. This may also be necessary in general, if admin priviliges are needed to write to files on your OS. Further notes on R under Ubuntu are <a href="https://github.com/John-Brandon/PBR-Tier-System/blob/master/PBR%20Netbeans/Readme_Ubuntu_R.md" target="_blank">available here</a>. 

Further notes below:   

1. **PBR_FileIO.R**
  * The `write_inits()` function is relied on heavily for batching runs of simulations.
  
2. **PBR_create_input_files.R**
  * See the Notes in the header comments for the naming convention of files for each trial.

3. **PBR_batch.R**
  * Contains functions for running batches of simulations given lists of input file names. 

4. **PBR_fortran_output.R**
  * Contains functions for compiling depletion statistics and plotting. 

5. **PBR_Fortran_controller.R**
  * Runs scripts above, compiles depletion statistics and creates plots.
  