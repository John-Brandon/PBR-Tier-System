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
  * Module with wrapper functions for eigen analysis (calculating lambda_max).

11. **input.par**  
  * This file contains the list of input variables (e.g., survey frequency, etc.) that differ between trials. 

## R Files 

These files are included as an example work-flow for running trials as part of a management strategy evaluation. These files were developed under Mac OS 10.9.5. Some modifications would be needed to run these scripts under Windows. They involve `system()` calls, which differ between operating systems. Further notes below:   

Typical `source()` order for running scripts would be:

1. **PBR_FileIO.R**
  * The `write_inits()` function is relied on heavily for batching runs of simulations.
  
2. **PBR_create_input_files.R**
  * See the Notes in the header comments for the naming convention of input files for each trial.
  * Please note the comment regarding directory structure next to the initialization of the `code_dir` variable. This could be improved, but for now, it will require modification by the user to match where the code resides in their personal directories. 

3. **PBR_batch.R**
  * Contains functions for running batches of simulations given lists of input file names. 
  * Note: These functions were developed to run under Mac OS X, and the `system()` calls would likely need to be modified under a different OS. 
  * Future versions of these functions could be made more cross-platform friendly by making the `system()` calls conditional on the value returned from, for example `.Platform$OS.type`

4. **PBR_fortran_output.R**
  * Contains functions for compiling depletion statistics and plotting. 

5. **PBR_Fortran_controller.R**
  * Runs scripts above, compiles depletion statistics and creates some plots.








