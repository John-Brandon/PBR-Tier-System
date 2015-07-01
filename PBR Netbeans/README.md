PBR-Tier-System
===============

# Brief file descriptions:
* `main.f90` This is the main program file. At the top of the file are a list of "use" statements. The use statements include the following modules of code: 

  * `Declare_variables_module.f90` Module containing variable declarations, which are potentially global in scope.

  * `PBR_FileIO_Module.f90` Module containing file input / output routines

  * `PBR_Errorcheck_module.f90` Module for error checking input file 

  * `Initialize_pop_module.f90` Module for initializing age structure

  * `BRENT.f90` Code for finding a root using Brent's method (Press et al. Num Recipes)

  * `PBR_calcs_module.f90` Module containing routines for calculating PBR

  * `A_Random_module.f90` Module that seeds random number generator (given user input)

  * `Generate_random_numbers_module.f90` Module with random number generating procedures

  * `Eigen_module.f90` Module with wrapper functions for eigen analysis (calculating lambda_max).

* `input.par` This file contains the list of input variables for running the model. 

* `nbproject` Contains Makefiles created by the Netbeans IDE (warning: likely outdated).

* `TEST COMMIT` TESTING4





