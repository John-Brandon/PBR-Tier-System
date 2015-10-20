program main        
!#######################################################    
! File:     main.f90
! Author:   John R. Brandon
! Contact:  jbrandon at gmail             
! Date:     Summer 2015
! Purpose:  Run PBR Tier System simulations 
! License:  The MIT License
!#######################################################              
! Modules of code contain: subroutines, functions and possibly variable / format declarations.
!  Each module of code is contained in a separate file (e.g. Declare_variables_module.f90).
!  The code in each module file needs to be compiled and linked with the compiled main program to produce an executable.
!  Note: The order in which these modules are compiled relative to the main program may matter. (using the gfortran compiler, at least) 
!   Assume that the modules need to be compiled before they are linked with the main program to form an executable. 
!
!   This can be done in one line from the shell prompt (denoted by $). Note the file order in the command-line example below:
!    $ gfortran module1.f90 module2.f90 main.f -o desired_executable_name_here
!
!   Whereas, this shell command might not link properly, because the main program is compiled before module2.f90:
!    $ gfortran module1.f90 main.f module2.f90 -o desired_executable_name_here
!             
!   This can be a source of maddening errors during compiling if you're developing code in an integrated development environment (IDE). 
!   As an example, using NetBeans IDE 8.0.2 under Mac OS 10.9.5, whether or not the first letter of the file name 
!     is upper or lower case splits a tie in alphabetic file order (e.g. "Module2.f90" is linked before "module1.f90").             
!    
!   If you can't get things to compile and link in an IDE (e.g. 'Error can't find .mod file') you can try and figure out the settings the IDE
!    uses when ordering its list of files to the Fortran compiler, and then change those settings to make sure the module files are compiled in proper order before linking. 
!    Alternatively, you can use the shell prompt to manually control the correct order of file compilation (as in the first $ example above).   
!  
!   Here's an example compiler command to try on Mac OS X (should produce a 'main' executable)
!    gfortran A_Random_module.f90 BRENT.f90 Declare_variables_module.f90 Eigen_module.f90 Generate_random_numbers_module.f90 Initialize_pop_module.f90 PBR_Errorcheck_module.f90 PBR_FileIO_Module.f90 PBR_calcs_module.f90 main.f90 -o main -fbounds-check -framework accelerate
!#######################################################
!   General comments:
!         (a) The code for the population dynamics model would need to be revised to take into account reproductive senescence 
!              (as might be expected for at least some "black-fish", e.g. killer whales).
!         (b) 
!         (c) 
!   (ii)             
!#######################################################
  use Declare_variables_module       ! Declares global variables accessible by the main program here: Declare_variables_module.f90
  use PBR_FileIO_Module              ! Reading initial values from files, and for writing output : PBR_FileIO_Module.f90
  use initialize_pop                 ! Initialization of life history and age structure : Initialize_pop_module.f90
  use calcs                          ! Routines for various calculations (e.g. calculating N_min) : PBRmodule.f
  use random, only : random_normal, qnorm   ! Routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : Random_module.f90
  use Generate_random_numbers_module ! Determine if seed for RNG is user defined (for reproducible results) or if seed is based on CPU clock (for different psuedo random variates each time program runs): Generate_random_numbers_module.f90
  use PBR_Errorcheck_module          ! Contains function 'error_check_input' to do error checking on input values [Beta]
  use eigen_module                   ! Contains calls to DGEEV for calculating the eigenvalues and eigenvectors of the projection matrix
!---> Turns off implicit typing by Fortran; all variables must be explicitly declared by type
  implicit none 
!---> Constant parameters
  integer(kind = 4), parameter :: stock_1 = 1         ! Indexes for stock structure array
  integer(kind = 4), parameter :: stock_2 = 2         
  integer(kind = 4), parameter :: all_areas = 0
  integer(kind = 4), parameter :: area_1 = 1
  integer(kind = 4), parameter :: area_2 = 2
  integer(kind = 4), parameter :: area_3 = 3    
  integer(kind = 4), parameter :: area_4 = 4 
  integer(kind = 4), parameter :: n_area = 4          ! Counter for do loops -- could instead read this as input from file?
  integer(kind = 4), parameter :: female = 1          
  integer(kind = 4), parameter :: male = 2            
!---> Local variables                 
  real(kind = 8), allocatable :: f_init_ii(:)         ! Initial human caused mortality rates for each stock  
  real(kind = 8), allocatable :: f_yr_stock(:,:)      ! Human caused mortality rate each year (rows) by stock (columns)
  real(kind = 8), allocatable :: b_init_ii(:)         ! Initial birth rate for each stock 
  real(kind = 8), allocatable :: b_yr_stock(:,:)      ! Birth rate each year (rows) by stock (columns)    
  real(kind = 8), allocatable :: b_yr_stock_sim(:,:,:)      ! Birth rate each year by stock and simulation
  real(kind = 8), allocatable :: depl_yr_stock_sim(:,:,:)   ! Depletion each year (rows) by stock (columns)
  real(kind = 8), allocatable :: depl_yr_stock(:,:)
  real(kind = 8), allocatable :: sorted_final_depl(:)
  real(kind = 8), allocatable :: transition_matrix_tmp(:, :) ! Tmp matrix to pass to eigen function (is changed by eigen on return)
  real(kind = 8), allocatable :: movement_matrix(:, :, :) ! rows = ages; cols = areas & stock. Values = proportion of stock in each area

! Array of numbers-at-sex and age for each stock in each sub-area by year of projection        
! Main pop array: age, sex, sub-area, stock, year, & simulation number 
  real(kind = 8), allocatable :: N_age_sex_area_stock_yr_sim(:, :, :, :, :, :)  
  real(kind = 8), allocatable :: N_tot_sex_area_stock_yr_sim(:, :, :, :, :)
  real(kind = 8), allocatable :: area_stock_prop(:,:)    ! Percentage of each stock in each area
!  real(kind = 8), allocatable :: eigv(:)                 ! Eigen vector of transition matrix (not currently used)
  integer(kind = 4), allocatable :: is_surv_yr( : )      ! Matrix containing 1s if element is survey year, zeros otherwise
  real(kind = 8), allocatable :: n_hat_yr_sim(:, :)      ! Estimate of abundance each year
  real(kind = 8), allocatable :: n_tier_yr_sim(:,:)      ! Estimate of abundance based on tier approach
  real(kind = 8), allocatable :: N_plus_area123(:, :, :) ! Total age 1+ abundance in the survey area, by stock and simulation.  
  real(kind = 8), allocatable :: N_tot_area123(:, :, :)  ! Total (age 0+) abundance in the survey area, by stock and simulation.
  real(kind = 8), allocatable :: N_plus_yr_stock_sim(:, :, :) ! Total age 1+ abundance by yr, stock, and simulation across all areas.  
  real(kind = 8), allocatable :: N_tot_yr_stock_sim(:, :, :) ! Total (age 0+) abundance by yr, stock, and simulation across all areas.  
  real(kind = 8), allocatable :: N_mature_yr_stock_sim(:, :, :)
  
! Variables associated with human caused mortality   
  real(kind = 8), allocatable :: pbr_yr_sim(:, :)        ! PBR by year and simulation
  real(kind = 8), allocatable :: sigma_pbr_yr_sim(:, :)  ! Standard deviation of human caused mortality estimate (given PBR and CV_Mortality)  
  real(kind = 8), allocatable :: M_yr_sim(:, :)          ! Realized human caused mortality each year ~N(mu = PBR, sigma = CV_mortality * PBR)
  real(kind = 8), allocatable :: M_age_sex_area_stock_yr_sim(:, :, :, :, :, :)  
  real(kind = 8), allocatable :: M_scale_yr_sim(:,:)     ! Quantity to scale (prorate) mortality across areas (denominator of spatial allocation equation)
  real(kind = 8), allocatable :: omega_yr_area(:,:)      ! Relative vulnerability of animals in area j to human caused mortality each yr
  real(kind = 8), allocatable :: selectivity_norm(:)     ! Standardized selectivity at age: selectivity(age) / sum(selectivity(age))

  real(kind = 8) :: cv_n_tier_out                        ! CV of N_hat returned from tier system calculations (e.g. weighted average)
  real(kind = 8) :: NPR_mature                           ! Numbers mature per female recruit. Used to calculate b_eq (birth rate at K)
  real(kind = 8) :: NPR_0_F                              ! Recruits (age 0) as function of F {Left side: Eqn 16 Punt 1999}
  real(kind = 8) :: NPR_sum_recd                         ! Total number of recruited animals per recruit
  real(kind = 8) :: NPR_sum_unrecd                       ! Total number of unrecruited animals per recruit
  real(kind = 8) :: objf_lambda                          ! Objective function for finding juvenile survival that results in R_max
  real(kind = 8) :: objf_f_init                          ! Objective function for finding f_init resulting in stable age structure
  real(kind = 8) :: brent                                ! Function brent() :: file = Brent.f90 
  real(kind = 8) :: lambda                               ! Dominant real eigen value of the transition matrix
  integer(kind = 4) :: io_error                          ! Error flag for checking initial values in input.par
  integer(kind = 4) :: sim_ii, ii, jj, aa, ss, yr        ! Counters for indexing loops  
  integer(kind = 4) :: stock_ii
  integer(kind = 4) :: lower_ii, upper_ii                ! index for 5th and 95th percentiles (elements of vector)
  integer(kind = 4) :: flag
  
  real(kind = 8), allocatable :: foo_vector(:)           ! DEBUGGING  
  real(kind = 8) :: foo, foo1, foo2                      ! DEBUGGING 
  real(kind = 8) :: test_sort(1:100) 
  
  character*4 :: percentile
  real(kind = 8) :: norm_deviate
  real(kind = 8) :: lower_depl, median_depl, upper_depl
  integer(kind = 4) :: ifault
!---> Read initial values and do error checking
  call read_inits()              ! Read initial values from 'input.par' file. PBR_FileIO_Module contains subroutine 'read_inits()' 
  io_error = 0                   ! Check for input value errors (values out of bounds, etc) 
  io_error = error_check_input() ! Function that contains error checking code. Located in file: PBR_Errorcheck_module.f90
  if(io_error .ne. 0) then       ! Check error code
    print *, "Error Code: ", io_error  
    stop                         ! If error, exit program
  end if                         ! TODO? Move this error checking to be called at the end of read_inits() procedure? 
!---> Given input, allocate array dimensions at run-time (local variables)
  allocate(f_init_ii(1:n_stocks))     ! Stock specific initial human caused mortality rates (can have stock specific init_depl's)
  allocate(f_yr_stock(0:yr_max, 1:n_stocks)) ! Human caused mortality rate each year (rows) by stock (columns)
  allocate(b_init_ii(1:n_stocks))        ! Stock specific initial birth rate, given stock specific initial depletion levels
  allocate(b_yr_stock(0:yr_max, 1:n_stocks)) ! Birth rate each year (rows) by stock (columns)
  allocate(b_yr_stock_sim(0:yr_max, 1:n_stocks, 1:n_sims))
  allocate(depl_yr_stock_sim(0:yr_max, 0:n_stocks, 0:n_sims)) ! Depletion each year (rows) by stock (columns)
  allocate(depl_yr_stock(0:yr_max, 0:n_stocks))
  allocate(sorted_final_depl(1:n_sims))
  allocate(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 0:4, 1:n_stocks, 0:yr_max, 0:n_sims)) ! Main pop array, sim = 0 is reference case with no mortality    
  allocate(N_tot_sex_area_stock_yr_sim(1:2, 1:n_area, 1:n_stocks, 0:yr_max, 0:n_sims)) ! 
  allocate(area_stock_prop(1:4, 1:n_stocks))   ! TODO: Hard-coded four areas here - TODO: soft-code
  allocate(transition_matrix(0:age_x, 0:age_x))! Transition matrix, added here while developing methods for population projections and solving for juvenile survival (given lamda_max)
  allocate(transition_matrix_tmp(0:age_x, 0:age_x)) ! Place-holder for transition matrix. Gets sent to eigen and is changed on return
  allocate(movement_matrix(0:age_x, 0:n_area, 0:n_stocks))
  allocate(N_plus_area123(0:yr_max, 0:n_stocks, 0:n_sims)) ! Total age 1+ abundance in the survey area, by stock and simulation. Stock = 0 will hold the sum of both stocks.  
  allocate(N_tot_area123(0:yr_max, 0:n_stocks, 0:n_sims))  ! Total (age 0+) abundance in the survey area, by stock and simulation. Stock = 0 will hold the sum of both stocks.   
  allocate(N_plus_yr_stock_sim(0:yr_max, 0:n_stocks, 0:n_sims))
  allocate(N_tot_yr_stock_sim(0:yr_max, 0:n_stocks, 0:n_sims))
  allocate(N_mature_yr_stock_sim(0:yr_max, 0:n_stocks, 0:n_sims))
  
  allocate(pbr_yr_sim(0:yr_max, 0:n_sims))       ! PBR each year for each simulation
  allocate(M_yr_sim(0:yr_max, 0:n_sims))         ! Realized human caused mortality each year ~N(mu = PBR, sigma = CV_mortality * PBR)
  allocate(M_age_sex_area_stock_yr_sim(0:age_x, 1:2, 0:4, 1:n_stocks, 0:yr_max, 0:n_sims)) ! Mortality array 
  allocate(M_scale_yr_sim(0:yr_max, 0:n_sims))
  allocate(omega_yr_area(0:yr_max, 1:n_area))    ! Relative vulnerability of animals (independent of stock) by area.
  allocate(sigma_pbr_yr_sim(0:yr_max, 0:n_sims)) ! Standard deviation of mortality (currently constant through time -- this is redundant) 
  allocate(selectivity_norm(0:age_x))            ! Selectivity at age, standardized so that the vector sums to 1.0
!  allocate(eigv(0:age_x))                       ! Eigen vector of transition matrix, not currently used
  allocate(is_surv_yr(0:yr_max))                 ! Matrix containing 1s if element is survey year, zeros otherwise 
  allocate(n_hat_yr_sim(0:yr_max, 0:n_sims))     ! Abundance estimates
  allocate(n_tier_yr_sim(0:yr_max, 0:n_sims))! 
  allocate(N_min_yr_sim(0:yr_max, 0:n_sims))

  allocate(S_age(0:age_x))       ! Survival at age vector
  allocate(prop_mat_age(0:age_x))! Proportion mature at age
  allocate(selectivity(0:age_x)) ! Selectivity at age -- currently assuming knife-edge selectivity at age a_r (same selectivity for each stock)
  allocate(NPR_age(0:age_x))     ! Numbers at age per recruit vector
  allocate(NPR_age_tmp(0:age_x)) ! Numbers at age per recruit vector, used to solve for initial human caused mortality rate 

! These three below are assigned in calc_NPR_age(), but currently not utilized    
  allocate(prop_NPR(0:age_x))    ! Rescaled numbers at age per recruit vector (sums to 1.0 over ages)
  allocate(Nage_imm_0(0:age_x))  ! Immature numbers at age vector
  allocate(Nage_mat_0(0:age_x))  ! Mature numbers at age vector
  allocate(N_age(0:age_x))       ! Numbers at age vector    
  allocate(N_plus(0:yr_max))     ! Vector of age 1+ population size over projection years -- summed across stocks in survey area
  allocate(N_tot(0:yr_max))      ! Total (0+) population size each year of projection -- summed across stocks in survey area
  allocate(N_calf(0:yr_max))     ! Vector of calf production for each projection year -- summed across stocks in survey area
  allocate(Female_age(0:age_x))  ! Females at age vector -- summed across stocks in survey area
  allocate(Male_age(0:age_x))    ! Males at age vector -- summed across stocks in survey area
  allocate(foo_vector(0:age_x))  ! DEBUGGING
!   allocate(NPR_age_ii(0:age_x, 1:n_stocks))      ! Numbers at age per female recruit (row) for stock i (column)
!   allocate(NPR_oneplus_ii(1:n_stocks))           ! Total number of age 1+ per female recruit for stock i 

!---> Initialize those variables declared in Declare_variables_module -- set to zero.
  Call initialize_global_vars()         ! Does not initialize (i.e. overwrite) those variables with values read from input.par file

!---> Initialize those variables declared in main program (above) -- set to zero.
  call initialize_local_vars()          ! This subroutine is contained in the main program (at bottom)

!---> Set seed for RNG -- see comment after 'use Generate_random_numbers_module' statement above
  Call set_random_seed() ! Set seed based on input.par, either: (a) given # (reproducible results), or (b) based on CPU clock

!---> Initialize vectors: survival, selectivity and proportion mature at age   
  Call assign_par_vectors(a_r, a_m, a_t, age_x, s_adult, s_juv, &  
                          S_age, selectivity, prop_mat_age) ! Currently these three vectors are identical for each stock  

!---> Solve for juvenile survival rate that results in specified r_max
  transition_matrix(:,:) = 0.0d0        ! Need to initialize this matrix to zero, because evidently, if we don't, LAPACK is not nice.
  fecundity_max = b_max * b_sex_ratio   ! Define fecundity in terms of female calves per female for eigen analysis of matrix
! Assign non-zero values to transition matrix    
  transition_matrix = assign_transition_matrix(a_m, a_t, age_x, fecundity_max, S_age, prop_mat_age) 
  print *,
  print *, "Initial transition matrix (showing fecundity): "  ! Check
  do aa = 0, age_x        
      write (*, "(400f8.3)") (transition_matrix(aa, jj), jj = 0, age_x)    
! The 400f... is a bit of a hack. Works if <= 400 columns, age / stage classes to be printed 
  end do

!---> Take into account any bias in assumed R_max value, e.g. this where set value for Wade's (1998) trial 3  
  r_max = r_bias * r_max
  
!---> Test LAPACK procedure DGEEV for finding eigenvalues and eigenvectors of real nonsymetric matrix
  transition_matrix_tmp = transition_matrix ! Assign transition matrix to temp matrix, destroyed on return from eigen()
  print *, "Calling eigen() for initial transition matrix: "
  call eigen(transition_matrix_tmp, (age_x + 1), lambda) ! Calculate Lambda_max for input life history values: s_juv, etc.
  print *, "Initial lambda: ", lambda

!---> Solve for juvenile survival rate that corresponds with the user specified lambda_max
  objf_lambda = BRENT(ax = 0.01d0, bx = 0.98d0, cx = 0.99d0, func = calc_lambda, & ! See Brent.f90 for details on arguments
                          tol = 0.0000001d0, xmin = s_juv)                         ! Note calc_lambda() assigns projection matrix
  print *, "objf_lambda: ", objf_lambda                 ! Check, should be very close to "tol", i.e. nearly zero
  print *, "Solution for juvenile survival: ", s_juv    ! Check

!---> Initialize vectors: survival, selectivity and proportion mature at age   
  Call assign_par_vectors(a_r, a_m, a_t, age_x, s_adult, s_juv, &  ! Slightly repetitious, but need to assign s_juv solution to vector
                          S_age, selectivity, prop_mat_age)        ! Currently these three vectors are identical for each stock   

!---> Re-assign s_juv (etc) to transition matrix and output for checking -- assigning fecundity_max instead of b_max
  transition_matrix = assign_transition_matrix(a_m, a_t, age_x, b_max * b_sex_ratio, S_age, prop_mat_age) ! Assign non-zero values to transition matrix
  print *, "New transition matrix (showing max fecundity): "
  do aa = 0, age_x ! Print the transition matrix with lambda_max = 1.04 (Check above) to the screen
      write (*, "(400f8.3)") (transition_matrix(aa, jj), jj = 0, age_x) ! The 400f format is a hack. Works if <= 400 columns to be printed
  end do
  Call eigen(transition_matrix, (age_x + 1), lambda)
  print *, "Solution for lambda: ", lambda
  print *, "r_max true: ", r_max ! Check
  print *, ""

!---> Reset R_max for PBR calculations below
  r_max = r_max / r_bias
  print *, "Reset assumed R_max: ", r_max

! Return normal variate given lower tail (e.g. 1.96 for the 2.5%-tile) used when calculating N_min
  Call qnorm(p = lower_tail, normal_dev = norm_deviate, ifault = ifault) 
          
  if (ifault .ne. 0) then ! Check to make sure call to qnorm returns without error
    print *, "ERROR from `qnorm()` during simulation: ", sim_ii, " and year: ", yr
    stop  ! Stop program execution
  end if

! Use absolute value to be consistent with calc_n_min() method from Wade (1998) Eqn 4.          
  norm_deviate = ABS(norm_deviate) 
          
!---> Assign values to the matrix with percentage of each stock in each area 
  area_stock_prop = assign_area_stock_prop(p_a1_s1, p_a2_s1, p_a2_s2, p_a3_s2, p_a4_s2) 
    
!---> Calculate Numbers per female recruit, with no human caused mortality (f_init = 0.0)
  Call calc_NPR_age(f_rate = 0.0d0, &                        ! d0 suffix for double precision to match argument type in function                                             
      N_recruits = b_sex_ratio, N_age_tmp = NPR_age, &       ! Calc NPR_age, NPR_oneplus and NPR_mature (F = 0)
      sum_1plus = NPR_oneplus, sum_mature = NPR_mature, NPR_sum_recd = NPR_sum_recd, NPR_sum_unrecd = NPR_sum_unrecd) 
  print *, "Finished calling, calc_NPR_age(f_rate = 0.0): " 
  print *, "NPR_age (just females):" 
  print *, NPR_age
  print *, "NPR_age (females & males):"
  print *, NPR_age * 2.d0

! *  
! DEBUGGING - Comparing with Andre's `Yield Function.R` script - but this resets NPR_mature and hence b_eq incorrectly   
!  print *, ""
!  Call calc_NPR_age(f_rate = 0.01d0, &           ! d0 suffix for double precision to match argument type in function                                             
!      N_recruits = b_sex_ratio, N_age_tmp = NPR_age, & ! Calc NPR_age, NPR_oneplus and NPR_mature (F = 0)
!      sum_1plus = NPR_oneplus, sum_mature = NPR_mature, NPR_sum_recd = NPR_sum_recd, NPR_sum_unrecd = NPR_sum_unrecd) 
!  print *, "Finished call #2, calc_NPR_age(f_rate = 0.0): " 
!  print *, "NPR_age (just females):" 
!  print *, NPR_age
!  print *, "NPR_age (females & males):"
!  print *, NPR_age * 2.d0
!  print *, "MAXLOC(NPR_age): "
!  print *, MAXLOC(NPR_age)   ! This is like the `which(ii == max(x))` function in R, for vector x and index ii
! *
  
!---> Calculate pre-exploitation equilibrium birth rate (at carrying capacity when F = 0)
  b_eq = 1 / NPR_mature               ! Equilibrium birth rate. Equal for both stocks under assumption of identical life histories                
  print *, "Equilibrium birth rate"
  print *, b_eq

  ! * * * * * *
  ! DEBUGGGING
!  NPR_0_F = calc_recruit_F(NPR = NPR_age, b_eq = b_eq, selex_a = selectivity, surv_a = S_age)
               
  ! *
  
!---> Calculate the initial age structure and distribute across areas for each stock         
  do ii = 1, n_stocks    ! Initial age structures for each stock can differ, e.g. initial depletion may not be equal
! *
! TODO? : Move this loop into a subroutine in Initialize_pop_module 
! *        
    init_depl_i = init_depl(ii)     ! New value for global variable init_depl_i. Used by Initial_F(), as called from brent()
    b_init = b_eq + (b_max - b_eq) * (1 - (init_depl_i ** theta))  ! Initial birth rate for stock i
    b_init_ii(ii) = b_init                                         ! Storing initial birth rates in vector by stock 
    b_yr_stock(0, ii) = b_init          ! Store initial birth rate 
    depl_yr_stock(0, ii) = init_depl_i  ! Store initial depletion level        
    prop_NPR = NPR_age / sum(NPR_age)   ! Proportion at age relative to numbers per female recruit
    
!---> Calculate initial human caused mortality rate (f_init) that results in stable age-structure at initial depletion                
    objf_f_init = BRENT(ax = 0.0d0, bx = 0.10d0, cx = 1.0d0, func = initial_F, & ! See Brent.f90 for details on arguments
                        tol = 0.0000001d0, xmin = f_init_ii(ii)) 
    
    print *, ""                    
    print *, "f_init_ii(ii = 1)"
    print *, f_init_ii(1)
    print *, ""
    print *, "NPR_age_tmp(F = f_init)"
    print *, NPR_age_tmp
!    stop 
    
    Call rescale_NPR(k_1plus_tmp = k_1plus(ii), initial_oneplus_tmp = NPR_oneplus, & ! Scale NPR to Numbers at age
                    N_age_unscaled = NPR_age_tmp, N_age_scaled = N_age) ! N_age returned as scaled numbers of females at age 
    
! * 
! TODO : If initial conditions are the same for every simulation, can this be done more efficiently without a loop over sims? 
! In other words, do calculations once, then just set array values across sims equal for year 0. Can probably get rid of sim_ii loop. 
! * 
    do sim_ii = 0, n_sims ! Initial conditions are deterministic and hence identical across simulations. Make it so.
            
      N_age_sex_area_stock_yr_sim(:, female, 0, ii, 0, sim_ii) = N_age  ! Assign scaled numbers at age for this stock to main array
      N_age_sex_area_stock_yr_sim(:, male, 0, ii, 0, sim_ii) = N_age    ! Note: area = 0 represents the sum of numbers across all areas                        

      do jj = 1, n_area ! Allocate abundance across areas for each stock
        do ss = 1,2     ! males and females

          N_age_sex_area_stock_yr_sim( : , ss, jj, ii, 0, :) = &
              N_age_sex_area_stock_yr_sim( : , ss, all_areas, ii, 0, :) * area_stock_prop(jj, ii)

          N_tot_sex_area_stock_yr_sim(ss, jj, ii, 0, sim_ii) = & ! Sum total abundance by sex for each stock in each area
              sum(N_age_sex_area_stock_yr_sim(0:age_x, ss, jj, ii, 0, sim_ii))    

          end do  ! End loop over sexes
      end do      ! End loop over areas
      
      depl_yr_stock_sim(0, ii, sim_ii) = depl_yr_stock(0, ii) ! set initial depletion in year zero for each stock
      
    end do        ! End loop over simulations
! *
! TODO : Check can we set right hand side to sim_ii instead of zero ?? (Note the TODO in the main loop below as well)
! *
    N_plus_area123(0, ii, 0) = sum(N_age_sex_area_stock_yr_sim(1:age_x, 1:2, 1:3, ii, 0, 0)) ! N_plus index = (yr, stock, sim) 
    N_tot_area123(0, ii, 0) = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, ii, 0, 0))  ! Abundance in survey area for each stock 
    
    N_plus_yr_stock_sim(0, ii, 0) = sum(N_age_sex_area_stock_yr_sim(1:age_x, 1:2, 1:4, ii, 0, 0)) ! Sum stock size across areas
    N_tot_yr_stock_sim(0, ii, 0) = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:4, ii, 0, 0))  
    N_mature_yr_stock_sim(0, ii, 0) = sum(N_age_sex_area_stock_yr_sim(a_m:age_x, 1:2, 1:4, ii, 0, 0))
  end do          ! end loop over stocks
  
!---> Some final accounting in recording initial conditions - "N_plus_area123" numbers at age 1+ in survey area (1,2,3)
  N_plus_area123(0, 0, 0) = sum(N_age_sex_area_stock_yr_sim(1:age_x, 1:2, 1:3, 1:n_stocks, 0, 0)) ! yr, stock, sim
  N_tot_area123(0, 0, 0) = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, 1:n_stocks, 0, 0))  
  N_mature_yr_stock_sim(0, 0, 0) = sum(N_age_sex_area_stock_yr_sim(a_m:age_x, 1:2, 1:4, 1:n_stocks, 0, 0))  
!---> Assign relative vulnerabilities by area (inter-annual variability not implemented)
! *
! TODO : Extend to allow for inter-annual variability in relative vulnerabilities by area -- multiply omega by random variable to add noise 
! *  
  do jj = 1, n_area ! Areas 
    omega_yr_area(: , jj) = omega(jj)
  end do 

!====== +++ === === +++ === === +++ === 
!---- START SIMULATIONS ------------------------------------------------------>
!====== +++ === === +++ === === +++ ===
  flag = 0              ! dummy variable to indicate whether or not this is first simulation 
    
  do sim_ii = 0, n_sims ! Note the zero'th sim is the reference projection without any human caused mortality    
    
! Year zero :: Assign initial abundance from numbers at age for each simulation
! * TODO : Do this summation one time (instead of n_sims times) then assign remaining (1 to n_sims) initial conditions equal to first summation    
    N_plus_area123(0, 0, sim_ii) = sum(N_age_sex_area_stock_yr_sim(1:age_x, 1:2, 1:3, 1:n_stocks, 0, sim_ii)) ! yr, stock, sim
    N_tot_area123(0, 0, sim_ii) = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, 1:n_stocks, 0, sim_ii))  ! 
    N_mature_yr_stock_sim(0, 0, sim_ii) = sum(N_age_sex_area_stock_yr_sim(a_m:age_x, 1:2, 1:4, 1:n_stocks, 0, sim_ii))
      
    do stock_ii = 1, n_stocks 
      N_plus_area123(0, stock_ii, sim_ii) = sum(N_age_sex_area_stock_yr_sim(1:age_x, 1:2, 1:3, stock_ii, 0, sim_ii)) ! yr, stock, sim
      N_tot_area123(0, stock_ii, sim_ii) = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, stock_ii, 0, sim_ii))  ! get abundance in survey area for each stock 
      N_mature_yr_stock_sim(0, stock_ii, sim_ii) = &
        sum(N_age_sex_area_stock_yr_sim(a_m:age_x, 1:2, 1:4, stock_ii, 0, sim_ii))
    end do
    
    do yr = 1, yr_max     ! Years -- starting at year one, because year zero is in the books and pop has been initialized     
    
      do ii = 1, n_stocks ! Stocks            

! Calculate density dependent annual birth rate for each stock ----------------       
        if (depl_yr_stock(yr - 1, ii) .le. 0.d0 .or. depl_yr_stock(yr - 1, ii) .ge. 1.d0) then
          b_yr_stock(yr, ii) = 0.d0
        else
          b_yr_stock(yr, ii) = b_eq + (b_max - b_eq) * (1 - depl_yr_stock(yr - 1, ii)**theta) 
        end if
        
        do jj = 1, n_area ! Areas    
          
          do ss = 1, 2    ! Sexes
            
!---> Reproduction and natural mortality --------------------------------------
! Note this procedure works with the total abundance (i.e. area = 0) for each stock          
            Call pop_projection(f_rate = 0.0d0, b_rate = b_yr_stock(yr, ii), &           ! Project males and females separately,
              N_age_old = N_age_sex_area_stock_yr_sim(:, ss, 0, ii, yr - 1, sim_ii),  &  !  because could potentially have different selectivities, sex ratios at birth, etc. 
              N_age_new = N_age_sex_area_stock_yr_sim(:, ss, 0, ii, yr, sim_ii))         !  in a future version of operating model         
                      
!            do aa = 0, age_x ! Ages (Don't think actually need to loop over ages, have just inserted ":" into array instead
!! Anything needed inside a loop over ages?       

! Redistribute each stock (by age) across areas -------------------------------              
            N_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) = & 
                N_age_sex_area_stock_yr_sim(:, ss, 0, ii, yr, sim_ii) * area_stock_prop(jj, ii)

!            end do ! End loop over ages

! Sum total by sex for this area for this stock this year ---------------------
            N_tot_sex_area_stock_yr_sim(ss, jj, ii, yr, sim_ii) = &
                sum(N_age_sex_area_stock_yr_sim(0:age_x, ss, jj, ii, yr, sim_ii)) 
                
          end do ! End loop over sexes
        end do   ! End loop over areas         
      end do     ! End loop over stocks 
      
! Check if this is a survey year --------------------------------------------------------- 
! mod() is intrinsic function for the remainder. If zero, it's a survey year
      is_surv_yr(yr) = mod(yr - 1, surv_freq) 
      
! If this is a survey year, calculate PBR
      if (is_surv_yr(yr) == 0) then 

! Check that population is not extinct already        
        if (sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, 1:n_stocks, yr, sim_ii)) > 0.d0) then
          
! Generate survey estimate of abundance --------------------------------------------------
! Estimates generated with CV = cv_n_true. Where cv_true not necessarily equal to cv_n. 
! cv_n is the estimated sampling error used to calculate N_min.          
! cv_n_true > cv_n represents a biased underestimate of the true sampling error
! Function gen_survey_estimate() located in PBR_calcs_module.f90
! Assumes surveys apply to ages 0+ (cf. Wade 1998) 
 
          n_hat_yr_sim(yr, sim_ii) = gen_survey_estimate( &
            true_abundance = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, 1:n_stocks, yr, sim_ii)), &
            cv = cv_n_true)  

! Account for potential bias in abundance estimates -------------------------------------- 
! If n_bias = 1.0, then abundance estimates are unbiased. 
! n_bias > 1.0 for scenario(s) where abundance estimates positively biased. 
! n_bias specified in input.par file.            
          n_hat_yr_sim(yr, sim_ii) = n_bias * n_hat_yr_sim(yr, sim_ii) 

! Apply any tier specific methods (e.g. averaging) to calculate current best estimate of abundance      
          call calc_n_hat(tier = tier, &
                          yr = yr, &
                          n_hat_yr_in = n_hat_yr_sim(:, sim_ii), &
                          cv_n_in = cv_n, &
                          n_hat_out = n_tier_yr_sim(yr, sim_ii), &
                          cv_n_out = cv_n_tier_out)
       
! Calculate N_min -------------------------------------------------------------		  		  
          N_min_yr_sim(yr, sim_ii) = calc_n_min(n_hat = n_tier_yr_sim(yr, sim_ii), &
                                                cv = cv_n_tier_out, &
                                                z_score = norm_deviate)
! Calculate PBR ---------------------------------------------------------------		                                                
          if (determ_pbr .eq. "N") then
            pbr_yr_sim(yr, sim_ii) = N_min_yr_sim(yr, sim_ii) * 0.5d0 * r_max * F_r(1)
          else if (determ_pbr .eq. "Y") then
! DEBUGGING -- set PBR equal to number of mature animals of box stocks in survey area ----            
            pbr_yr_sim(yr, sim_ii) = N_mature_yr_stock_sim(yr, 0, sim_ii) * F_rate
          end if
          
! ****************************************************************************************          
        else ! If population extinct, set these to zero 
          n_hat_yr_sim(yr, sim_ii) = 0.d0
          N_min_yr_sim(yr, sim_ii) = 0.d0
          pbr_yr_sim(yr, sim_ii) = 0.d0
        end if

      else 
! If not a survey year, set PBR (expectation of human caused mortality) equal to previous year        
        pbr_yr_sim(yr, sim_ii) = pbr_yr_sim(yr - 1, sim_ii) 

      end if      

! Set human caused mortality based on stochastic normal random variable (Wade 1998) or deterministically equal to PBR
      if (determ_pbr .eq. "N") then ! Stochastic
! Given CV for uncertainty in human caused mortality (CV_M), calculate the SD of human caused mortality         
        sigma_pbr_yr_sim(yr, sim_ii) = cv_mortality(1) * pbr_yr_sim(yr, sim_ii) 

! Generate normal random deviate with PBR as expectation of human caused mortality -- follows approach of Wade (1998) p.9 step 4       
        M_yr_sim(yr, sim_ii) = random_normal(mean = real(pbr_yr_sim(yr, sim_ii),  4), & ! Convert from 8 to 4 bytes w real() : TODO modify random_normal() to 8 byte
          sd = real(sigma_pbr_yr_sim(yr, sim_ii), 4)) 
      else ! Deterministic
        M_yr_sim(yr, sim_ii) = pbr_yr_sim(yr, sim_ii) ! Realized human caused mortality exactly equal to PBR in this case
      end if
! Account for potential bias in mortality estimates. 
! If m_bias = 1.0, then human caused mortality estimates are unbiased. 
! If m_bias > 1.0 then scenario(s) where mortality estimates negatively biased, e.g. actual mortality is twice PBR.       
      M_yr_sim(yr, sim_ii) = m_bias * M_yr_sim(yr, sim_ii) ! m_bias read from input file (e.g. input.par)
      
! If populations reach low numbers, negative mortality (i.e. zombies) possible. Alert user.
      if (M_yr_sim(yr, sim_ii) < 0.0d0) then 
          print *, "Negative mortality!: ", "yr: ", yr, "sim_ii: ", sim_ii, "M_yr_sim(yr, sim_ii): ", M_yr_sim(yr, sim_ii)
          print *, "Resetting human cause mortality to zero."
          M_yr_sim(yr, sim_ii) = 0.0d0
      end if

!====== +++ === === +++ === === +++ ===          
! Allocate spatial age- and sex-structured human caused mortality --------------      
!====== +++ === === +++ === === +++ ===   
      if (sim_ii > 0) then        ! First simulation do not (calculate?) or subtract mortality
!  These loops needed for the summation in the denominator of spatial mortality allocation equation
        M_scale_yr_sim(yr, sim_ii) = 0.d0 ! Should be initialized already, but just to be sure...
        
        do ii = 1, n_stocks ! Stocks 
          do jj = 1, n_area ! Areas    
            do ss = 1, 2    ! Sexes

! Numerator from spatial allocation of mortality Eqn XX                                 
              M_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) = M_yr_sim(yr, sim_ii) * omega_yr_area(yr, jj) &
                * selectivity(:) * N_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) 

! Denominator from spatial allocation of mortality Eqn XX                 
              M_scale_yr_sim(yr, sim_ii) = M_scale_yr_sim(yr, sim_ii) + &
                sum(omega_yr_area(yr, jj) * selectivity(:) * N_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii)) 

            end do ! Sexes
          end do   ! Areas
        end do     ! Stocks

! Make another pass to (i) calculate mortality by stock, area, sex and age classes; and (ii) subtract spatial mortality by age/sex class        
        do ii = 1, n_stocks ! Stocks 
          do jj = 1, n_area ! Areas    
            do ss = 1, 2    ! Sexes
!---> Prevent division by zero with extinction (else get NaN)               
              if ( M_scale_yr_sim(yr, sim_ii) == 0.d0) then 
                M_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) = 0.d0
              else  
                M_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) = &
                  M_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) / M_scale_yr_sim(yr, sim_ii)
              end if  

! Subtract spatial mortality for this stock in this area by age and sex class
                N_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) = &
                  N_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii) - M_age_sex_area_stock_yr_sim(:, ss, jj, ii, yr, sim_ii)
                
! Check for negative abundance at age                  
                do aa = 0, age_x
                  if (N_age_sex_area_stock_yr_sim(aa, ss, jj, ii, yr, sim_ii) < 0.d0) then
                      
                    N_age_sex_area_stock_yr_sim(aa, ss, jj, ii, yr, sim_ii) = 0.d0

! *
! TODO: Allocate surplus mortality at age 
! *                    
                    ! Calculate extra mortality 
                    ! Reset abundance for this age-class to zero
                    ! Redistribute additional unaccounted for mortality to age-class that have positive abundance 
                    
                  end if ! Finished check on negative abundance
                end do   ! Ages
            end do       ! Sexes
          end do         ! Areas
        end do           ! Stocks      
      end if             ! Finished allocating kills for sim_ii > 0
   
! Total numbers at age and sex for each stock (i.e. sum across areas and assign totals to 'area 0'        
      do ii = 1, n_stocks    ! Stocks 
          do ss = 1, 2       ! Sexes
            do aa = 0, age_x ! Ages
              if (ii == 1) then             
                N_age_sex_area_stock_yr_sim(aa, ss, 0, ii, yr, sim_ii) = &
                  sum(N_age_sex_area_stock_yr_sim(aa, ss, 1:2, ii, yr, sim_ii))

! Set negative abundance at age equal to zero for stock 1                     
                if (N_age_sex_area_stock_yr_sim(aa, ss, 0, ii, yr, sim_ii) <= 0.d0) then
                  N_age_sex_area_stock_yr_sim(aa, ss, 0, ii, yr, sim_ii) = 0.d0
                end if

              else if (ii == 2) then 
                N_age_sex_area_stock_yr_sim(aa, ss, 0, ii, yr, sim_ii) = &
                  sum(N_age_sex_area_stock_yr_sim(aa, ss, 2:4, ii, yr, sim_ii))

! Set negative abundance at age equal to zero for stock 2 (if two-stock model)                     
                if (N_age_sex_area_stock_yr_sim(aa, ss, 0, ii, yr, sim_ii) <= 0.d0) then
                  N_age_sex_area_stock_yr_sim(aa, ss, 0, ii, yr, sim_ii) = 0.d0 
                end if  

              end if 
            end do   ! Ages
          end do     ! Sexes
                    
        N_plus_area123(yr, ii, sim_ii) = sum(N_age_sex_area_stock_yr_sim(1:age_x, 1:2, 1:3, ii, yr, sim_ii)) ! by stock 
        N_tot_area123(yr, ii, sim_ii) = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, ii, yr, sim_ii)) 
        N_mature_yr_stock_sim(yr, ii, sim_ii) = sum(N_age_sex_area_stock_yr_sim(a_m:age_x, 1:2, 1:4, ii, yr, sim_ii))
! *
! TODO : Could minimize summing below, e.g. replace with N_plus(yr, ii, sim_ii) = stock size summed over areas
! *     
        if (dd_component .eq. "1") then  
          if (ii == 1) then   
            depl_yr_stock(yr, ii) = sum(N_age_sex_area_stock_yr_sim(1:age_x, :, 1:2, ii, yr, sim_ii)) / k_1plus(ii) ! In terms of ages 1+
          else if (ii == 2) then
            depl_yr_stock(yr, ii) = sum(N_age_sex_area_stock_yr_sim(1:age_x, :, 2:4, ii, yr, sim_ii)) / k_1plus(ii) ! In terms of ages 1+ 
          end if
        else if (dd_component .eq. "M") then  
          if (ii == 1) then   
            depl_yr_stock(yr, ii) = sum(N_age_sex_area_stock_yr_sim(a_m:age_x, :, 1:2, ii, yr, sim_ii)) / k_1plus(ii) ! In terms of ages 1+
          else if (ii == 2) then
            depl_yr_stock(yr, ii) = sum(N_age_sex_area_stock_yr_sim(a_m:age_x, :, 2:4, ii, yr, sim_ii)) / k_1plus(ii) ! In terms of ages 1+ 
          end if
        end if

        
        depl_yr_stock_sim(yr, ii, sim_ii) = depl_yr_stock(yr, ii) ! keep track of annual depletion by simulation            
        
      end do       ! Stocks
      N_plus_area123(yr, 0, sim_ii) = sum(N_age_sex_area_stock_yr_sim(1:age_x, 1:2, 1:3, 1:n_stocks, yr, sim_ii)) ! yr, stock, sim
      N_tot_area123(yr, 0, sim_ii) = sum(N_age_sex_area_stock_yr_sim(0:age_x, 1:2, 1:3, 1:n_stocks, yr, sim_ii))  ! add abundance across stocks, this is abundance available to survey       
      
    end do           ! End loop over years
      
    print *, "Finished Simulation: ", sim_ii
    
  end do          ! End loop over number of simulations
!====== +++ === === +++ === === +++ ===    
  print *, "Finished simulations"
  print *, "Writing output to files"
!====== +++ === === +++ === === +++ ===

!  write(percentile, '(f4.3)') lower_tail ! transfer from type real to type character 
!  print *, "N_min percentile: ", percentile    ! check
  
! Test output to try and reduce explicit do loops on file output and see if can speed that process
!  print *, "N_plus_area123(: , 1, 1) "
!  print *, N_plus_area123(: , 1, 1)
  
!====== +++ === === +++ === === +++ ===          
! Write results to output files -----------------------------------------------
!====== +++ === === +++ === === +++ ===   
  open(unit = 6, file = "depl_by_trial.out")
  write(6, "(4(a10))") "trial", "lower", "median", "upper"
  write(6, "(a10, 3(f10.4))") REF, sorted_final_depl(lower_ii), median_depl, sorted_final_depl(upper_ii)
  close(unit = 6)
  
!  open(unit = 5, file = "N_aggregated" // percentile // ".out")
  open(unit = 5, file = "N_aggregated.out")  
  write(5, "(11(a15))") "sim", "yr", "stock", "N_tot_area123", "N_plus_area123", "n_hat_yr", &
    "Nmat_yr_ii_sim", "n_tier_yr_sim", "depl_yr_stock", "pbr_yr_sim", "M_yr_sim"
    
20  format(3(i15), 4(f15.4)) !
  do sim_ii = 0, n_sims
    do yr = 0, yr_max
      do ii = 0, n_stocks
        
        write(5, 20, advance = 'no') sim_ii, yr, ii, N_tot_area123(yr, ii, sim_ii), &
          N_plus_area123(yr, ii, sim_ii), n_hat_yr_sim(yr, sim_ii), N_mature_yr_stock_sim(yr, ii, sim_ii)
        
        write(5, "(f15.4)", advance = 'no') n_tier_yr_sim(yr, sim_ii)
        
        write(5, "(f15.4)", advance = 'no') depl_yr_stock_sim(yr, ii, sim_ii) 
        
        write(5, "(f15.4)", advance = 'no') pbr_yr_sim(yr, sim_ii) 
        
        write(5, "(f15.4)") M_yr_sim(yr, sim_ii)
        
      end do
    end do
  end do
  
! Deprecated output below (N_array files get very large, etc.) 
!====== +++ === === +++ === === +++ ===    
! Writing results of major array with spatial-temporal age structure to output file 
! Also writing results of summing over ages in major array (i.e. the slightly smaller array that is pertinent to PBR_calcs)    
!====== +++ === === +++ === === +++ === TODO: Move this into a function in the File IO module       
!  open(unit = 3, file = "N_array.out")
!  open(unit = 4, file = "N_tot_sex_area_stock_yr_sim.out")    
!  write(3, "(13(a15))") "sim", "yr", "stock", "age", "sex", "all_areas", "area_1", "area_2", "area_3", "area_4", &
!      "N_plus_area123", "N_tot_area123", "n_hat_yr"  
!  write(4, "(8(a15))") "sim", "yr", "stock", "sex", "area1", "area2", "area3", "area4"
!30  format(5(i15), 5(f15.4)) ! 
!40  format(4(i15), 5(f15.4))
!  do sim_ii = 0, n_sims
!    do yr = 0, yr_max
!      do ii = 1, n_stocks
!        do ss = 1, 2
!          do aa = 0, age_x ! Note the implicit do loop over areas in next line -- and advance = 'no' for suppressing new line 
!
!            write(3, 30,  advance='no') sim_ii, yr, ii, aa, ss, (N_age_sex_area_stock_yr_sim( aa , ss, jj , ii, yr, sim_ii), &
!                  jj = 0, n_area)
!
!            write(3, "(3(f15.4))") N_plus_area123(yr, ii, sim_ii), N_tot_area123(yr, ii, sim_ii), n_hat_yr_sim(yr, sim_ii)
!
!            write(4, 40) sim_ii, yr, ii, ss, (N_age_sex_area_stock_yr_sim( aa , ss, jj , ii, yr, sim_ii), &
!                  jj = 1, n_area)
!!           write(4, "(5(f15.4))") (N_tot_sex_area_stock_yr_sim(ss, jj, ii, 0, sim_ii), jj = 1, n_area) 
!
!          end do  ! End ouput for this age
!        end do ! End ouput for this sex
!      end do ! End output this stock
!    end do ! End ouput for this year
!  end do ! End output for this simulation
!  !write (1, "(A7,I3)") "hello", 10    ! example of conversion and concantination of character string in fortran    
!  close(unit = 3) ! close output file
!  close(unit = 4) 
  
  print *, "Trial REF: ", REF  ! Reference case 

!$         print *, "Compiled with -fopenmp"    ! This is a test for compiling with OpenMP (parallel processor directive <- !$)
! OpenMP not available for earlier versions of gfortran (including gfortran 4.2)
    
  print *, "Closing down"
  return     
! *
  contains
! ----------------------------------------------------------------------------------------      
  subroutine initialize_local_vars()
! ----------------------------------------------------------------------------------------          
    print *, "Hello from initialize_local_vars()"
    f_init_ii = 0.d0
    f_yr_stock = 0.d0
    b_init_ii = 0.d0
    b_yr_stock = 0.d0
!    depl_yr_stock = 0.d0
    depl_yr_stock_sim = 0.d0
    transition_matrix_tmp = 0.d0
    movement_matrix = 0.d0
    N_age_sex_area_stock_yr_sim = 0.d0
    area_stock_prop = 0.d0
    M_scale_yr_sim = 0.d0
!    eigv = 0.d0                        ! Eigen vector of transition matrix, not currently used
    n_hat_yr_sim = 0.0d0
    is_surv_yr = 0
    pbr_yr_sim = 0.d0
    NPR_mature = 0.d0 
    NPR_age = 0.0d0
    NPR_oneplus = 0.0d0

    foo_vector = 0.d0     
    foo = 0.d0 
    foo1 = 0.d0  
    
    objf_lambda = 0.d0 
    objf_f_init = 0.d0 
    lambda = 0.d0 
    io_error = 0
    ii = 0
    jj = 0
    aa = 0
    ss = 0
    yr = 0
    print *, "Goodbye from initialize_local_vars()"
        !
    return
  end subroutine initialize_local_vars
! *
! ----------------------------------------------------------------------------------------        
  subroutine sort(arr)
! ----------------------------------------------------------------------------------------          
  implicit none
  real(kind = 8), dimension(:), intent(inout) :: arr
! Sorts an array arr into ascending numerical order by Shell’s method (diminishing increment
! sort). arr is replaced on output by its sorted rearrangement.
! From Press et al. Numerical Recipes
  integer(kind = 4) :: i,j,inc,n
  real(kind = 8) :: v
  n = size(arr)
  inc = 1
! Determine the starting increment
  do 
    inc = 3 * inc +1
    if (inc > n) exit
  end do
! *
! Loop over the partial sorts.  
  do 
    inc = inc / 3
    do i = inc + 1, n ! Outer loop of straight insertion.
      v = arr(i)
      j=i
      do              ! Inner loop of straight insertion.
        if (arr(j - inc) <= v) exit
        arr(j) = arr(j - inc)
        j = j - inc
        if (j <= inc) exit
      end do
      arr(j) = v
    end do
    if (inc <= 1) exit
  end do
  
  end subroutine sort
! *
! ----------------------------------------------------------------------------------------        
  real(kind = 8) function median(x)
! ----------------------------------------------------------------------------------------        
! calculate median of a *sorted* vector of double precision real numbers
  implicit none
  integer(kind = 4) :: n  ! length of vector
  real(kind = 8), intent(in) :: x(:)
  
  n = size(x) 
  
  if (MOD(n, 2) == 0) then           ! interpolate if length of vector is odd
    median = (x(n / 2) + x(n / 2 + 1)) / 2.0d0
  else
    median = x(n / 2 + 1)
  end if
  
  end function median
! ----------------------------------------------------------------------------------------        
end program main

