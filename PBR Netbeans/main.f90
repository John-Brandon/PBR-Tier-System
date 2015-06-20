program main        
!#######################################################    
! File:     main.f
! Author:   John R. Brandon
! Contact:  jbrandon at gmail             
! Date:     Spring 2015
! Purpose:  Run PBR Tier System simulations            
!#######################################################              
! Modules of code contain: subroutines, functions and possibly variable / format declarations.
!  Each module of code is contained in a separate file (e.g. Declare_variables_module.f90).
!  The code in each module file needs to be compiled and linked with the compiled main program to produce an executable.
!  Note: The order in which these modules are compiled relative to the main program may matter. (using the gfortran compiler, at least) 
!   Assume that the modules need to be compiled before they are linked with the main program to form an executable. 
!
!   This can be done in one line from the shell prompt ($); but again, note the file order in the command-line example below:
!    $ gfortran module1.f90 module2.f90 main.f -o desired_executable_name_here
!
!   Whereas, this shell command might not link properly, because the main program is compiled before module2.f90:
!    $ gfortran module1.f90 main.f module2.f90 -o desired_executable_name_here
!             
!   This can be a source of maddening errors during compiling if you're developing code in an integrated development environment (IDE). 
!   As an example, using NetBeans IDE 8.0.2 under Mac OS 10.9.5, whether or not the first letter of the file name is upper or lower case splits 
!    a tie in alphabetic file order (e.g. "Module2.f90" is linked before "module1.f90").             
!   If you can't get things to compile and link in an IDE (e.g. 'Error can't find .mod file') you can try and figure out the settings the IDE
!    uses when ordering its list of files to the Fortran compiler, and then change those settings to make sure the module files are compiled in proper order before linking. 
!    Alternatively, you can use the shell prompt to manually control the correct order of file compilation (as in the first $ example above).              
!#######################################################
!   General comments:
!   (i) At present, the population dynamics model assumes this order for entering female adult-hood: 
!       (1st) Maturity (ovulation) on her a_m'th birthday -> (2nd) Adult mortality rate applied during her a_m'th year -> (3rd) First partuition (one year after her a_m'th birthday, if she survives that first year of adult-hood).
!         (a) a_m is the age at which females mature (males are assumed to be a non-limiting factor for reproduction / birth rates)
!         (b) The code for the population dynamics model would need to be revised to take into account reproductive senescence (as might be expected for at least some "black-fish", e.g. killer whales).
!         (c) Gestation length after reaching a_m is assumed to be one year. TODO: make sure birth rates bounded between 0.0 and 1.0 
!   (ii)             
!#######################################################
    use Declare_variables_module       ! Declares global variables accessible by the main program here: Declare_variables_module.f90
    use PBR_FileIO_Module              ! Reading initial values from files, and for writing output : PBR_FileIO_Module.f90
    use initialize_pop                 ! Initialization of life history and age structure : Initialize_pop_module.f90
    use calcs                          ! Routines for various calculations (e.g. calculating N_min) : PBRmodule.f
    use random, only : random_normal   ! Routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : Random_module.f90
    use Generate_random_numbers_module ! Determine if seed for RNG is user defined (for reproducible results) or if seed is based on CPU clock (for different psuedo random variates each time program runs): Generate_random_numbers_module.f90
    use PBR_Errorcheck_module          ! Contains function 'error_check_input' to do error checking on input values [Very Beta]
    use eigen_module                   ! Contains calls to DGEEV for calculating the eigenvalues and eigenvectors of the projection matrix
!====== +++ === === +++ === === +++ === ! Turns off implicit typing by Fortran; all variables must be explicitly declared by type
    implicit none 
!====== +++ === === +++ === === +++ ===                 ! Constant parameters
    integer(kind = 4), parameter :: stock_1 = 1         ! Indexes for stock structure array
    integer(kind = 4), parameter :: stock_2 = 2         ! TODO? Move parameters into Main_pars_module?
    integer(kind = 4), parameter :: all_areas = 0
    integer(kind = 4), parameter :: area_1 = 1
    integer(kind = 4), parameter :: area_2 = 2
    integer(kind = 4), parameter :: area_3 = 3    
    integer(kind = 4), parameter :: area_4 = 4 
    integer(kind = 4), parameter :: n_area = 4          ! Counter for do loops
    integer(kind = 4), parameter :: female = 1          
    integer(kind = 4), parameter :: male = 2            ! TODO : Move this list of variable declarations into a module (e.g. main_vars_module)
!====== +++ === === +++ === === +++ ===                 ! Local variables 
    real(kind = 8), allocatable :: f_init_ii(:)         ! Initial human caused mortality rates for each stock  
    real(kind = 8), allocatable :: f_yr_stock(:,:)      ! Human caused mortality rate each year (rows) by stock (columns)
    real(kind = 8), allocatable :: b_init_ii(:)         ! Initial birth rate for each stock 
    real(kind = 8), allocatable :: b_yr_stock(:,:)      ! Birth rate each year (rows) by stock (columns)    
    real(kind = 8), allocatable :: depl_yr_stock(:,:) ! Depletion each year (rows) by stock (columns)
    real(kind = 8), allocatable :: transition_matrix_tmp(:, :) ! Tmp matrix to pass to eigen function (is changed by eigen on return)
    real(kind = 8), allocatable :: movement_matrix(:, :, :) ! rows = ages x cols = areas x stock. Values = proportion of stock in each area
! Array of numbers-at-sex and age for each stock in each sub-area by year of projection        
    real(kind = 8), allocatable :: N_age_sex_area_stock_yr(:,:,:,:,:)  ! Main pop array: age, sex, sub-area, stock, yr
    real(kind = 8), allocatable :: area_stock_prop(:,:)    ! Percentage of each stock in each area
    real(kind = 8), allocatable :: eigv(:)                 ! Eigen vector of transition matrix, e.g. returned by power_method()
    integer(kind = 4), allocatable :: seq_yrs(:)           ! Used to create a consecutive sequence of year numbers, e.g. 1-100 
    integer(kind = 4), allocatable :: is_surv_yr(:, :)     ! Matrix containing 1s if element is survey year, zeros otherwise
    real(kind = 8), allocatable :: pbr_yr_stock(:, :)  ! PBR each year (rows) for each stock (columns)
    real(kind = 8), allocatable :: foo_vector(:)      ! DEBUGGING
    real(kind = 8) :: NPR_mature                        ! Numbers mature per female recruit. Used to calculate b_eq (birth rate at K)
!    real(kind = 8) :: sum_1plus_tmp                     ! DEBUGGING
    real(kind = 8) :: foo, foo1                         ! DEBUGGING
!    real(kind = 8) :: start, finish                     ! For timing / optimizing code
    real(kind = 8) :: objf_lambda                       ! Objective function for finding juvenile survival that results in R_max
    real(kind = 8) :: objf_f_init                       ! Objective function for finding f_init resulting in stable age structure
    real(kind = 8) :: brent    ! Function brent() :: file = Brent.f90 (TODO? : Create and add this to a Roots_and_Mins_module.f90)
    real(kind = 8) :: lambda                            ! Dominant real eigen value of the transition matrix
    integer(kind = 4) :: io_error                       ! Error flag for checking initial values in input.par
    integer(kind = 4) :: ii, jj, aa, ss, yr             ! Counters for indexing loops  
!    integer(kind = 4) :: it_num                         ! Iteration number returned from power_method() for calculating lambda
!    character(len = 10) :: stock_name                   ! For printing output tables
!====== +++ === === +++ === === +++ === ! Read initial values and do error checking
    call read_inits()              ! Read initial values from 'input.par' file. PBR_FileIO_Module contains subroutine 'read_inits()' 
    io_error = 0                   ! Check for input value errors (values out of bounds, etc) 
    io_error = error_check_input() ! Function that contains error checking code. Located in file: PBR_Errorcheck_module.f90
    if(io_error .ne. 0) then       ! Check error code
        print *, "Error Code: ", io_error  ! If error, print code
        stop                       ! If error, exit program
    end if                         ! TODO? Move this error checking to be called at the end of read_inits() procedure? 
!====== +++ === === +++ === === +++ === ! Initialize
!    age_x = a_m + 1                     ! Might eventually want to change this to take user input for age_x (?)                   
!====== +++ === === +++ === === +++ === ! Given input, allocate array dimensions at run-time (local variables)
    allocate(f_init_ii(1:n_stocks))     ! Stock specific initial human caused mortality rates (can have stock specific init_depl's)
    allocate(f_yr_stock(0:yr_max, 1:n_stocks)) ! Human caused mortality rate each year (rows) by stock (columns)
    allocate(b_init_ii(1:n_stocks))        ! Stock specific initial birth rate, given stock specific initial depletion levels
    allocate(b_yr_stock(0:yr_max, 1:n_stocks)) ! Birth rate each year (rows) by stock (columns)
    allocate(depl_yr_stock(0:yr_max, 1:n_stocks)) ! Depletion each year (rows) by stock (columns)
    allocate(N_age_sex_area_stock_yr(0:age_x, 1:2, 0:4, 1:n_stocks, 0:yr_max)) ! Main pop array
    allocate(area_stock_prop(1:4, 1:n_stocks))  ! TODO: Hard-coded four areas here - TODO: soft-code
    allocate(transition_matrix(0:age_x, 0:age_x))! Transition matrix, added here while developing methods for population projections and solving for juvenile survival (given lamda_max)
    allocate(transition_matrix_tmp(0:age_x, 0:age_x)) ! Place-holder for transition matrix. Gets sent to eigen and is changed on return
    allocate(movement_matrix(0:age_x, 0:n_area, 0:n_stocks))
    allocate(eigv(0:age_x))         ! Eigen vector of transition matrix, returned by power_method()  
!    allocate(seq_yrs(0:yr_max))                  ! Used to create a consecutive sequence of year numbers, e.g. 1-100 
    allocate(is_surv_yr(0:yr_max, 1:n_stocks))     ! Matrix containing 1s if element is survey year, zeros otherwise 
    allocate(pbr_yr_stock(0:yr_max, 1:n_stocks))  
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
    allocate(Nplus(0:yr_max))      ! Vector of age 1+ population size over projection years 
    allocate(Ntot(0:yr_max))       ! Total (0+) population size each year of projection 
    allocate(N_calf(0:yr_max))     ! Vector of calf production for each projection year
    allocate(Female_age(0:age_x))  ! Females at age vector
    allocate(Male_age(0:age_x))    ! Males at age vector
    allocate(foo_vector(0:age_x))  ! DEBUGGING
!   allocate(depl_yr_stock(0:yr_max, 1:n_stocks))  ! Depletion in year t (row) of stock i (column) 
!   allocate(NPR_age_ii(0:age_x, 1:n_stocks))      ! Numbers at age per female recruit (row) for stock i (column)
!   allocate(NPR_oneplus_ii(1:n_stocks))           ! Total number of age 1+ per female recruit for stock i 
!====== +++ === === +++ === === +++ === ! Initialize those variables declared in Declare_variables_module 
    Call initialize_global_vars()       ! Does not initialize ('rewrite') those variables with values read from input.par file
!====== +++ === === +++ === === +++ === ! Initialize those variables declared in main program (above) 
    call initialize_local_vars()        ! This subroutine is contained in the main program (at bottom)
!====== +++ === === +++ === === +++ === ! Set seed for RNG -- see comment after 'use Generate_random_numbers_module' statement above
    Call set_random_seed() ! Set seed based on input.par, either: (a) given # (reproducible results), or (b) based on CPU clock
!====== +++ === === +++ === === +++ === ! Initialize vectors: survival, selectivity and proportion mature at age   
    Call assign_par_vectors(a_r, a_m, a_t, age_x, s_adult, s_juv, &  
                            S_age, selectivity, prop_mat_age) ! Currently these three vectors are identical for each stock  
!!====== +++ === === +++ === === +++ === ! Solve for juvenile survival rate that results in specified r_max
    transition_matrix = 0.0d0       ! Need to initialize this matrix to zero, because evidently, if we don't, LAPACK(?) is not nice.
! Specifically, I (JRB) had horrible buggy trouble with Mac OS X's 'accelerate' framework (LAPACK library call to function DGEEV).    
    fecundity_max = b_max * b_sex_ratio ! Define fecundity in terms of female calves per female for eigen analysis of matrix
    transition_matrix = assign_transition_matrix(a_m, a_t, age_x, fecundity_max, S_age, prop_mat_age) ! Assign non-zero values to transition matrix
    print *,
    print *, "Initial transition matrix: "  ! Check
    do aa = 0, age_x
        write (*, "(100f8.3)") (transition_matrix(aa, jj), jj = 0, age_x)    ! The 100f... is a bit of a hack. Works if <= 100 columns to be printed
    end do
!====== +++ === === +++ === === +++ === ! Test LAPACK procedure DGEEV for finding eigenvalues and eigenvectors of real nonsymetric matrix
    transition_matrix_tmp = transition_matrix ! Assign transition matrix to temp matrix, destroyed on return from eigen()
!
    print *, "Calling eigen() for initial transition matrix: "
!    
    call eigen(transition_matrix_tmp, (age_x + 1), lambda) ! Calculate Lambda_max for input life history values: s_juv, etc.
    print *, "Initial lambda: ", lambda
!====== +++ === === +++ === === +++ === ! Solve for juvenile survival rate that corresponds with the user specified lambda_max
    objf_lambda = BRENT(ax = 0.01d0, bx = 0.98d0, cx = 0.99d0, func = calc_lambda, & ! See Brent.f90 for details on arguments
                            tol = 0.0000001d0, xmin = s_juv)                         ! Note calc_lambda() assigns projection matrix
!
    print *, "objf_lambda: ", objf_lambda ! Check, should be very close to "tol", i.e. nearly zero within machine precision
    print *, "Solution for juvenile survival: ", s_juv    ! Check
!====== +++ === === +++ === === +++ === ! Initialize vectors: survival, selectivity and proportion mature at age   
    Call assign_par_vectors(a_r, a_m, a_t, age_x, s_adult, s_juv, &  ! Slightly repetitious, but need to re-assign s_juv rates
                            S_age, selectivity, prop_mat_age)        ! Currently these three vectors are identical for each stock   
!====== +++ === === +++ === === +++ === ! Re-assign s_juv (etc) to transition matrix and output for checking
    transition_matrix = assign_transition_matrix(a_m, a_t, age_x, b_max, S_age, prop_mat_age) ! Assign non-zero values to transition matrix
    print *, "New transition matrix: "
    do aa = 0, age_x ! Print the transition matrix with lambda_max = 1.04 (Check above) to the screen
        write (*, "(100f8.3)") (transition_matrix(aa, jj), jj = 0, age_x)    ! The 100f format is a hack. Works if <= 100 columns to be printed
    end do
    print *, ""
!====== +++ === === +++ === === +++ === ! Initialize the matrix with percentage of each stock in each area 
    area_stock_prop = 0.0d0 ! Initialize the movement matrix, before assigning values in next line
    area_stock_prop = assign_area_stock_prop(p_a1_s1, p_a2_s1, p_a2_s2, p_a3_s2, p_a4_s2) 
!====== +++ === === +++ === === +++ === ! Calculate Numbers per female recruit, with no human caused mortality (f_init = 0.0)
    NPR_age = 0.0d0
    NPR_oneplus = 0.0d0
    NPR_mature = 0.0d0
    print *, "NPR_age: ", NPR_age
    print *, "b_sex_ratio: ", b_sex_ratio
    foo_vector = NPR_age
! THIS NEXT CALL IS CAUSING A SEGMENTATION FAULT WHEN n_stocks = 1, using NetBeans IDE - John    
    Call calc_NPR_age(f_rate = 0.0d0, &                        ! d0 suffix for double precision to match argument type in function                                             
        N_recruits = b_sex_ratio, N_age_tmp = NPR_age, &       ! Calc NPR_age, NPR_oneplus and NPR_mature (F = 0)
        sum_1plus = NPR_oneplus, sum_mature = NPR_mature) 
    print *, "Finished calling, calc_NPR_age(): " 
    print *, "NPR_age :", NPR_age
!====== +++ === === +++ === === +++ === ! Calculate equilibrium birth rate (at carrying capacity)
    b_eq = 1 / NPR_mature               ! Equilibrium birth rate. Equal for both stocks under assumption of identical life histories                
!====== +++ === === +++ === === +++ === ! Calculate the initial age structure for each stock         
! TODO : Move this loop into a subroutine in Initialize_pop_module 
    do ii = 1, n_stocks                 ! Initial age structures for each stock can differ, e.g. initial depletion may not be equal
        init_depl_i = init_depl(ii)     ! New value for global variable init_depl_i. Used by Initial_F(), as called from brent()
        b_init = b_eq + (b_max - b_eq) * (1 - (init_depl_i ** theta))  ! Initial birth rate for stock i
        b_init_ii(ii) = b_init                                         ! Storing initial birth rates in vector by stock 
        b_yr_stock(0, ii) = b_init          ! Store initial birth rate 
        depl_yr_stock(0, ii) = init_depl_i  ! Store initial depletion level 
        print *, 
        print *, "Stock: ", ii                        ! DEBUGGING  
!        print *, "init_depl_i", init_depl_i
!        print *, "b_init: ", b_init
        print *, "prop_NPR :", prop_NPR
        print *, "NPR_age :", NPR_age
        print *, "sum(NPR_age) :", sum(NPR_age)
        print *, "NPR_age / sum(NPR_age) :", NPR_age / sum(NPR_age)
        prop_NPR = NPR_age / sum(NPR_age)
        print *, prop_NPR
! Calculate initial human caused mortality rate (f_init) that results in stable age-structure at initial depletion                
        objf_f_init = BRENT(ax = 0.0d0, bx = 0.10d0, cx = 1.0d0, func = initial_F, & ! See Brent.f90 for details on arguments
                            tol = 0.0000001d0, xmin = f_init_ii(ii)) 
!        
!        print *, "f_init_ii: ", f_init_ii(ii)        ! DEBUGGING
!        print *, "NPR_age_tmp: ", NPR_age_tmp
        print *, "NPR_oneplus: ", NPR_oneplus        
        Call rescale_NPR(k_1plus_tmp = k_1plus(ii), initial_oneplus_tmp = NPR_oneplus, & ! Scale NPR to Numbers at age
                        N_age_unscaled = NPR_age_tmp, N_age_scaled = N_age) ! N_age returned as scaled numbers of females at age 
                        
!        print *, "N_age: ", N_age                
        print *, "Debugging"
        N_age_sex_area_stock_yr(:, female, 0, ii, 0) = N_age  ! Assign scaled numbers at age for this stock to main array
        N_age_sex_area_stock_yr(:, male, 0, ii, 0) = N_age    ! Note: area = 0 represents the sum of numbers across all areas

        do jj = 1, n_area                             ! Allocate stock abundance across areas
            N_age_sex_area_stock_yr( : , female, jj, ii, 0) = &
              N_age_sex_area_stock_yr( : , female, all_areas, ii, 0) * area_stock_prop(jj, ii)
              
            N_age_sex_area_stock_yr( : , male, jj, ii, 0) = & 
              N_age_sex_area_stock_yr( : , male, all_areas, ii, 0) * area_stock_prop(jj, ii)
        end do 
!        print *, "N_age_sex_area_stock_yr( : , male, 4 , ii, 0)", N_age_sex_area_stock_yr( : , male, 4 , ii, 0)  ! DEBUGGING
    end do  ! end loop over stocks
    
!====== +++ === === +++ === === +++ === ! Try with characteristic equation from Punt (1999)
!    foo = characteristic_eq(lambda_tmp = 1.04d0)   ! Broken function. Do not use unless fixed.
!====== +++ === === +++ === === +++ === Check
!    print *, "Sum stock 1 ages 1+ across areas: "
!    print *, sum(N_age_sex_area_stock_yr(1:age_x, :, 0, 1, 0)) 
!    print *, "N_age_sex_area_stock_yr(1:age_x, :, :, 1, 0)"
!    print *, N_age_sex_area_stock_yr(1:age_x, :, 0, 1, 0)
!====== +++ === === +++ === === +++ === Generate random variates before start looping over years    
!    seq_yrs = (/(yr, yr = 1, yr_max)/) ! Create a sequence of years, e.g. 1 - 100 (cf, seq_yrs = 1:100 in R)  
!
    is_surv_yr = assign_surv_yrs_stock(yr_max, n_stocks, surv_freq) ! Function in PBR_calcs_module

!    print *, 
!    do yr = 0, yr_max
!        print *, yr, (is_surv_yr(yr, jj), jj = 1, n_stocks) 
!    end do
!    print *, "Number of surveys planned for stock 1: ", sum(is_surv_yr(:, 1))
!    print *, "Number of surveys planned for stock 2: ", sum(is_surv_yr(:, 2))

!====== +++ === === +++ === === +++ === Do projections through time with no human caused mortality as reference case    
    yr = 1 ! DEBUGGING
!    b_t = b_eq+(b_max-b_eq)*(1-pow(depl(Year),z))
    b_yr_stock(yr, :) = b_eq + (b_max - b_eq) * (1 - depl_yr_stock(yr - 1, :)**theta) ! Annual birth rate for each stock

    print *, "b_yr_stock(0, :) : ", b_yr_stock(0, :)
    print *, "b_yr_stock(1, :) : ", b_yr_stock(1, :)
!    print *, "depl_yr_stock(yr, :) : ", depl_yr_stock(yr - 1, :)
    print *, "theta: ", theta
    print *, "k_1plus: ", k_1plus
!    Call pop_projection(f_rate = 0.0d0, b_rate = 0.50d0, &
!        N_age_old = N_age_sex_area_stock_yr(:, female, 0, 1, 0),  &
!        N_age_new = N_age_sex_area_stock_yr(:, female, 0, 1, 1))
 
!    print *, "prop_mat_age: "
!    print *, prop_mat_age
!====== +++ === === +++ === === +++ === Do projections through time with no human caused mortality as reference case        
    do yr = 1, yr_max ! Years -- starting at year one, because year zero is in the books and pop has been initialized
        
        do ii = 1, n_stocks ! Stocks
            
            b_yr_stock(yr, ii) = b_eq + (b_max - b_eq) * (1 - depl_yr_stock(yr - 1, ii)**theta) ! Annual birth rate for each stock
!            do jj = 1, n_area ! Areas
                
                do ss = 1, 2  ! Sexes
        
                    Call pop_projection(f_rate = 0.0d0, b_rate = b_yr_stock(yr, ii), &
                      N_age_old = N_age_sex_area_stock_yr(:, ss, 0, ii, yr - 1),  &
                      N_age_new = N_age_sex_area_stock_yr(:, ss, 0, ii, yr))                    
!                    Call pop_projection(f_rate = 0.0d0, b_rate = 0.50d0, &
!                      N_age_old = N_age_sex_area_stock_yr(:, ss, 0, ii, yr - 1),  &
!                      N_age_new = N_age_sex_area_stock_yr(:, ss, 0, ii, yr))                    
                      
                    do aa = 0, age_x ! Ages
                        
                        
                    end do ! End loop over ages
                depl_yr_stock(yr, ii) = sum(N_age_sex_area_stock_yr(1:age_x, :, 0, ii, yr)) / k_1plus(ii)        
                end do ! End loop over sexes
                
!            end do ! End loop over areas
! Determine if this is a survey year and generate abundance estimate if so (surveys assumed to occur at end of year, after birth/death)
! This checks if it's a survey year for each stock, because they are allowed to have different survey intervals            
        ! If(is_surv_yr(yr, n_stocks) == 1) n_best = gen_survey_estimate(true_abundance = 1000.0d0, cv_n = 0.20d0)
        
        end do     ! End loop over stocks 
        
    end do        ! End loop over years
                  ! End loop over number of simulations
!====== +++ === === +++ === === +++ ===    
! DEBUGGING: Look at calculated depletion through time -- ideally want this written to file for plotting in R etc.
!====== +++ === === +++ === === +++ === 
    print *, "Depletion_yr  : "
    do yr = 0, yr_max
        print *, (depl_yr_stock(yr, ii), ii = 1, n_stocks)
    end do 
    
    
!====== +++ === === +++ === === +++ ===    
! Writing results of main population array with spatial-temporal age structure to output file 
!====== +++ === === +++ === === +++ === TODO: Move this into a function in the File IO module       
    open(unit = 1, file = "N_array.out")
    write(1, "(9(a15))") "yr", "stock", "age", "sex", "all_areas", "area_1", "area_2", "area_3", "area_4"  
30  format(4(i15), 5(f15.4)) !       
    do yr = 0, yr_max
        do ii = 1, n_stocks
            do ss = 1, 2
                do aa = 0, age_x ! Note the implicit do loop over areas in next line
                    write(1, 30) yr, ii, aa, ss, (N_age_sex_area_stock_yr( aa , ss, jj , ii, yr), jj = 0, n_area)
                end do  ! End ouput over ages
            end do ! End ouput over sex
        end do ! End output over stock
    end do ! End ouput over years
    !write (1, "(A7,I3)") "hello", 10    ! example of conversion and concantination of character string in fortran    
    close(unit = 1) ! close output file
!====== +++ === === +++ === === +++ ===        
! DEVELOPING
! Test some random number generation
!   open(unit = 1, file = "z_variate.out")
!    write(1, format2) "draw_ID", "z_variate"
!    write(*,*) 'z_variates from random_normal()'
!    do jj = 1, 100
!!            z_variate = random_normal(mean = real(S_adult, kind = 4), sd = real(S_juv, kind = 4)) ! Function located in Random_module.f90
!        z_variate = random_normal(mean = 2., sd = 4.) ! Function located in Random_module.f90
!        write(1, format4) jj, z_variate  ! "format4" declared in file = Format_module.f90
!       if (jj .eq. 1) print *, jj, z_variate
!    end do
! END DEVELOPING    
    
!$         print *, "Compiled with -fopenmp"    ! This is a test for compiling with OpenMP (parallel processor directive <- !$)
         
    print *, "Closing down"
    return     
!====== +++ === === +++ === === +++ ===   
!###### +++ ### ### +++ ### ### +++ ###             
    contains
!====== +++ === === +++ === === +++ ===   
!###### +++ ### ### +++ ### ### +++ ###             
    subroutine initialize_local_vars()
        f_init_ii = 0.d0
        f_yr_stock = 0.d0
        b_init_ii = 0.d0
        b_yr_stock = 0.d0
        depl_yr_stock = 0.d0
        transition_matrix_tmp = 0.d0
        movement_matrix = 0.d0
        N_age_sex_area_stock_yr = 0.d0
        area_stock_prop = 0.d0
        eigv = 0.d0
        seq_yrs = 0
        is_surv_yr = 0
        pbr_yr_stock = 0.d0
        foo_vector = 0.d0 
        NPR_mature = 0.d0 
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
        !
        return
    end subroutine initialize_local_vars

end program main

