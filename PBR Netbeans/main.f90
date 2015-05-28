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
    use Format_module                  ! Fortran format statements (everyone's favorite) declared as type character (bit of a hack): Format_module.f90
    use random, only : random_normal   ! Routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : Random_module.f90
    use Generate_random_numbers_module ! Determine if seed for RNG is user defined (for reproducible results) or if seed is based on CPU clock (for different psuedo random variates each time program runs): Generate_random_numbers_module.f90
    use PBR_Errorcheck_module          ! Contains function 'error_check_input' to do error checking on input values [Very Beta]
    ! use debug
!====== +++ === === +++ === === +++ === ! Turns off implicit typing by Fortran; all variables must be explicitly declared by type
    implicit none 
!====== +++ === === +++ === === +++ ===                 ! Constant parameters
    integer(kind = 4), parameter :: stock_1 = 1         ! Indexes for stock structure array
    integer(kind = 4), parameter :: stock_2 = 2
    integer(kind = 4), parameter :: all_areas = 0
    integer(kind = 4), parameter :: area_1 = 1
    integer(kind = 4), parameter :: area_2 = 2
    integer(kind = 4), parameter :: area_3 = 3    
    integer(kind = 4), parameter :: area_4 = 4 
    integer(kind = 4), parameter :: n_area = 4       ! For counters
    integer(kind = 4), parameter :: female = 1          
    integer(kind = 4), parameter :: male = 2   
!====== +++ === === +++ === === +++ ===                 ! Local variables 
    real(kind = 8), allocatable :: f_init_ii(:)         ! Initial human caused mortality rates for each stock   
    real(kind = 8), allocatable :: b_init_ii(:)         ! Initial birth rate for each stock 
! Array of numbers-at-sex and age for each stock in each sub-area by year of projection    
    real(kind = 8), allocatable :: N_yr_ii_jj_mf_age(:,:,:,:,:)  ! (yr, stock, sub-area, sex, age) 
    real(kind = 8), allocatable :: N_age_mf_jj_ii_yr(:,:,:,:,:)  ! (age, sex, sub-area, stock, yr
    real(kind = 8), allocatable :: area_stock_prop(:,:)    ! Percentage of each stock in each area
    real(kind = 8), allocatable :: transition_matrix(:, :) ! Transition matrix (with survival and birth rates)
    real(kind = 8), allocatable :: eigv(:)                 ! Eigen vector of transition matrix, returned by power_method()    
    real(kind = 8) :: NPR_mature                        ! Numbers mature per female recruit. Used to calculate b_eq (birth rate at K)
!    real(kind = 8) :: sum_1plus_tmp                     ! DEBUGGING
    real(kind = 8) :: foo, foo1                         ! DEBUGGING
    real(kind = 8) :: start, finish                     ! For timing / optimizing code
    real(kind = 8) :: objf_f_init                       ! Objective function for finding f_init resulting in stable age structure
    real(kind = 8) :: brent    ! Function brent() :: file = Brent.f90 (TODO? : Create and add this to a Roots_and_Mins_module.f90)
    real(kind = 8) :: lambda                            ! Dominant real eigen value of the transition matrix
    integer(kind = 4) :: io_error                       ! Error flag for checking initial values in input.par
    integer(kind = 4) :: ii, jj, kk, aa, yr, ss         ! Counters for indexing loops  
    integer(kind = 4) :: it_num                         ! Iteration number returned from power_method() for calculating lambda
    character(len = 10) :: stock_name                   ! For printing output tables
!====== +++ === === +++ === === +++ === ! Read initial values and do error checking
    call read_inits()              ! Read initial values from 'input.par' file. PBR_FileIO_Module contains subroutine 'read_inits()' 
    io_error = 0                   ! Check for input value errors (values out of bounds, etc) 
    io_error = error_check_input() ! Function that contains error checking code. Located in file: PBR_Errorcheck_module.f90
    if(io_error .ne. 0) then       ! Check error code
        print *, "Error Code: ", io_error  ! If error, print code
        stop                       ! If error, exit program
    end if                         ! TODO? Move this error checking to be called at the end of read_inits() procedure? 
!====== +++ === === +++ === === +++ === ! Initialize
    age_x = a_m + 1                     ! Might eventually want to change this to take user input for age_x (?)                   
!====== +++ === === +++ === === +++ === ! Given input, allocate array dimensions at run-time (local variables)
    allocate(f_init_ii(1:n_stocks))     ! Stock specific initial human caused mortality rates (can have stock specific init_depl's) 
    allocate(b_init_ii(1:n_stocks))     ! Stock specific initial birth rate, given stock specific initial depletion levels
!    allocate(N_yr_ii_jj_mf_age(0:yr_max, 1:n_stocks, 0:4, 1:2, 0:age_x)) ! Array of numbers-at-age and sex for each stock by year of projection
    allocate(N_age_mf_jj_ii_yr(0:age_x, 1:2, 0:4, 1:n_stocks, 0:yr_max))
    allocate(area_stock_prop(1:4, 1:n_stocks))  ! Hard-coded four areas here, not great
    allocate(transition_matrix(0:age_x, 0:age_x))! Transition matrix, added here while developing methods for population projections and solving for juvenile survival (given lamda_max)
    allocate(eigv(0:age_x))         ! Eigen vector of transition matrix, returned by power_method()    
!====== +++ === === +++ === === +++ === ! Given input, allocate array dimensions at run-time (global variables)    
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
!   allocate(depl_ii_tt(0:yr_max, 1:n_stocks))     ! Depletion in year t (row) of stock i (column) 
!   allocate(NPR_age_ii(0:age_x, 1:n_stocks))      ! Numbers at age per female recruit (row) for stock i (column)
!   allocate(NPR_oneplus_ii(1:n_stocks))           ! Total number of age 1+ per female recruit for stock i 
!====== +++ === === +++ === === +++ === ! Set seed for RNG -- see comment after 'use Generate_random_numbers_module' statement above
    Call set_random_seed() ! Uses input from read_inits() to set seed based on user input (reproducible), or based on CPU clock
!====== +++ === === +++ === === +++ === ! Initialize vectors: survival, selectivity and proportion mature at age   
    Call assign_par_vectors(a_r, a_m, age_x, S_adult, S_juv, &  
                            S_age, selectivity, prop_mat_age) ! Currently these three vectors are identical for each stock                        
!====== +++ === === +++ === === +++ === ! Solve for juvenile survival rate that results in specified r_max
    transition_matrix = assign_transition_matrix(a_m, age_x, b_max, S_age, prop_mat_age) ! Assign values to transition matrix
    print *, "transition_matrix: "
    do aa = 0, age_x
        write (*,"(100f4.3)") (transition_matrix(aa, jj), jj = 0, age_x)    ! The 100f... is a bit of a hack. Works if <= 100 columns to be printed
    end do
!====== +++ === === +++ === === +++ === ! Test power method for calculating dominant eigen value (and eigen vector -- is this right evector, i.e. stable age distribution??)
    print *, 
    print *, "Calling power method to calculate Lambda and eigv: "
    eigv = 1.d0
    eigv = eigv / sum(eigv) ! Debugging: See if this helps with occasional "sticky" cycle in trying to solve for lambda below
    print *, "eigv: "
    print *, eigv
    print *, "size(transition_matrix): "
    print *, size(transition_matrix)
    Call power_method(n =(age_x + 1), a = transition_matrix, y = eigv, & ! n, a, y, it_max, tol, lambda, it_num 
            it_max = 200, tol = 0.00000001d0, lambda = lambda, it_num = it_num) 
    print *, "Power method iterations: ", it_num
    print *, "Eigen Vector: "
    print *, eigv
    print *, "Lambda: ", lambda
    print *,
    if (it_num > 200) then
        print *, "Error with power method finding Lambda"
        stop
    end if
!====== +++ === === +++ === === +++ === ! Assign percentage of each stock to areas 
    area_stock_prop = assign_area_stock_prop(p_a1_s1, p_a2_s1, p_a2_s2, p_a3_s2, p_a4_s2) 
!    print *, "area_stock_prop: "
!    do jj=1, n_area
!        print *, (area_stock_prop(jj,ii), ii = 1, n_stocks)
!    end do
!====== +++ === === +++ === === +++ === ! Calculate Numbers per female recruit, with no human caused mortality (f_init = 0.0)
    Call calc_NPR_age(f_rate = 0.0d0, &                         ! d0 suffix for double precision to match argument type                                             
        N_recruits = b_sex_ratio, N_age = NPR_age, &            ! Calc NPR_age, NPR_oneplus and NPR_mature (F = 0)
        sum_1plus = NPR_oneplus, sum_mature = NPR_mature) 
!====== +++ === === +++ === === +++ === ! Calculate equilibrium birth rate (at carrying capacity)
    b_eq = 1 / NPR_mature               ! Equilibrium birth rate. Equal for both stocks under assumption of identical life histories                
!====== +++ === === +++ === === +++ === ! Calculate the initial age structure for each stock         
    do ii = 1, n_stocks                 ! Initial age structures for each stock can differ, e.g. initial depletion may not be equal
        init_depl_i = init_depl(ii)     ! New value for global variable init_depl_i. Used by Initial_F(), as called from brent()
        b_init = b_eq + (b_max - b_eq) * (1 - (init_depl_i ** theta))  ! Initial birth rate for stock i
        b_init_ii(ii) = b_init                                         ! Storing initial birth rates in vector by stock 
! Calculate initial human caused mortality rate (f_init) that results in stable age-structure at initial depletion        
        print *, 
        print *, "Stock: ", ii                        ! DEBUGGING  
        print *, "b_init: ", b_init
        print *, "prop_NPR: "
        prop_NPR = NPR_age / sum(NPR_age)
        print *, prop_NPR
        objf_f_init = BRENT(ax = 0.0d0, bx = 0.10d0, cx = 1.0d0, func = initial_F, & ! See Brent.f90 for details on arguments
                            tol = 0.0000001d0, xmin = f_init_ii(ii)) 
        print *, "f_init_ii: ", f_init_ii(ii)        ! DEBUGGING
        Call rescale_NPR(k_1plus_tmp = k_1plus(ii), initial_oneplus_tmp = NPR_oneplus, & ! Could turn this into a function
                        N_age_unscaled = NPR_age_tmp, N_age_scaled = N_age) ! N_age returned as scaled numbers of females at age 
!        N_yr_ii_jj_mf_age(0, ii, 0, female, :) = N_age ! Numbers(yr, stock, sub-area, sex, age) 
!        N_yr_ii_jj_mf_age(0, ii, 0, male,  :) = N_age   ! Assuming a 50:50 sex ratio at birth
        N_age_mf_jj_ii_yr(:, female, 0, ii, 0) = N_age  ! This approach to accessing array (with ":" on left) should be faster
        N_age_mf_jj_ii_yr(:, male, 0, ii, 0) = N_age  ! This approach to accessing array (with ":" on left) should be faster        
        print *, "init_depl_i: ", init_depl_i        ! DEBUGGING        
        print *, "N_age female 1+ after scaling"            ! DEBUGGING
        print *, sum(N_age(1 : age_x))               ! DEBUGGING
        do jj = 1, n_area                            ! Allocate stock abundance across areas
            N_age_mf_jj_ii_yr( : , female, jj, ii, 0) = N_age_mf_jj_ii_yr( : , female, all_areas, ii, 0) * area_stock_prop(jj, ii)
            N_age_mf_jj_ii_yr( : , male, jj, ii, 0) = N_age_mf_jj_ii_yr( : , male, all_areas, ii, 0) * area_stock_prop(jj, ii)
!            N_yr_ii_jj_mf_age(0, ii, jj, female, :) = N_age ! Numbers(yr, stock, sub-area, sex, age) 
!            N_yr_ii_jj_mf_age(0, ii, jj, male,  :) = N_age   ! Assuming a 50:50 sex ratio at birth
        end do 
    end do  
    
    print *, "prop_mat_age: " 
    print *, prop_mat_age
    print *, "S_age: "
    print *, S_age
!====== +++ === === +++ === === +++ === ! Start projections over years
    do yr = 1, yr_max
        
    end do
    
!====== +++ === === +++ === === +++ ===    
! Test writing results to a formatted output file 
! Eventually move this into the File IO module    
!====== +++ === === +++ === === +++ ===    
    open(unit = 1, file = "Ntest.out")
    write(1, "(7(a10))") "stock", "age", "all_areas", "area_1", "area_2", "area_3", "area_4"
30  format(i10, i10, 5(f10.4))  
    do ii = 1, n_stocks
        do kk = 0, age_x
            write(1, 30) ii, kk, (N_age_mf_jj_ii_yr( kk , female, jj , ii, 0), jj = 0, n_area)
        end do  
    end do
    write (1, "(A7,I3)") "hello", 10    ! example of conversion and concantination of character string in fortran    
    close(unit = 1) ! close output file
!====== +++ === === +++ === === +++ ===
    
    print *, "Sum stock 1 ages 1+ across areas: "
    print *, sum(N_age_mf_jj_ii_yr(1:age_x, 1:n_area, 1, 1, 0)) 
! Take a look with some output
    
!    print *, "N_age:"
!    print *, N_age
!    print *, "area_stock_prop(1,1): "
!    print *, area_stock_prop(1,1)
!    print *, "N_age_mf_jj_ii_yr(:, female, 1, 1, 0): "
!    print *, N_age_mf_jj_ii_yr(:, female, 1, 1, 0)
    !! DO PROJECTION from year 0 to year 1
    
!    print *, "N_yr_ii_jj_mf_age(0, 1, male, 0, :)"
!    print *, N_yr_ii_jj_mf_age(all_areas, stock_1, 0, male, :)
!    print *, "N_yr_ii_jj_mf_age(0, 2, male, 0, :)"
!    print *, N_yr_ii_jj_mf_age(all_areas, stock_2, 0, male, :)    
    
!    Call rescale_NPR(k_1plus_tmp = k_1plus(1), initial_oneplus_tmp = NPR_oneplus, &
!                     N_age_unscaled = NPR_age_tmp, N_age_scaled = N_age)
!    print *, "NPR_oneplus"           ! DEBUGGING
!    print *, NPR_oneplus                               ! DEBUGGING        
!    print *, "init_depl(1)"           ! DEBUGGING
!    print *, init_depl(1)                               ! DEBUGGING    
!    print *, "N_age, after scaling"           ! DEBUGGING
!    print *, N_age                               ! DEBUGGING
    
!    print *, "Initial depletion by stock"            ! DEBUGGING
!    print *, init_depl                               ! DEBUGGING
!    print *, "Initial birth rate by stock"           ! DEBUGGING
!    print *, b_init_ii                               ! DEBUGGING
!    print *, "Carrying capacity for each stock"           ! DEBUGGING
!    print *,  k_1plus                              ! DEBUGGING
!    print *, "NPR_age from main: " ! DEBUGGING
!    print *, NPR_age ! DEBUGGING
!    print *, "NPR_oneplus from main: " ! DEBUGGING
!    print *, NPR_oneplus ! DEBUGGING
!    print *, "sum_mature_tmp from main: " ! DEBUGGING
!    print *, sum_mature_tmp ! DEBUGGING 
       
!   print *, "Calling rescale_NPR()"
!   Call rescale_NPR(k_1plus(1), init_depl(1), initial_oneplus(1)) ! Rescale the numbers per recruit to initial population size, i.e. scale up to the initial numbers at age vector

!   print *, "Calling BRENT" ! DEBUGGING 
!!   ! ax,bx,cx,func,tol,xmin
!   call cpu_time(start)
!   objf_brent = BRENT(ax = 0.0d0, bx = 0.10d0, cx = 1.0d0, func = initial_F, &
!        tol = 0.0000001D0, xmin = f_init_tmp(1)) ! Calculate the initial human caused mortality rate 
!   call cpu_time(finish)
!   print '("Time with calc_NPR_age = ",f8.7," seconds.")',finish-start   
!    print *, "f_init_tmp(1): "
!    print *, f_init_tmp(1)
!
!   call cpu_time(start)
!   objf_brent = BRENT(ax = 0.0d0, bx = 0.10d0, cx = 1.0d0, func = initial_F_oldcode, &
!        tol = 0.0000001D0, xmin = f_init_tmp(1)) ! Calculate the initial human caused mortality rate 
!   call cpu_time(finish)
!   print '("Time with initial_F_oldcode = ",f8.7," seconds.")',finish-start 
!    print *, "f_init_tmp(1): "
!    print *, f_init_tmp(1)   
!   
!   print *, "f_init_tmp(1): ", f_init_tmp(1)
!   print *, "objf_brent: ", objf_brent    
        
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
end program main

