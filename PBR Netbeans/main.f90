program main        
!#######################################################    
! File:     main.f
! Author:   John R. Brandon
! Contact:  jbrandon at gmail             
! Date:     Spring 2015
! Purpose:  Run PBR Tier System simulations            
!#######################################################              
! Modules of code to use -- these modules contain subroutines, functions and possibly variable / format declarations.
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
    use Declare_variables_module       ! This module declares variables accessible by the main program here: File = Declare_variables_module.f90
    use PBR_FileIO_Module              ! Routines for reading initial values from files, and for writing output : File = PBR_FileIO_Module.f90
    use initialize_pop                 ! Includes initialization of age structure : File = initialize_Mod.f90
    use calcs                          ! Routines for various calculations (e.g. calculating N_min) : File = PBRmodule.f
    use Format_module                  ! Fortran format statements (everyone's favorite) declared as type character (bit of a hack): File = Format_module.f90
    use random, only : random_normal   ! Module with routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : File = Random_module.f90
    use Generate_random_numbers_module ! Determine if seed for RNG is user defined (for reproducible results) or if seed is based on CPU clock (for different psuedo random variates each time program runs): File = Generate_random_numbers_module.f90
    use PBR_Errorcheck_module          ! [Beta] Contains function 'error_check_input' to do error checking on input values 
    ! use debug
    
    implicit none ! Turns off implicit typing by Fortran; now all variables must be explicitly declared by type
!====== +++ === === +++ === === +++ ===                 ! Local variables 
    real(kind = 8), allocatable :: f_init_ii(:)         ! Initial human caused mortality rates for each stock   
    real(kind = 8), allocatable :: b_init_ii(:)         ! Initial birth rate for each stock 
    real(kind = 8) :: sum_1plus_tmp, sum_mature_tmp     ! DEBUGGING
    real(kind = 8) :: foo, foo1                         ! DEBUGGING
    real(kind = 8) :: start, finish                     ! For timing / optimizing code
    real(kind = 8) :: objf_brent
    real(kind = 8) :: brent    ! Function brent() :: file = Brent.f90 (TODO : Make this a Roots_and_Mins_module.f90
    integer(kind = 4) :: io_error     
!====== +++ === === +++ === === +++ === 
    call read_inits()              ! Read initial values from 'input.par' file. PBR_FileIO_Module contains subroutine 'read_inits()' 
    io_error = 0                   ! Check for input value errors (values out of bounds, etc) 
    io_error = error_check_input() ! Function that contains error checking code. Located in file: PBR_Errorcheck_module.f90
    if(io_error .ne. 0) then       ! Check error code
        print *, "Error Code: ", io_error  ! If error, print code
        stop                       
    end if
!====== +++ === === +++ === === +++ === ! Initialize
    age_x = a_m + 1                     ! Might eventually want to change this to take user input for age_x (?)                   
!====== +++ === === +++ === === +++ === ! Given input, allocate array dimensions at run-time 
    allocate(f_init_ii(1:n_stocks))     ! Stock specific initial human caused mortality rates (can have stock specific init_depl's) 
    allocate(b_init_ii(1:n_stocks))     ! Stock specific initial birth rate
    
    allocate(Female_age(0:age_x))  ! Females at age vector
    allocate(Male_age(0:age_x))    ! Males at age vector
    allocate(S_age(0:age_x))       ! Survival at age vector
    allocate(prop_mat_age(0:age_x))! Proportion mature at age
    allocate(NPR_age(0:age_x))     ! Numbers at age per recruit vector
    allocate(Nage(0:age_x))        ! Numbers at age per recruit vector    
    allocate(NPR_age_tmp(0:age_x)) ! Numbers at age per recruit vector, used to solve for initial human caused mortality rate         
    allocate(Nage_imm_0(0:age_x))  ! Immature numbers at age vector
    allocate(Nage_mat_0(0:age_x))  ! Mature numbers at age vector
    allocate(prop_NPR(0:age_x))    ! Rescaled numbers at age per recruit vector (sums to 1.0 over ages)
    allocate(selectivity(0:age_x)) ! Selectivity at age -- currently assuming knife-edge selectivity at age a_r (same selectivity for each stock)
                                   ! Note that all of the above start with index 0, to be consistent with age 0 notation
    allocate(Nplus(0:yr_max))      ! Vector of age 1+ population size over projection years 
    allocate(Ntot(0:yr_max))       ! Total (0+) population size each year of projection
    allocate(N_calf(0:yr_max))     ! Vector of calf production for each projection year

!   allocate(depl_ii_tt(0:yr_max, 1:n_stocks))     ! Depletion in year t (row) of stock i (column) 
!   allocate(NPR_age_ii(0:age_x, 1:n_stocks))      ! Numbers at age per female recruit (row) for stock i (column)
!   allocate(NPR_oneplus_ii(1:n_stocks))           ! Total number of age 1+ per female recruit for stock i 

    Call set_random_seed() ! Uses input from read_inits() to set seed for RNG -- see comment after 'use Generate_random_numbers_module' statement above
    
    Call assign_par_vectors(a_r, a_m, age_x, S_adult, & ! Initialize vectors: survival, selectivity and proportion mature at age
                            S_juv, S_age, selectivity, prop_mat_age) ! Currently these three vectors assumed to be identical for each stock                        

    Call calc_NPR_age(f_rate = 0.0d0, &                 ! Calculate Numbers per female recruit, with no human caused mortality
        N_recruits = b_sex_ratio, N_age = NPR_age, &    ! d0 suffix for double precision to match argument type                                             
        sum_1plus = NPR_oneplus, sum_mature = sum_mature_tmp) 
        
    b_eq = 1 / sum_mature_tmp                           ! Calculate equilibrium birth rate (at K); same for both stocks
        
    do ii = 1, n_stocks                         ! Calculate the initial age structure for each stock; can differ, e.g. initial depletions not equal 
        init_depl_i = init_depl(ii)             ! New value for global variable init_depl_i, that is used by function Initial_F()
        b_init = b_eq + (b_max - b_eq) * (1 - (init_depl_i ** theta))  ! Initial birth rate for stock i
        b_init_ii(ii) = b_init                                         ! Can be different for each stock if they have different init_depl_ii  
        objf_brent = BRENT(ax = 0.0d0, bx = 0.10d0, cx = 1.0d0, func = initial_F, &
            tol = 0.0000001d0, xmin = f_init_ii(ii)) ! Calculate the initial human caused mortality rate that results in stable age-structure at initial depletion
        print *, "Stock ", ii            
        print *, "f_init_ii: ", f_init_ii(ii)         
    end do
    
    print *, "Initial depletion by stock"
    print *, init_depl    
    print *, "Initial birth rate by stock"    
    print *, b_init_ii    
         
!$         print *, "Compiled with -fopenmp"    ! This is a test for compiling with OpenMP (parallel processor directive <- !$)
         
    print *, "Closing down"
    return      
end program main

