         program main        
!#######################################################    
! File:   main.f
! Author: John R. Brandon
! Contact: jbrandon at gmail             
!
! Created on April 23, 2015, 6:38 PM
!#######################################################              
! Modules of code to use -- these modules contain subroutines, functions and possibly variable / format declarations.
!  Each module of code is contained in a separate file (e.g. Declare_variables_module.f90).
!  The code in each module file needs to be compiled and linked with the compiled main program to produce an executable.
!  Note: The order in which these modules are compiled relative to the main program may matter. (using the gfortran compiler, at least) 
!   Assume that the modules need to be compiled before they are linked with the main program to form an executable. 
!
!   This can be done in one line from the shell prompt ($); but again, note the file order in the example below:
!    $ gfortran module1.f90 module2.f90 main.f -o desired_executable_name_here
!
!   Whereas, the shell command below might not link properly, because the main program is compiled before module2.f90:
!    $ gfortran module1.f90 main.f module2.f90 -o desired_executable_name_here
!             
!   This can be a source of maddening errors during compiling if you're developing code in an integrated development environment (IDE). 
!   If you can't get things to compile and link in an IDE (e.g. 'Error can't find .mod file') you can try and figure out the settings the IDE
!    uses when ordering its list of files to the Fortran compiler, and then change those settings to make sure the module files are compiled in proper order before linking. 
!    Alternatively, you can use the shell prompt to manually control for the correct order of file compilation (as in the first $ example above).              
!#######################################################
!   General comments:
!   (i) At present, the population dynamics model assumes this order for entering female adult-hood: 
!       (1st) Maturity (ovulation) on her a_m'th birthday -> (2nd) Adult mortality rate applied during her a_m'th year -> (3rd) First partuition (one year after her a_m'th birthday, if she survives that first year of adult-hood).
!           (a) a_m is the age at which females mature (males are assumed to be a non-limiting factor for reproduction / birth rates)
!            (b) The code for the population dynamics model would need to be revised to take into account reproductive senescence (as might be expected for at least some "black-fish", e.g. killer whales).
!             (c) Minimum gestation length after reach a_m is assumed to be one year. TODO: make sure birth rates bounded between 0.0 and 1.0 
!   (ii)             
!#######################################################
             
         use Declare_variables_module       ! This module declares the global variables (like a 'Common' block in FORTRAN 77): File = Declare_variables_module.f90
         use PBR_FileIO_Module              ! Routines for reading initial values from files, and for writing output : File = PBR_FileIO_Module.f90
         use initialize_pop                 ! Includes initialization of age structure : File = initialize_Mod.f90
         use calcs                          ! Routines for various calculations (e.g. calculating N_min) : File = PBRmodule.f
         use Format_module                  ! Format statements declared as type character (bit of a hack): File = Format_module.f90
         use random, only : random_normal   ! Module with routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : File = Random_module.f90
         use Generate_random_numbers_module ! Determine if seed for RNG is user defined (for reproducible results) or if seed is based on CPU clock (for different psuedo random variates each time program runs): File = Generate_random_numbers_module.f90
!         use debug
         
         implicit none ! Turns off implicit typing by Fortran; now all variables must be explicitly declared by type

         real(kind = 8) :: foo ! DEBUGGING
         
         call read_inits() ! Read initial values from 'input.par' file. PBR_FileIO_Module contains subroutine 'read_inits()' 

         age_x = a_m + 1    ! At present, population dynamics model assumes this order for entering adult-hood: maturity (ovulation) -> adult mortality rate -> partuition
         allocate(Female_age(0:age_x))  ! allocate array dimensions based on input file 
         allocate(Male_age(0:age_x))    ! Males at age vector
         allocate(S_age(0:age_x))       ! Survival at age vector
         allocate(prop_mat_a(0:age_x))  ! Proportion mature at age
         allocate(NPR_age(0:age_x))     ! Numbers at age per recruit vector
         allocate(Nage(0:age_x))        ! Numbers at age per recruit vector, used to solve for initial human caused mortality rate         
         allocate(Nage_imm_0(0:age_x))  ! Immature numbers at age vector
         allocate(Nage_mat_0(0:age_x))  ! Mature numbers at age vector
         allocate(prop_NPR(0:age_x))    ! Rescaled numbers at age per recruit vector (sums to 1.0 over ages)
         allocate(selectivity(0:age_x)) ! Selectivity at age -- values initialized in by Initial_F() routine in 'initialize_pop' module: File = initialize_Mod.f90
                                        ! Note that all of the above start with index 0, to be consistent with age 0 notation
         
         call set_random_seed() ! Uses input from read_inits() to set seed for RNG -- see comment after 'use Generate_random_numbers_module' statement above
         
         delt_s = S_adult - S_juv ! Calculate difference between juvenile and adult survival
         print *, "delt_s: ", delt_s
         print *, "S_adult: ", S_adult
         print *, "S_juv: ", S_juv         
! For DEBUGGING (uncomment "use debug", in linked modules above if want to test / debug a routine here that is being developed in module 'debug')        
!         write(*,*) "a_m = : ", a_m
!         yy = debugz(a_m)
!         write(*,*) "yy = : ", yy

! DEVELOPING
         print *, "init_age_distribution" ! This is translated from Do_NPR() function from JRB's ADMB gray whale code
         Call initialize_age_struc() !(a_m, npr, S_adult, delt_s) -- initializes age structured vector based on numbers per recruit approach

         print *, "Calling newtons_root()"         
         foo = newtons_root()           ! 
         print *, "initial_F: ", foo ! DEBUGGING

!         foo = initial_F() ! Calculate the initial human caused mortality rate 
!         print *, "initial_F: ", foo ! DEBUGGING
         
         print *, "rescale_NPR()"
         Call rescale_NPR() ! Rescale the numbers per recruit to initial population size, i.e. scale up to the initial numbers at age vector
! DEVELOPING
! Test some random number generation
!	 open(unit = 1, file = "z_variate.out")
!	 write(1, format2) "draw_ID", "z_variate"
!         print *, real(S_adult, kind = 4)
!         write(*,*) 'z_variates from random_normal()'
!         do jj = 1, 100
!!            z_variate = random_normal(mean = real(S_adult, kind = 4), sd = real(S_juv, kind = 4)) ! Function located in Random_module.f90
!             z_variate = random_normal(mean = 2., sd = 4.) ! Function located in Random_module.f90
!             write(1, format4) jj, z_variate  ! "format4" declared in file = Format_module.f90
!            if (jj .eq. 1) print *, jj, z_variate
!         end do
! END DEVELOPING
         
!$         print *, "Compiled with -fopenmp"
         write(*,*) "Closing down"
         return
         end program main