         program main        
!#######################################################    
! File:   main.f
! Author: John R. Brandon
! Contact: jbrandon at gmail             
!
! Created on April 23, 2015, 6:38 PM
!#######################################################    
             
! Modules of code to use -- these modules contain subroutines, functions and possibly variable declarations.
!  The modules are contained in seperate files (e.g. PBR_FileIO_Module.f90).
!  The module files need to be compiled and linked with the compiled main program to produce an executable.         
         use Declare_variables_module
         use PBR_FileIO_Module ! routines for reading initial values from files and writing output : PBR_FileIO_Module.f90
         use initialize_pop ! includes initialization of age structure : initialize_Mod.f90
         use calcs    ! routines for various calculations (e.g. calculating N_min) : PBRmodule.f
         use Format_module ! module with various format statements
         use random, only : random_normal ! Module with routines for random number generators -- only using random_normal at this stage

         implicit none ! turns off implicit typing by Fortran; now all variables must be explicitly declared by type
         
         integer, allocatable :: seed(:)
         integer :: n, clock
         
         call cpu_time(start)
! put code to test here
! Want to benchmark matrix multiplication in FORTRAN vs. that in R (for age-structured model)
	  
	 call cpu_time(finish)
	 print '("Time = ",f6.3," seconds.")',finish-start

         call random_seed(size = n)
         allocate(seed(n))
         print *, "seed(n): ", seed ! in gfortran 4.2 compiler on Mac OS X 10.9, seed is vector of length eight (zeros at this stage)
         print *, "n: ", n
         call system_clock(count = clock)
         print *, "clock: ", clock
         do ii = 1, n
             seed(ii) = clock + 37*(ii-1) ! recommended function for setting seeds of RNG because 37 is prime number and RNG likes primes
             print *, "ii / seed(ii): ", ii, seed(ii)
         enddo
         call random_seed(put = seed) ! initialize seed for random number generator
         
         call read_inits()! (CV_N, CV_MORTALITY, THETA, R_MAX, F_R, INIT_DEPL, &
         !LOWER_TAIL, YR_MAX, SURV_FREQ, KK, IPAR)
         
         delt_s = S_adult - S_juv ! calculate difference between juvenile and adult survival
         ! maxage = 59 ! 
       
!         write(*,*) "a_m = : ", a_m
!         yy = debugz(a_m)
!         write(*,*) "yy = : ", yy

         print *, "init_age_distribution"
         Call initialize_age_struc(a_m, npr, S_adult, delt_s)
         
! Test some random number generation
	 open(unit = 1, file = "z_variate.out")
	 write(1, format2) "draw_ID", "z_variate"
         
         write(*,*) 'z_variates from random_normal()'
         do jj = 1, 1000
!            z_variate = r8_normal_01(12345) ! 
            z_variate = random_normal() ! Function located in 
            !write(*,*) z_variate
            write(1, format4) jj, z_variate 
         end do

         
         write(*,*) "Closing down"
         return
         end program main
	
