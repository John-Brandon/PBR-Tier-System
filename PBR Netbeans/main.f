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
         use debug    ! testing

         implicit none ! turns off implicit typing by Fortran; now all variables must be explicitly declared by type
               
10	 format(a10 f5.3 a10 f8.3)
20	 format(a10 a10) ! header for CV_N output file
30	 format(f10.3 f10.3) ! header for CV_N output file

	 call random_seed() ! initialize seed for random number generator
         
         call read_inits()! (CV_N, CV_MORTALITY, THETA, R_MAX, F_R, INIT_DEPL, &
         !LOWER_TAIL, YR_MAX, SURV_FREQ, KK, IPAR)
         
         delt_s = S_adult - S_juv ! calculate difference between juvenile and adult survival
         ! maxage = 59 ! 
       
         write(*,*) "a_m = : ", a_m
         yy = debugz(a_m)
         write(*,*) "yy = : ", yy

         print *, "init_age_distribution"
         Call initialize_age_struc(a_m, npr, S_adult, delt_s)
         
! Test some random number generation
!	 do jj = 1, 10
!            call random_number(z_variate) ! intrinsic uniform (0,1) Random Number Generator
!            write(*,*) 'z_variate RANDOM_NUMBER(r)'
!            write(*,*) z_variate
!         end do

         
         write(*,*) "Closing down"
         return
         end program main
	   


	
! *************************************************************************
! *************************************************************************      
	subroutine test_cpu_time(start, finish, aa, bb)
	
	 real start, finish
	 real aa, bb
	 call cpu_time(start)
! put code to test here
! Want to benchmark matrix multiplication in FORTRAN vs. that in R (for age-structured model)
	  
	 call cpu_time(finish)
	 print '("Time = ",f6.3," seconds.")',finish-start
	return
	end subroutine test_cpu_time
	
