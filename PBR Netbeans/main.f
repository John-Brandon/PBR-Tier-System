         program main        
!#######################################################    
! File:   main.f
! Author: John R. Brandon
! Contact: jbrandon at gmail             
!
! Created on April 23, 2015, 6:38 PM
!#######################################################
! Calculate the stable age structure and equilibrium birth rate based on numbers-per-recruit calculations        
!#######################################################
             
! modules of code to use         
         use calcs 
         use debug    
         use Mod_initialize ! includes initialization of age structure
         
         implicit none ! turns off implicit typing; all variables must be declared by type
         
	 integer(kind = 4) :: ii, jj, kk, ll ! counters
	 integer(kind = 4) :: xx, nn
         integer(kind = 4) :: a_m, npr
         
	 real(kind = 8) :: aa, bb
	 real(kind = 8) :: start, finish
	 real(kind = 8) :: N_best, CV_N, z_variate
	 real(kind = 8) :: N_min
	 real(kind = 8) :: x_min, x_max, step_CV
         real(kind = 8) :: yy
         real(kind = 8) :: S_juv, S_adult, delt_s ! delt_s is the difference between adult and juvenile survival rates
         
         !real(kind = 8), external :: debugz ! debugz 
	 character*72 ccon
	 integer ic
	 
10	 format(a10 f5.3 a10 f8.3)
20	 format(a10 a10) ! header for CV_N output file
30	 format(f10.3 f10.3) ! header for CV_N output file

	 CALL RANDOM_SEED() ! initialize seed for random number generator
		 
!		write(*,*) 'Calling cpu_time'
!		CALL test_cpu_time(start, finish, aa, bb)
         a_m = 3 ! of course will want to move these hard-coded values to be read from input file
         S_adult = 0.95
         S_juv = 0.80
         delt_s = S_adult - S_juv
         npr = -99
         write(*,*) "a_m = : ", a_m
         yy = debugz(a_m)
         write(*,*) "yy = : ", yy

         print *, "init_age_distribution"
         Call init_age_distribution(a_m, npr, S_adult, delt_s)

         
! Test some random number generation
	 write(*,*) "Random standard uniform deviate"
!	 z_variate = RAN3(-1)

	 z_variate = RAND(1) ! RAND is an intrinsic GNU fortran compiler function
	 write(*,*) 'z_variate'
	 write(*,*) z_variate
	 write(*,*) RAND(0), RAND(0), RAND(0), RAND(0), RAND(0)
	 write(*,*) RAND(1), RAND(1), RAND(1), RAND(1)

	 do 12 jj = 1, 10
		 CALL RANDOM_NUMBER(z_variate) ! a better intrinsic uniform (0,1) Ran Num Generator
		 write(*,*) 'z_variate RANDOM_NUMBER(r)'
		 write(*,*) z_variate
12	 continue

	 	 
! Set up code for file I/O, e.g. reading N_best, CV_N, z from input file

	 N_best = -99.99
	 CV_N = -99.99
	 z_variate = -99.99

	 write(*,*) 'z_variate'
	 write(*,*) z_variate

! open file
	 open(unit = 1, file = "par_input_file.txt")
! read file
!  here: N_best = 1042., CV_N = 0.65, z_variate = 0.842
	 read(1,*)
	 read(1,*) ccon,ic
	 write(*,*) 'ccon: ', ccon
	 write(*,*) 'ic: ', ic
	 read(1,*)
	 read(1, '(/f8.3 /// f8.3 /// f8.3)') N_best, CV_N, z_variate ! read from input file
	 close(unit = 1) ! close input file
	 
	 write(*,*) 'z_variate after read'
	 write(*,*) z_variate

	 write(*,*) 'N_best'
	 write(*,*) N_best
	 write(*,*) 'CV_N'
	 write(*,*) CV_N

	 N_min = CALC_N_min(N_best, CV_N, z_variate)
	 write(*,*) 'N_min'
	 write(*,*) N_min
	 	 
! Loop over a range of values for CV_N and output N_min, given fixed N_best
	 x_min = 0.1
	 x_max = 1.5
	 step_CV = 0.1

	 nn = nint((x_max - x_min)/step_CV)
	 write(*,*) 'nn'
	 write(*,*) nn

	 yy = CALC_N_min(N_best, x_min, z_variate)
	 write(*,*) 'yy'
	 write(*,*) yy
	 	 

! Test writing results to an output file, which will be read by R code
	 open(unit = 1, file = "Nmin.out")
	 write(1, 20) "CV_N", "N_min"
	 
	 Do 11 ii = 1, nn+1
	  yy = CALC_N_min(N_best, step_CV*ii, z_variate)
	  write(*,10) "CV_N: ", step_CV*ii, " : N_min = ", yy !
	 write(1, 30) step_CV*ii, yy 
11	 continue 

	 close(unit = 1) ! close output file
         
         write(*,*) "Closing down"
         return
         end program main
	   

! *************************************************************************      
! Function to open and read input file with parameters for this run
!  Testing at this stage by reading (i) N_best (ii) CV and (iii) standard normal variate, z
! *************************************************************************            
	subroutine read_par_file(par_file_name, N_best, CV_N, z_variate)
	
	 real N_best, CV_N, z_variate
	 integer in
     
	 in = 1
     
! put code to open input par file 
	 open(unit = 1, file = "par_input_file.txt")
! put code to read input par file 
! will need a format statement
	 read(1, '(/f8.3 // f8.3 // f8.3)') N_best, CV_N, z_variate
	 close(UNIT = 1)

	 write(*,*) 'N_best'
	 write(*,*) N_best

	 
	return
	end subroutine read_par_file 
	
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
	
