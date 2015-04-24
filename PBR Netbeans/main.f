         program main        

         
         use specs
         use debug    
         !Include 'Initialize_AgeDistribution.f90'
         
         implicit none
         
	 integer ii, jj, kk, ll ! counters
	 integer xx, nn
         
	 real aa, bb
	 real start, finish
	 real N_best, CV_N, z_variate
	 real N_min
	 real x_min, x_max, step_CV
         real a_m
         real yy

	 character*72 ccon
	 integer ic
	 
10	 format(a10 f5.3 a10 f8.3)
20	 format(a10 a10) ! header for CV_N output file
30	 format(f10.3 f10.3) ! header for CV_N output file

	 CALL RANDOM_SEED() ! initialize seed for random number generator
		 
!		write(*,*) 'Calling cpu_time'
!		CALL test_cpu_time(start, finish, aa, bb)
         a_m = 3.
         write(*,*) "a_m = : ", a_m
         yy = debugz(a_m)
         write(*,*) "yy = : ", yy
         
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
	end
	   

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
	
