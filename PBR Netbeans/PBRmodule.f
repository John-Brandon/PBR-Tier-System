      module calcs

         !implicit none
         
	 contains
! *************************************************************************      
! Function to calculate N_min, given (i) CV of abundance estimate (ii) N_best and (iii) standard normal variate, z
! Uses equation (4) of Wade(1998 Mar Mamm Sci)
! Note that N_best is the expectation (mean) of a log-normal distribution, not the median.        
! *************************************************************************            
	REAL(kind = 8) FUNCTION CALC_N_min(N_best, CV_N, z_variate)

	 real(kind = 8) :: CV_N, N_best, z_variate ! Argument declarations	

	 CALC_N_min = log(1 + CV_N * CV_N) ! start by calculating denominator
	 CALC_N_min = sqrt(CALC_N_min)
	 CALC_N_min = z_variate * CALC_N_min
	 CALC_N_min = exp(CALC_N_min)
	 CALC_N_min = N_best / CALC_N_min ! divide N_best by denominator
	 
	return
	end function CALC_N_min
	   
      end module calcs