module calcs

  use random, only : random_normal   ! Routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : Random_module.f90
  use Declare_variables_module       ! Access to some global variables, like selectivity patterns
  implicit none
         
  contains
!====== +++ === === +++ === === +++ ===   
!###### +++ ### ### +++ ### ### +++ ###         
  subroutine pop_projection(f_rate, b_rate, N_age_old, N_age_new)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===     
! Project population one year into the future, given human caused mortality (f_rate), selectivity pattern etc.            
! Assumes access to some global variables, like selectivity patterns at age (via, e.g. "use Declare_variables_module"      
! Presently assumes b_sex_ratio is 50:50, if this is not true, then need to modify this code      
!====== +++ === === +++ === === +++ === 
        real(kind = 8), intent(in) :: f_rate                ! Human caused mortality rate -- input, calculated prior to projection
        real(kind = 8), intent(in) :: b_rate                ! The birth rate at the start of year t+1; calculated prior to projection
        real(kind = 8), intent(in) :: N_age_old(0:age_x)    ! Numbers at age at previous time step
        real(kind = 8), intent(out) :: N_age_new(0:age_x)   ! Numbers at age at next time step (ouput)
        integer(kind = 4) :: ii                             ! Counter for looping over ages
! Calculate calves before mortality
        N_age_new(0) = b_sex_ratio * b_rate * sum(N_age_old(:) * prop_mat_age(:))  

        do ii = 1, (age_x - 1)                                              ! NPR calculations from Age 1 to Age x-1
            N_age_new(ii) = N_age_old(ii - 1) * S_age(ii - 1)                       ! Subtract natural mortality 
            N_age_new(ii) = N_age_new(ii) * (1 - selectivity(ii - 1) * f_rate)      ! Subtract any human caused mortality 
        end do
! Calculate numbers in the plus-group        
        N_age_new(age_x) = N_age_old(age_x - 1) * S_age(age_x - 1)                  ! Natural mortality into plus-group
        N_age_new(age_x) = N_age_new(age_x) * (1 - selectivity(age_x - 1) * f_rate) ! Human caused mortality into plus-group         
        N_age_new(age_x) = N_age_new(age_x) + N_age_old(age_x) * S_age(age_x) * (1 - selectivity(age_x) * f_rate) ! Plus-group mortality

        return
    end subroutine pop_projection   
! *    
!====== +++ === === +++ === === +++ ===
!###### +++ ### ### +++ ### ### +++ ###       
  function assign_surv_yrs_stock(yr_max, n_stocks, surv_freq)
!###### +++ ### ### +++ ### ### +++ ###           
!====== +++ === === +++ === === +++ ===         
! Given the interval between surveys for each stock (surv_freq(1:n_stocks), 
!  Assign a value of 1 to the is_surv_yr(0:yr_max, 1_n_stocks) matrix, if the year (row) is
!  a survey year for that stock (each stock can have different survey intervals), 
!  or a zero to that element otherwise. 
!====== +++ === === +++ === === +++ === 
  ! TODO DECLARE VARIABLES, MAKE SURE FUNCTION DECLARED AS RETURNING A MATRIX
    integer(kind = 4):: ii, yr ! Counters
    integer(kind = 4), intent(in) :: yr_max
    integer(kind = 4), intent(in) :: n_stocks  
    integer(kind = 4), intent(in) :: surv_freq(1:n_stocks)
    integer(kind = 4) :: assign_surv_yrs_stock(0:yr_max, 1:n_stocks)
    
    assign_surv_yrs_stock = 0       ! Set all elements in this matrix to zero

    do ii = 1, n_stocks ! Debugging. Can move this test for survey year inside main loop -- or a function call
        assign_surv_yrs_stock(1, :) = 1 ! Year one is first survey year
        do yr = 2, yr_max
            assign_surv_yrs_stock(yr, ii) = mod(yr - 1, surv_freq(ii))   ! See if this is a survey year for this stock (subtract one from year because first survey defined to be in year one)
            if(assign_surv_yrs_stock(yr, ii) == 0) then
                assign_surv_yrs_stock(yr, ii) = 1
            else
                assign_surv_yrs_stock(yr, ii) = 0
            end if
        end do ! End loop over years
    end do     ! End loop over stocks
    return
  end function assign_surv_yrs_stock
!====== +++ === === +++ === === +++ === 
!###### +++ ### ### +++ ### ### +++ ###       
  real(kind = 8) function gen_survey_estimate(true_abundance, cv_n)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===         
! Given 'true' (operating model) abundance, generate a survey estimate of abundance, given survey sampling error (CV)
! The sampling error of the survey is assumed to be log-normally distributed, and the true abundance is assumed
!    to be equal to the expectation (mean) of that sampling error.      
! This function will use a standard normal random deviate, and is based on Eqn 3 of Wade(1998, Mar Mamm Sci)    
!====== +++ === === +++ === === +++ === 
        real(kind = 8), intent(in) :: true_abundance ! 'True' abundance from operating model 
        real(kind = 8), intent(in) :: cv_n           ! Coefficient of variation of survey sampling error
        real(kind = 8) :: z_variate                  ! Standard normal random variate 
        
        z_variate = random_normal(mean = 0., sd = 1.)
        
        gen_survey_estimate = z_variate
        return
    end function gen_survey_estimate
        
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

