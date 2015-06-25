module calcs
! *
  use random, only : random_normal   ! Routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : Random_module.f90
  use Declare_variables_module       ! Access to some global variables, like selectivity patterns
  implicit none
! *         
  contains
! *  
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
! *  
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
        real(kind = 8) :: surv_tmp                   ! Temporary variable for calculations
        real(kind = 8) :: z_variate                  ! Standard normal random variate 
        
        z_variate = random_normal(mean = 0., sd = 1.) ! Generate standard normal deviate
        surv_tmp = cv_n * cv_n
        surv_tmp = log(1 + surv_tmp)
        surv_tmp = sqrt(surv_tmp)
        surv_tmp = z_variate * surv_tmp
        gen_survey_estimate = true_abundance / sqrt(1 + cv_n * cv_n)
        gen_survey_estimate = log(gen_survey_estimate)
        gen_survey_estimate = gen_survey_estimate + surv_tmp
        gen_survey_estimate = exp(gen_survey_estimate)
        
        return
    end function gen_survey_estimate
! *
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===    
  real(kind = 8) function calc_n_min(n_best, cv_n, z_score)
!   real(kind = 8) function calc_nn_min(foo)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ === 
! Function to calculate N_min, given (i) CV of abundance estimate (ii) N_best and (iii) standard normal z_score
! Uses equation (4) of Wade(1998 Mar Mamm Sci)
! Notes: 
! N_best is the expectation (mean) of a log-normal distribution, not the median.  
! CV is the coefficient of variation of the abundance estimate, gets transformed to the standard deviation in log-space
! Positive values for z-score return lower tail the way that Eqn 4 of Wade (1998) is derived. 
! For example, for the lower 2.5th percentile, z_score = 1.96 (not -1.96) 
! If you're more familiar with R, this function is equivalent to:
!  <- qlnorm(p = z_score, meanlog = log(n_best), sdlog = sqrt(log(1 + cv_n * cv_n)))  
!====== +++ === === +++ === === +++ === 
    real(kind = 8), intent(in) :: n_best
    real(kind = 8), intent(in) :: cv_n
    real(kind = 8), intent(in) :: z_score
!   real(kind = 8), intent(in) :: foo
    
    print *, "Hello from calc_n_min"
    print *, "n_best: ", n_best
    print *, "cv_n: ", cv_n
    print *, "z_score: ", z_score
    
    calc_n_min = log(1 + cv_n * cv_n)   ! Start by calculating denominator
    calc_n_min = sqrt(calc_n_min)       ! calc_n_min is now the standard deviation in log-space
    calc_n_min = z_score * calc_n_min   ! z_score of -1.96 for lower 2.5th percentile; -0.842 for lower 20th percentile
    calc_n_min = exp(calc_n_min)
    calc_n_min = n_best / calc_n_min ! divide N_best by denominator

    return
    end function calc_n_min
! *
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===    
end module calcs

