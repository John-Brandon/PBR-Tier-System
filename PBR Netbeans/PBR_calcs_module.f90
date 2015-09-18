module calcs
! *
  use random, only : random_normal   ! Routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : Random_module.f90
  use Declare_variables_module       ! Access to some global variables, like selectivity patterns
  implicit none
! *         
  contains
! ----------------------------------------------------------------------------------------
!* subroutine pop_projection(f_rate, b_rate, N_age_old, N_age_new) *    
! ----------------------------------------------------------------------------------------
! Project population one year into the future, given human caused mortality (f_rate), selectivity pattern etc.            
! Assumes access to some global variables, like selectivity patterns at age (via, e.g. "use Declare_variables_module"      
! Presently assumes sex ratio at birth (b_sex_ratio) is 50:50, if this is not true, then need to modify this code       
!
! ----------------------------------------------------------------------------------------
!* gen_survey_estimate(true_abundance, cv) *
! ----------------------------------------------------------------------------------------
! Given 'true' (operating model) abundance, generate a survey estimate of abundance, given survey sampling error (CV)
! The sampling error of the survey is assumed to be log-normally distributed, and the true abundance is assumed
!    to be equal to the expectation (mean) of that sampling error.      
! This function will use a standard normal random deviate, and is based on Eqn 3 of Wade(1998, Mar Mamm Sci)    
!  
! ----------------------------------------------------------------------------------------
!* function calc_n_min(n_hat, cv, z_score) * 
! ----------------------------------------------------------------------------------------
! Function to calculate N_min, given (i) CV of abundance estimate (ii) N_best and (iii) standard normal z_score
! Uses equation (4) of Wade(1998 Mar Mamm Sci)
! Notes: 
! N_best is the expectation (mean) of a log-normal distribution, not the median.  
! CV is the coefficient of variation of the abundance estimate, gets transformed to the standard deviation in log-space
! Positive values for z-score return lower tail the way that Eqn 4 of Wade (1998) is derived. 
! For example, for the lower 2.5th percentile, z_score = 1.96 (not -1.96) 
! If you're more familiar with R, this function is equivalent to:
!  <- qlnorm(p = z_score, meanlog = log(n_best), sdlog = sqrt(log(1 + cv_n * cv_n)))  
!  
! ----------------------------------------------------------------------------------------
!* subroutine  weight_avg(n_yrs, N_hat, CV_N, N_avg, CV_avg) *
! ----------------------------------------------------------------------------------------
! Calculate a weighted average abundance estimate. 
! Average is weighted by the precision (inverse of the variance) of each estimate.
! Methods follow NMFS (2005) pp. 12-13 
!  NMFS. 2005. Revisions to Guidelines for Assessing Marine Mammal Stocks. 24 pp.
!   Available at: http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf  
!  
! ----------------------------------------------------------------------------------------
!* calc_n_mean(cv, n_hat, cv_mean, n_mean)
! ----------------------------------------------------------------------------------------
! Function to calculate the CV of weighted average. 
! This is similar to the `weight_avg()` subroutine, but differs in calculation of N_mean
! The assumption here for N_mean follows from the approach taken by,
!   Wade and DeMaster. 1999. Determining the optimum interval for abundance surveys
! That assumption is that the CV's are equal for the simulated survey estimates,
!   so they just used the arithmetic mean (N_mean) to plug into calculations of N_min
! I think the geometric mean might be a closer approximation, given assumption of log-normality? -jbrandon
!  
! ----------------------------------------------------------------------------------------
!* subroutine calc_n_hat(tier, pbr, n_hat, cv_n)   
! ----------------------------------------------------------------------------------------
! Calculate estimate of current abundance (and associated CV) based on data tier
!  
! ----------------------------------------------------------------------------------------
!* function assign_surv_yrs_stock(yr_max, n_stocks, surv_freq) *
! ----------------------------------------------------------------------------------------
! Not used. Assigns vector with dummy variables indicating if year is a survey year.  
!  
! ----------------------------------------------------------------------------------------      
  subroutine pop_projection(f_rate, b_rate, N_age_old, N_age_new)
! ----------------------------------------------------------------------------------------    
! Project population one year into the future, given human caused mortality (f_rate), selectivity pattern etc.            
! Assumes access to some global variables, like selectivity patterns at age (via, e.g. "use Declare_variables_module"      
! Presently assumes b_sex_ratio is 50:50, if this is not true, then need to modify this code      
! ----------------------------------------------------------------------------------------
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
! ----------------------------------------------------------------------------------------      
  real(kind = 8) function gen_survey_estimate(true_abundance, cv)
! ----------------------------------------------------------------------------------------        
! Given 'true' (operating model) abundance, generate a survey estimate of abundance, given survey sampling error (CV)
! The sampling error of the survey is assumed to be log-normally distributed, and the true abundance is assumed
!    to be equal to the expectation (mean) of that sampling error.      
! This function will use a standard normal random deviate, and is based on Eqn 3 of Wade(1998, Mar Mamm Sci)    
! ----------------------------------------------------------------------------------------
        real(kind = 8), intent(in) :: true_abundance ! 'True' abundance from operating model 
        real(kind = 8), intent(in) :: cv             ! Coefficient of variation of survey sampling error
        real(kind = 8) :: surv_tmp                   ! Temporary variable for calculations
        real(kind = 8) :: z_variate                  ! Standard normal random variate 
        
        z_variate = random_normal(mean = 0., sd = 1.) ! Generate standard normal deviate
        surv_tmp = cv * cv
        surv_tmp = log(1 + surv_tmp)
        surv_tmp = sqrt(surv_tmp)
        surv_tmp = z_variate * surv_tmp
        gen_survey_estimate = true_abundance / sqrt(1 + cv * cv)
        gen_survey_estimate = log(gen_survey_estimate)
        gen_survey_estimate = gen_survey_estimate + surv_tmp
        gen_survey_estimate = exp(gen_survey_estimate)
        
        return
    end function gen_survey_estimate
! *
! ----------------------------------------------------------------------------------------   
  real(kind = 8) function calc_n_min(n_hat, cv, z_score)
! ----------------------------------------------------------------------------------------
! Function to calculate N_min, given (i) CV of abundance estimate (ii) N_best and (iii) standard normal z_score
! Uses equation (4) of Wade(1998 Mar Mamm Sci)
! Notes: 
! N_best is the expectation (mean) of a log-normal distribution, not the median.  
! CV is the coefficient of variation of the abundance estimate, gets transformed to the standard deviation in log-space
! Positive values for z-score return lower tail the way that Eqn 4 of Wade (1998) is derived. 
! For example, for the lower 2.5th percentile, z_score = 1.96 (not -1.96) 
! If you're more familiar with R, this function is equivalent to:
!  <- qlnorm(p = z_score, meanlog = log(n_best), sdlog = sqrt(log(1 + cv_n * cv_n)))  
! ----------------------------------------------------------------------------------------
    real(kind = 8), intent(in) :: n_hat
    real(kind = 8), intent(in) :: cv
    real(kind = 8), intent(in) :: z_score
    
    calc_n_min = log(1 + cv * cv)   ! Start by calculating denominator
    calc_n_min = sqrt(calc_n_min)       ! calc_n_min is now the standard deviation in log-space
    calc_n_min = z_score * calc_n_min   ! z_score of -1.96 for lower 2.5th percentile; -0.842 for lower 20th percentile
    calc_n_min = exp(calc_n_min)
    calc_n_min = n_hat / calc_n_min ! divide N_best by denominator

    return
  end function calc_n_min
! *
! ----------------------------------------------------------------------------------------   
  subroutine weight_avg(n_yrs, N_hat, CV_N, N_avg, CV_avg)
! ----------------------------------------------------------------------------------------
! Calculate a weighted average abundance estimate. 
! Average is weighted by the precision (inverse of the variance) of each estimate.
! Methods follow NMFS (2005) pp. 12-13 
!  NMFS. 2005. Revisions to Guidelines for Assessing Marine Mammal Stocks. 24 pp.
!   Available at: http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf  
! ----------------------------------------------------------------------------------------
    integer(kind = 4), intent(in) :: n_yrs         ! Number of years to average over
    real(kind = 8), intent(in) :: N_hat(1:n_yrs)   ! Vector of abundance (point) estimates
    real(kind = 8), intent(in) :: CV_N(1:n_yrs)    ! Vector of corresponding CVs
    real(kind = 8) :: Var_N(1:n_yrs)               ! Variances of abundance estimates
    real(kind = 8) :: weights_tmp(1:n_yrs)         ! Inverse of the variances
    real(kind = 8) :: weights(1:n_yrs)             ! Relative precision
    real(kind = 8) :: weights_sq(1:n_yrs)          ! Square of the weights
    real(kind = 8) :: weight_tot                   ! Sum of weights
    real(kind = 8) :: N_weighted(1:n_yrs)          ! Weighted abundance estimates
    real(kind = 8) :: weighted_var(1:n_yrs)
    real(kind = 8) :: Var_avg
    real(kind = 8), intent(out) :: N_avg           ! Weighted average abundance estimate
    real(kind = 8), intent(out) :: CV_avg          ! CV of weighted average abundance estimate    
    
    Var_N = CV_N * N_hat   ! SD of abundance estimates 
    Var_N = Var_N * Var_N  ! Variance of abundance estimates
  
    weights_tmp = 1.d0 / Var_N    ! Calculate precision of each abundance estimate (note, vector operation)
    weight_tot = sum(weights_tmp) ! 
    weights = weights_tmp / weight_tot ! Calculate vector of weights
    
    N_weighted = weights * N_hat
    N_avg = sum(N_weighted)
    
    weights_sq = weights * weights
    weighted_var = weights_sq * Var_N
    Var_avg = sum(weighted_var)
    CV_avg = sqrt(Var_avg) / N_avg 
    
    return
  end subroutine weight_avg
! *  
! ---------------------------------------------------------------------------------------- 
  subroutine calc_n_mean(cv, n_hat, cv_mean, n_mean)
! ----------------------------------------------------------------------------------------
! Function to calculate the CV of weighted average. 
! This is similar to the `weight_avg()` subroutine, but differs in calculation of N_mean
! The assumption here for N_mean follows from the approach taken by,
!   Wade and DeMaster. 1999. Determining the optimum interval for abundance surveys
! That assumption is that the CV's are equal for the simulated survey estimates,
!   so they just used the arithmetic mean (N_mean) to plug into calculations of N_min
! I think the geometric mean might be a closer approximation, given assumption of log-normality? -jbrandon
! ----------------------------------------------------------------------------------------
    real(kind = 8), intent(in) :: cv(:)          ! vector of cv's
    real(kind = 8), intent(in) :: n_hat(:)       ! vector of abundance estimates (point estimates) 
    integer(kind = 4) :: n_surveys               ! number of survey estimates to average
    real(kind = 8), allocatable :: n_variance(:) ! vector of variances
    real(kind = 8) :: sum_variance               ! sum of the variances
    real(kind = 8), intent(out) :: cv_mean       ! cv of the averaged abundance estimate
    real(kind = 8), intent(out) :: n_mean        ! arithmetic mean of the abundance estimates
    
    n_surveys = size(cv)                         ! retrieve length of vector
    allocate(n_variance(1:n_surveys))            ! allocate size of variance vector
    
    n_variance = cv * n_hat                      ! standard deviation
    n_variance = n_variance * n_variance         ! variance
    
    sum_variance = sum(n_variance)               ! sum variances
    sum_variance = sum_variance / (n_surveys * n_surveys)
    
    n_mean = sum(n_hat) / n_surveys              ! arithmetic mean of the abundance estimates
    cv_mean = sqrt(sum_variance) / n_mean        ! cv of the averaged abundance estimate  
        
    return
  end subroutine calc_n_mean

! ----------------------------------------------------------------------------------------
  subroutine calc_n_hat(tier, yr, n_hat_yr_in, cv_n_in, n_hat_out, cv_n_out) 
! ----------------------------------------------------------------------------------------       
! Return an estimate of abundance and associated CV
! The approach to calculate N_hat will depend on data tier
! ----------------------------------------------------------------------------------------
    integer(kind = 4), intent(in) :: tier
    integer(kind = 4), intent(in) :: yr
    integer(kind = 4) :: yr_ii, yr_diff
    integer(kind = 4) :: n_surveys
    real(kind = 8), intent(in) :: n_hat_yr_in(0:yr_max), cv_n_in 
    real(kind = 8), intent(out) :: n_hat_out, cv_n_out 
    real(kind = 8) :: var_weighted 
    real(kind = 8) :: term1, term2, term3, term4
    real(kind = 8) :: n_hat
    real(kind = 8) :: n_hat_tmp(1:n_yrs_avg)
    real(kind = 8) :: cv_n_tmp(1:n_yrs_avg)
    real(kind = 8) :: var_n_tmp(1:n_yrs_avg)
    real(kind = 8), allocatable :: survey_ii(:)
    real(kind = 8), allocatable :: var_ii(:)
    real(kind = 8), allocatable :: tau_ii(:)
    real(kind = 8), allocatable :: weight_ii(:)

    term1 = 0.d0 ! initialize
    term2 = 0.d0
    term3 = 0.d0
    term4 = 0.d0
            
    select case (tier)
! this tier is for an index of abundance (not absolute abundance) ------------------------    
      case (1)    
! *         
! TODO         
! *                

! 
! The standard approach of using only last estimate of abundance as n_hat ----------------
      case (2) 
        
        n_hat_out = n_hat_yr_in(yr) 
        cv_n_out = cv_n_in   
!
! Simple weighted average tier (window to average = `n_yrs_avg`) -------------------------  
      case (3)       
! Don't have full set of years to average over early in time series
        if (yr .lt. n_yrs_avg) then 

! count estimates of abundance, everything after n_yrs_avg will be zero
! `pack()` is intrinsic function returning number of elements meeting certain condition           
          n_surveys = size(pack(n_hat_yr_in, n_hat_yr_in .gt. 0.d0))
          
          allocate(survey_ii(1:n_surveys)) ! allocate vectors
          allocate(var_ii(1:n_surveys))
          allocate(tau_ii(1:n_surveys))
          allocate(weight_ii(1:n_surveys))      
                       
          survey_ii(:) = pack(n_hat_yr_in, n_hat_yr_in .gt. 0.d0) 
          var_ii(:) = cv_n_in * survey_ii(:)        ! sd
          var_ii(:) = var_ii(:) * var_ii(:)         ! variance
          tau_ii(:) = 1.d0 / var_ii(:)              ! precision
          weight_ii(:) = tau_ii(:) / sum(tau_ii(:)) ! weights
                   
          n_hat_out = sum(weight_ii(:) * survey_ii(:))
          var_weighted = sum(weight_ii(:) * weight_ii(:) * var_ii(:))            
          cv_n_out = sqrt(var_weighted) / n_hat_out

        else ! yr .ge. n_yrs_avg

! Count number of abundance estimates
          
! Retrieve abundance estimates from the last `n_yrs_avg` number of years           
          n_hat_tmp(:) = n_hat_yr_in((yr - n_yrs_avg + 1):yr)
          n_surveys = size(pack(n_hat_tmp, n_hat_tmp .gt. 0.d0)) 

          allocate(survey_ii(1:n_surveys)) ! allocate vectors
          allocate(var_ii(1:n_surveys))
          allocate(tau_ii(1:n_surveys))
          allocate(weight_ii(1:n_surveys))      

          survey_ii(:) = pack(n_hat_tmp, n_hat_tmp .gt. 0.d0) 
          var_ii(:) = cv_n_in * survey_ii(:)        ! sd
          var_ii(:) = var_ii(:) * var_ii(:)         ! variance
          tau_ii(:) = 1.d0 / var_ii(:)              ! precision
          weight_ii(:) = tau_ii(:) / sum(tau_ii(:)) ! weights		      

          n_hat_out = sum(weight_ii(:) * survey_ii(:))

          var_weighted = sum(weight_ii(:) * weight_ii(:) * var_ii(:))

          cv_n_out = sqrt(var_weighted) / n_hat_out
                 
        end if

! Most data rich tier considered at present -----------------------------------        
      case (4) 

! Method to set PBR to time weighted-average follows IWC AWMP method for SLA's      
        if (yr .lt. n_yrs_avg) then                     

! Count estimates of abundance, everything after n_yrs_avg will be zero
          n_surveys = size(pack(n_hat_yr_in, n_hat_yr_in .gt. 0.d0))
          
          allocate(survey_ii(1:n_surveys)) ! allocate vectors
          allocate(var_ii(1:n_surveys))
          allocate(tau_ii(1:n_surveys))
          allocate(weight_ii(1:n_surveys))          

          do yr_ii = yr, 1, -1 ! looping backwards through years
          
            yr_diff = yr - yr_ii
            
            if (n_hat_yr_in(yr_ii) .gt. 0.d0) then

              term1 = term1 + 0.9 ** yr_diff * log(n_hat_yr_in(yr_ii)) / cv_n_in ** 2.0
              term2 = term2 + 0.9 ** yr_diff / cv_n_in ** 2.0
              term3 = term3 + 0.9 ** (2 * yr_diff) / cv_n_in ** 2.0
              term4 = term4 + 0.9 ** yr_diff / cv_n_in ** 2.0
            
            end if

          end do       
                                           
          n_hat_out = exp(term1 / term2)
          cv_n_out = sqrt(term3) / term4
                                                
        else ! yr .ge. n_yrs_avg        
          n_surveys = 0        
          do yr_ii = yr, (yr - n_yrs_avg + 1), -1 ! looping backwards through years
            yr_diff = yr - yr_ii
            if (n_hat_yr_in(yr_ii) .gt. 0.d0) then
              n_surveys = n_surveys + 1 
              term1 = term1 + 0.9 ** yr_diff * log(n_hat_yr_in(yr_ii)) / cv_n_in ** 2.0
              term2 = term2 + 0.9 ** yr_diff / cv_n_in ** 2.0
              term3 = term3 + 0.9 ** (2 * yr_diff) / cv_n_in ** 2.0
              term4 = term4 + 0.9 ** yr_diff / cv_n_in ** 2.0
            end if
          end do
          
          n_hat_out = exp(term1 / term2)
          cv_n_out = sqrt(term3) / term4
          
        end if
        
      case default ! Should never get here
        
        print *, "ERROR IN TIER: SHOULD BE EITHER: 0, 1, 2 or 3"
        
    end select 
              
    return
  end subroutine calc_n_hat

! ----------------------------------------------------------------------------------------
  subroutine avg_N(n_hat, cv_n, n_avg, cv_avg)
! ----------------------------------------------------------------------------------------
! Calculate a weighted average abundance estimate. 
! Average is weighted by precision (inverse of the variance) of each estimate.
! Methods follow NMFS (2005) pp. 12-13: 
!  NMFS. 2005. Revisions to Guidelines for Assessing Marine Mammal Stocks. 
!  Available at: http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf  
! ----------------------------------------------------------------------------------------  
    real(kind = 8), allocatable :: var_n(:)
    real(kind = 8), allocatable :: n_weighted(:)
    real(kind = 8) :: var_avg    
    real(kind = 8), intent(out) :: n_avg, cv_avg
    real(kind = 8) :: weights_tmp, weight_tot, weights, weights_sq, weighted_var
    real(kind = 8), intent(in) :: n_hat(:), cv_n(:)
    real(kind = 8), allocatable :: n_hat_tmp(:), cv_n_tmp(:)    
    
! debugging
    n_hat_tmp = [100.d0, 0.d0, 0.d0, 150.d0, 0.d0, 0.d0, 0.d0, 200.d0]  
    cv_n_tmp =  0.20d0
    
    print *, "Hello from `avg_N()`"
    print *, "n_hat_tmp: ", n_hat_tmp
    print *, "cv_n_tmp : ", cv_n_tmp
    stop    

  end subroutine avg_N
  
!! ----------------------------------------------------------------------------------------  
!  avg_N = function(N_hat, CV_N){
!#
!#  Calculate a weighted average abundance estimate. 
!#  Average is weighted by precision (inverse of the variance) of each estimate.
!#  Methods follow NMFS (2005) pp. 12-13 
!#   NMFS. 2005. Revisions to Guidelines for Assessing Marine Mammal Stocks. 
!#    Available at: http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf  
!! ----------------------------------------------------------------------------------------  
!  Var_N = NULL; N_weighted = NULL; N_avg = NULL
!  weights_tmp = NULL; weight_tot = NULL; weights = NULL
!  weights_sq = NULL; weighted_var = NULL
!  Var_avg = NULL; CV_avg = NULL
!  
!# DEBUGGING
!  N_hat = c(100, 200); CV_N = 0.30
!  
!  Var_N = CV_N * N_hat   # SD of abundance estimates 
!  Var_N = Var_N * Var_N  # Variance of abundance estimates
!
!# Calculate precision of each abundance estimate (note, vector operation)  
!  weights_tmp = 1. / Var_N    
!  weight_tot = sum(weights_tmp)      #  
!  weights = weights_tmp / weight_tot # Calculate vector of weights
!  
!  N_weighted = weights * N_hat
!  N_avg = sum(N_weighted)
!  
!  weights_sq = weights * weights
!  weighted_var = weights_sq * Var_N
!  Var_avg = sum(weighted_var)
!  CV_avg = sqrt(Var_avg) / N_avg 
!  return(data.frame(N_avg = N_avg, CV_avg = CV_avg))
!}
  
! *
! ----------------------------------------------------------------------------------------
!  function assign_surv_yrs_stock(yr_max, n_stocks, surv_freq)
! ----------------------------------------------------------------------------------------
!! Given the interval between surveys for each stock (surv_freq(1:n_stocks), 
!!  Assign a value of 1 to the is_surv_yr(0:yr_max, 1_n_stocks) matrix, if the year (row) is
!!  a survey year for that stock (each stock can have different survey intervals), 
!!  or a zero to that element otherwise. 
! ----------------------------------------------------------------------------------------
!  ! TODO DECLARE VARIABLES, MAKE SURE FUNCTION DECLARED AS RETURNING A MATRIX
!    integer(kind = 4):: ii, yr ! Counters
!    integer(kind = 4), intent(in) :: yr_max
!    integer(kind = 4), intent(in) :: n_stocks  
!    integer(kind = 4), intent(in) :: surv_freq(1:n_stocks)
!    integer(kind = 4) :: assign_surv_yrs_stock(0:yr_max, 1:n_stocks)
!    
!    assign_surv_yrs_stock = 0       ! Set all elements in this matrix to zero
!
!    do ii = 1, n_stocks ! Debugging. Can move this test for survey year inside main loop -- or a function call
!        assign_surv_yrs_stock(1, :) = 1 ! Year one is first survey year
!        do yr = 2, yr_max
!            assign_surv_yrs_stock(yr, ii) = mod(yr - 1, surv_freq(ii))   ! See if this is a survey year for this stock (subtract one from year because first survey defined to be in year one)
!            if(assign_surv_yrs_stock(yr, ii) == 0) then
!                assign_surv_yrs_stock(yr, ii) = 1
!            else
!                assign_surv_yrs_stock(yr, ii) = 0
!            end if
!        end do ! End loop over years
!    end do     ! End loop over stocks
!    return
!  end function assign_surv_yrs_stock
! *    
end module calcs

