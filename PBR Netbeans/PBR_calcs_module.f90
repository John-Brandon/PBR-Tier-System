module calcs
! *
  use random, only : random_normal   ! Routines for psuedo random number generators (RNG) -- only using random_normal() function at this stage : Random_module.f90
  use Declare_variables_module       ! Access to some global variables, like selectivity patterns
  implicit none
! *         
  contains
!====== +++ === === +++ === === +++ ===
!* subroutine pop_projection(f_rate, b_rate, N_age_old, N_age_new) *    
!====== +++ === === +++ === === +++ ===
! Project population one year into the future, given human caused mortality (f_rate), selectivity pattern etc.            
! Assumes access to some global variables, like selectivity patterns at age (via, e.g. "use Declare_variables_module"      
! Presently assumes sex ratio at birth (b_sex_ratio) is 50:50, if this is not true, then need to modify this code       
!
!====== +++ === === +++ === === +++ ===
!* gen_survey_estimate(true_abundance, cv) *
!====== +++ === === +++ === === +++ ===
! Given 'true' (operating model) abundance, generate a survey estimate of abundance, given survey sampling error (CV)
! The sampling error of the survey is assumed to be log-normally distributed, and the true abundance is assumed
!    to be equal to the expectation (mean) of that sampling error.      
! This function will use a standard normal random deviate, and is based on Eqn 3 of Wade(1998, Mar Mamm Sci)    
!  
!====== +++ === === +++ === === +++ ===
!* function calc_n_min(n_hat, cv, z_score) * 
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
!  
!====== +++ === === +++ === === +++ ===
!* subroutine  weight_avg(n_yrs, N_hat, CV_N, N_avg, CV_avg) *
!====== +++ === === +++ === === +++ ===
! Calculate a weighted average abundance estimate. 
! Average is weighted by the precision (inverse of the variance) of each estimate.
! Methods follow NMFS (2005) pp. 12-13 
!  NMFS. 2005. Revisions to Guidelines for Assessing Marine Mammal Stocks. 24 pp.
!   Available at: http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf  
!  
!====== +++ === === +++ === === +++ ===  
!* calc_n_mean(cv, n_hat, cv_mean, n_mean)
!====== +++ === === +++ === === +++ ===  
! Function to calculate the CV of weighted average. 
! This is similar to the `weight_avg()` subroutine, but differs in calculation of N_mean
! The assumption here for N_mean follows from the approach taken by,
!   Wade and DeMaster. 1999. Determining the optimum interval for abundance surveys
! That assumption is that the CV's are equal for the simulated survey estimates,
!   so they just used the arithmetic mean (N_mean) to plug into calculations of N_min
! I think the geometric mean might be a closer approximation, given assumption of log-normality? -jbrandon
!  
!====== +++ === === +++ === === +++ === 
!* subroutine calc_pbr(tier, pbr, n_hat, cv_n)   
!====== +++ === === +++ === === +++ ===    
! Calculate PBR based on data tier
!  
!====== +++ === === +++ === === +++ ===      
!* function assign_surv_yrs_stock(yr_max, n_stocks, surv_freq) *
!====== +++ === === +++ === === +++ ===      
! Not used. Assigns vector with dummy variables indicating if year is a survey year.  
!  
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
  real(kind = 8) function gen_survey_estimate(true_abundance, cv)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===         
! Given 'true' (operating model) abundance, generate a survey estimate of abundance, given survey sampling error (CV)
! The sampling error of the survey is assumed to be log-normally distributed, and the true abundance is assumed
!    to be equal to the expectation (mean) of that sampling error.      
! This function will use a standard normal random deviate, and is based on Eqn 3 of Wade(1998, Mar Mamm Sci)    
!====== +++ === === +++ === === +++ === 
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
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===    
  real(kind = 8) function calc_n_min(n_hat, cv, z_score)
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
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===    
  subroutine weight_avg(n_yrs, N_hat, CV_N, N_avg, CV_avg)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ === 
! Calculate a weighted average abundance estimate. 
! Average is weighted by the precision (inverse of the variance) of each estimate.
! Methods follow NMFS (2005) pp. 12-13 
!  NMFS. 2005. Revisions to Guidelines for Assessing Marine Mammal Stocks. 24 pp.
!   Available at: http://www.nmfs.noaa.gov/pr/pdfs/sars/gamms2005.pdf  
!====== +++ === === +++ === === +++ === 
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
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===    
  subroutine calc_n_mean(cv, n_hat, cv_mean, n_mean)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ === 
! Function to calculate the CV of weighted average. 
! This is similar to the `weight_avg()` subroutine, but differs in calculation of N_mean
! The assumption here for N_mean follows from the approach taken by,
!   Wade and DeMaster. 1999. Determining the optimum interval for abundance surveys
! That assumption is that the CV's are equal for the simulated survey estimates,
!   so they just used the arithmetic mean (N_mean) to plug into calculations of N_min
! I think the geometric mean might be a closer approximation, given assumption of log-normality? -jbrandon
!====== +++ === === +++ === === +++ === 
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
! *    
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===    
  subroutine calc_pbr(tier, pbr, n_hat, cv_n) 
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===        
! Calculate PBR based on data tier
!
!====== +++ === === +++ === === +++ ===     
    integer(kind = 4), intent(in) :: tier
    real(kind = 8), intent(in) :: n_hat ! TODO: make this an allocatable array, or else split this subroutine into multiple procedures -- one for each tier   
    real(kind = 8), intent(out) :: pbr
    real(kind = 8), intent(in) :: cv_n  ! TODO: make this an allocatable array, or else split this subroutine into multiple procedures -- one for each tier
    real(kind = 8) :: term1, term2, term3, term4
    real(kind = 8) :: n_min
    
    if (tier.eq.0) then
! this tier is for an index of abundance (not absolute abundance)
! TODO         
    else if (tier.eq.1) then ! standard pbr tier
      n_min = calc_n_min(n_hat = n_hat, cv = cv_n, z_score = -0.842d0)  
      pbr = 0.50 * r_max * n_min * F_r(1) ! TODO: input.par has two values for F_r(1:2), but this is not how model is set-up
    else if (tier.eq.2) then
! simple weighted average tier (window to average = 8 yrs?)
      
    else if (tier.eq.3) then ! most data rich tier
! Method to set PBR to time weighted-average follows IWC AWMP method for SLA's      
! Alternative (weighted) method (IYR is the year for which an estimate is needed)
      term1 = 0.d0
      term2 = 0.d0
      term3 = 0.d0
      term4 = 0.d0
!      DO IY = IYR-1,-39,-1 ! From (IYR-1) to -39, by -1
!        YrDiff = (IYR-IY)
!        if (CVX(IY).GT.0) then ! 0.90 is the down-weighting factor
!          term1 = term1 + 0.9 ** YrDiff * log(SIGHT(IY)) / CVX(IY) ** 2.0
!          term2 = term2 + 0.9 ** YrDiff / CVX(IY) ** 2.0
!          term3 = term3 + 0.9 ** (2 * YrDiff) / CVX(IY) ** 2.0
!          term4 = term4 + 0.9 ** YrDiff / CVX(IY) ** 2.0
!        end if
!      end do
!      ESTIMATE = EXP(TERM1/TERM2)
!      CVFINAL = SQRT(TERM3)/TERM4
!      SLA = 0.02 * ESTIMATE * EXP(-1.645 * CVFINAL) * IQUOTA ! Replace 0.02 with 0.5 * R_max 
!      CATCHQ = MIN(NEED, SLA)
    end if
    
    return
  end subroutine calc_pbr
! *
!====== +++ === === +++ === === +++ ===
!###### +++ ### ### +++ ### ### +++ ###       
!  function assign_surv_yrs_stock(yr_max, n_stocks, surv_freq)
!!###### +++ ### ### +++ ### ### +++ ###           
!!====== +++ === === +++ === === +++ ===         
!! Given the interval between surveys for each stock (surv_freq(1:n_stocks), 
!!  Assign a value of 1 to the is_surv_yr(0:yr_max, 1_n_stocks) matrix, if the year (row) is
!!  a survey year for that stock (each stock can have different survey intervals), 
!!  or a zero to that element otherwise. 
!!====== +++ === === +++ === === +++ === 
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

