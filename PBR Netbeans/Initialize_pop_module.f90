MODULE initialize_pop    
  use Declare_variables_module ! Access global variables 
  use eigen_module             ! Contains wrappers for LAPACK procedure DGEEV: Calculates eigenvalues/vectors of projection matrix
  
  implicit none                ! This rule now applies to each procedure below    
!-----------------------------------------------------------------------------------------
  contains        
!   * subroutine assign_par_vectors() *    
! Assign values to survival, proportion mature and selectivity at age vectors
!    
!   * subroutine calc_NPR_age() *    
! Does numbers-per-recruit calculations (NPR). 
! Used to calculate the unexploited and exploited age-structure in the first year of the projections                 
!    
!   * function Initial_F() *
! Given unexploited numbers per-recruit at equilibrium - apply initial human caused mortality rate (F_init)
! This function is called by Brent's method (file: BRENT.f90) to find f_init that results in a stable age structure, given 
!  the initial depletion level of the population relative to carrying capacity, i.e. the root of the objective function value
!  returned by Initial_F()    
!
!   * subroutine rescale_NPR() * 
! Rescale numbers at age per-recruit to initial numbers at age vector, given initial depletion 
!  Depletion is the abundance as a fraction of carrying capacity (e.g. if carrying capacity = 100, and abundance = 50, depletion = 0.50)  
!   
!   * function assign_area_stock_prop() * 
! Assign percentage of each stock (columns) to areas (rows)
!  
!   * function assign_transition_matrix() *
! Assign life history values to elements of population projection matrix. The projection matrix is used when solving for the 
!  juvenile survival rate that corresponds with a pre-specified maximum rate of population growth. The maximum rate of pop growth
!  is the eigenvalue of the projection matrix with birth rate at its maximum (lambda_max). Note: R_max = (Lambda_max - 1).   
!  
!   * function calc_lambda() *
! Given a value for juvenile survival calculate lambda_max_tmp, and return the squared difference between that and target lambda_max.
!  Brent's method is used to find the root of this objective function, the solution for juvenile survival.    
!  
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  real (kind = 8) function calc_recruit_F(NPR, b_eq, selex_a, surv_a)
!-----------------------------------------------------------------------------------------  
! Returns the number of recruits (age 0) given an initial exploitation rate, selectivity ogive etc.
! This is Eqn 16 of Punt, A.E. 1999. A full description of the BALEEN II Model ... JCRM 1 (suppl): 267-276
! Notes:  
! `NPR_oneplus` and `sum_NPR` are global in scope
! `NPR_mature` is declared in main program
!-----------------------------------------------------------------------------------------
    real (kind = 8), intent(in) :: selex_a(0:age_x), surv_a(0:age_x)
    real (kind = 8), intent(in) :: NPR(0:age_x) ! Numbers per recruit
!    real (kind = 8) :: recruit
    real (kind = 8) :: NPR_tmp(0:age_x), NPR_F0(0:age_x) 
    real (kind = 8) :: sum_NPR_1plus, sum_NPR_mat
    real (kind = 8) :: sum_NPR_1plus_F0, sum_NPR_mat_F0    
    real (kind = 8) :: NPR_sum_recd, NPR_sum_unrecd, NPR_sum_recd_F0, NPR_sum_unrecd_F0
    real (kind = 8) :: A
    real (kind = 8) :: b_eq ! Calculated in main program after calling NPR_age() = 1 / NPR_mature(F = 0)
    real (kind = 8) :: b_Finit_tmp 
    integer (kind = 4) :: Nsearch, ii
    real (kind = 8) :: ff
    real (kind = 8) :: N_0_Finit
    
! Initialize
    Nsearch = 500
    A = (b_max - b_eq) / b_eq ! b_eq calculated in main program given F = 0 (pre-exploitation equilibrium birth rate) 
    
    print *, "A parameter: "
    print *, A

! Caculate Numbers per female recruit at pre-exploitation equilibrium (F = 0)    
    Call calc_NPR_age(f_rate = 0.d0, N_recruits = b_sex_ratio, &
      N_age_tmp = NPR_F0, sum_1plus = sum_NPR_1plus_F0, sum_mature = sum_NPR_mat_F0, &
      NPR_sum_recd = NPR_sum_recd_F0, NPR_sum_unrecd = NPR_sum_unrecd_F0)

    print *, "ii; ff:"    
    
    do ii = 0, Nsearch
      
      ff = real(ii, kind = 8) / real(Nsearch, kind = 8) * 0.05
      
      Call calc_NPR_age(f_rate = ff, N_recruits = b_sex_ratio, &
        N_age_tmp = NPR_tmp, sum_1plus = sum_NPR_1plus, sum_mature = sum_NPR_mat, &
        NPR_sum_recd = NPR_sum_recd, NPR_sum_unrecd = NPR_sum_unrecd)
      
      b_Finit_tmp = 1.d0 / sum_NPR_mat
      
      N_0_Finit = (b_Finit_tmp / b_eq) - 1
      N_0_Finit = N_0_Finit / A
      N_0_Finit = 1.d0 - N_0_Finit
      N_0_Finit = N_0_Finit ** (1 / theta)
      
      if(N_0_Finit .gt. 0.d0) then 
        
      else
        N_0_Finit = 0.d0
      end if
      
!      print *, ii, ff, N_0_Finit
        
      if(dd_component .eq. "1") then
              
      else 
        print *, "ERROR in function calc_recruit_F(): `dd_component` not equal to '1'"
        print *, "Density dependence is currently only supported in terms of the age 1+ component"
        stop
      end if
      
    end do

! Debugging checking -- start with pre-exploitation equilibrium: F = 0.    
    Call calc_NPR_age(f_rate = 0.d0, N_recruits = b_sex_ratio, &
      N_age_tmp = NPR_tmp, sum_1plus = sum_NPR_1plus, sum_mature = sum_NPR_mat, &
      NPR_sum_recd = NPR_sum_recd, NPR_sum_unrecd = NPR_sum_unrecd)
      
    print *, "NPR_tmp"
    print *, NPR_tmp    
    print *, "sum_NPR_1plus"
    print *, sum_NPR_1plus            
    print *, "sum_NPR_mat"
    print *, sum_NPR_mat    
    print *, "NPR_sum_recd"
    print *, NPR_sum_recd        
    print *, "NPR_sum_unrecd"
    print *, NPR_sum_unrecd      
!    stop
    
    calc_recruit_F = 1.d0

! # This is apparently the function that need to code in Fortran -----------------
!    Recruit <- function(N, Apar, fec_eq, theta, K){
!      # Looks like Eqn 16 of Punt (1999)
!      # Apar = (f_max - f_eq) / f_eq
!      # Returns the number of recruits (age 0, indexed as age 1 in R code)
!      # Recruits are a function of fishing mortality / selectivity ogive 
!      InvF <- 1.0 / (fec_eq * N[6]) 
!      Top <- 1.0 - (InvF - 1) / Apar
!      if (Top > 0)
! # The "K/sumN[6]" bit corresponds to the density-dependent component    
! # So, sum(N[:]) needs to be same component as K (e.g. 1+ ; mature ; etc.)    
! # Rescales up from single recruit to total number of recruits, given K     
!       Recruit <- K / sum(N[6]) * (Top^(1/theta)) 
!      else
!       Recruit <- 0
!    }  
    
    return
  end function calc_recruit_F
  
!-----------------------------------------------------------------------------------------
  subroutine assign_par_vectors(a_r, a_m, a_t, age_x, S_adult, S_juv, S_age, selectivity_age, prop_mat_age)
!-----------------------------------------------------------------------------------------        
! Assign values to survival, proportion mature and selectivity at age vectors
!-----------------------------------------------------------------------------------------
    integer(kind = 4), intent(in) :: a_r    ! Age at recruitment to human caused mortality
    integer(kind = 4), intent(in) :: a_m    ! Age at first partuition
    integer(kind = 4), intent(in) :: a_t    ! Age at transition to adult survival
    integer(kind = 4), intent(in) :: age_x  ! Plus-group age        
    real(kind = 8), intent(in) :: S_adult
    real(kind = 8), intent(in) :: S_juv
    real(kind = 8), intent(out) :: S_age(0:age_x)        
    real(kind = 8), intent(out) :: selectivity_age(0:age_x)
    real(kind = 8), intent(out) :: prop_mat_age(0:age_x)
! Note vectorized style for assignments below
    S_age(0 : (a_t - 1) ) = S_juv           ! Assign juvenile survival rates to ages < a_t
    S_age( a_t : age_x) = S_adult           ! Assign adult survival rates to a_m <= ages = age_x                  
    prop_mat_age(0 : (a_m - 1)) = 0.0d0     ! Assign maturity-at-age vector (0 for immature ages, 1 for mature ages)
    prop_mat_age(a_m : age_x) = 1.0d0       ! Ages that have reached age of first partuition 
    selectivity_age(0 : (a_r - 1)) = 0.0d0  ! Assign selectivity at age
    selectivity_age(a_r : age_x) = 1.0d0    ! Assumes knife-edge selectivity. Could expand to include alternative functions.  
    return
  end subroutine assign_par_vectors   

!-----------------------------------------------------------------------------------------
  function assign_transition_matrix(a_m, a_t, age_x, b_rate, S_age, prop_mat_age) result(transition_matrix)
! Assign life history values to elements of population projection matrix 
!-----------------------------------------------------------------------------------------
    integer(kind = 4), intent(in) :: a_m, a_t, age_x
    integer(kind = 4) :: ii
    real(kind = 8), intent(in) :: b_rate
    real(kind = 8), intent(in) :: S_age(0:age_x)
    real(kind = 8), intent(in) :: prop_mat_age(0:age_x)
    real(kind = 8) :: transition_matrix(0:age_x, 0:age_x)

    transition_matrix(0, :) = b_rate * prop_mat_age     ! Assign birth rates to mature ages in first row of matrix

    if (a_t == age_x) then
      do ii = 0, (age_x - 1)
        transition_matrix(ii + 1, ii) = S_age(ii)       ! Assign juvenile survival rates
      end do        
    else
      do ii = 0, a_t  
        transition_matrix(ii + 1, ii) = S_age(ii)       ! Assign juvenile survival rates              
      end do
      do ii = (a_t + 1), (age_x - 1) 
        transition_matrix(ii+1, ii) = S_age(ii)         ! Assign adult survival rates              
      end do
    end if

    transition_matrix(age_x, age_x) = S_age(age_x)      ! Adult survival assumed for plus-group

    return
  end function assign_transition_matrix   

!-----------------------------------------------------------------------------------------
  subroutine calc_NPR_age(f_rate, N_recruits, N_age_tmp, sum_1plus, sum_mature, NPR_sum_recd, NPR_sum_unrecd)
!-----------------------------------------------------------------------------------------    
! Do numbers-per-recruit calculations (NPR). 
! Used to calculate the unexploited and exploited stable age-structure in the first year of the projections                 
!-----------------------------------------------------------------------------------------
    real(kind = 8), intent(in) :: f_rate                ! Human caused mortality rate
    real(kind = 8), intent(in) :: N_recruits            ! Female recruits that are age zero (female calves) for NPR calcs and Initial_F()
    real(kind = 8), intent(inout) :: N_age_tmp(0:age_x) ! Numbers at age vector   
    real(kind = 8), intent(out) :: sum_1plus            ! Size of 1+ component 
    real(kind = 8), intent(out) :: sum_mature           ! Size of mature component 
    real(kind = 8), intent(out) :: NPR_sum_recd, NPR_sum_unrecd 
    real(kind = 8) :: prop_N_age(0:age_x)               ! Proportion at age : sum(prop_NPR) = 1.0
    real(kind = 8) :: N_age_mat(0:age_x)                ! Numbers at age per that are mature
    real(kind = 8) :: N_age_immat(0:age_x)              ! Numbers at age per that are immature
    real(kind = 8) :: sum_N_age                         ! The sum of numbers across ages        
    integer(kind = 4) :: ii                             ! Counter for looping over ages
! Initialize to number of age zero recruits
    N_age_tmp = 0.0d0
    sum_1plus = 0.0d0
    sum_mature = 0.0d0

! Fraction of female offspring per recruit, assuming 50:50 sex ratio at birth
    N_age_tmp(0) = N_recruits    

! NPR calculations from Age 1 to Age x-1
    do ii = 1, (age_x - 1)                                              
        N_age_tmp(ii) = N_age_tmp(ii - 1) * S_age(ii - 1)                       ! Subtract natural mortality 
        N_age_tmp(ii) = N_age_tmp(ii) * (1 - selectivity(ii - 1) * f_rate)      ! Subtract any human caused mortality 
    end do
! Calculate numbers in the plus-group        
    N_age_tmp(age_x) = N_age_tmp(age_x - 1) * S_age(age_x - 1)                  ! Natural mortality into plus-group
    N_age_tmp(age_x) = N_age_tmp(age_x) * (1 - selectivity(age_x - 1) * f_rate) ! Human caused mortality into plus-group         
    N_age_tmp(age_x) = N_age_tmp(age_x) / (1 - S_age(age_x) * (1 - selectivity(age_x) * f_rate)) ! Plus-group mortality
! Note these are vectorized operations below 
    N_age_mat = N_age_tmp * prop_mat_age          ! Mature females
    N_age_immat = N_age_tmp * (1 - prop_mat_age)  ! Numbers immature at age per recruit
    sum_1plus = sum(N_age_tmp(1 : age_x))         ! Keep track of one_plus component (females)
    sum_mature = sum(N_age_mat)                   ! Sum mature females per female recruit
    sum_N_age = sum(N_age_tmp)                    ! Sum across ages in NPR vector
    prop_N_age = N_age_tmp / sum_N_age            ! Proportions at each age : sum(prop_N_age) = 1.0
    NPR_sum_recd = sum(N_age_tmp(a_r : age_x))    ! Numbers recruited per recruit
    NPR_sum_unrecd = sum(N_age_tmp(0 : a_r - 1))  ! Numbers un-recruited per recruit
    
    return
  end subroutine calc_NPR_age   
  
!-----------------------------------------------------------------------------------------
  function Initial_F(f_init) result(objf_f_init)
!-----------------------------------------------------------------------------------------
! Given numbers per-recruit at equilibrium - apply initial human caused mortality rate (F_init)
! Global Vars accessed:    
!        real(kind = 8) :: NPR_age                  ! Numbers per female recruit, see initialize_age_struc()      
!        real(kind = 8) :: NPR_age_tmp              ! Numbers at age conditioned on F_init
!        real(kind = 8) :: selectivity(age_x)       ! Selectivity at age
!        integer(kind = 4) :: a_r                   ! Age at recruitment (vulnerability) to human caused mortality
!        integer(kind = 4) :: age_x                 ! Plus group age (user input)
!        real(kind = 8) :: temp_mature          
!        real(kind = 8) :: temp_1plus
!        real(kind = 8) :: rec_init ! initial recruitment conditioned on Finit (in terms of females)
!        real(kind = 8) :: pred_rec
!        real(kind = 8) :: b_init
!        real(kind = 8) :: initial_oneplus
! b_sex_ratio : sex ratio at birth
!-----------------------------------------------------------------------------------------
    real(kind = 8), intent(inout) :: f_init ! Initial human caused mortality rate 
!    integer(kind = 4) :: ii, jj             ! Index counters for loops
    real(kind = 8) :: temp_mature           ! Size of mature component per female recruit given f_init = x_x
    real(kind = 8) :: temp_1plus            ! Size of 1+ component per female recruit given f_init = x_x
    real(kind = 8) :: foo                   ! calc_NPR_age returns sum_mature per recruit, which isn't necessary when finding f_init at the root of objf_f_init()
    real(kind = 8) :: rec_init              ! Initial recruitment conditioned on Finit (in terms of females)
    real(kind = 8) :: pred_rec              ! Predicted recruitment, used when finding f_init at the root of objf_f_init()
    real(kind = 8) :: initial_oneplus       ! Initial size of the 1+ component per female recruit (used in scaling)
    real(kind = 8) :: objf_f_init           ! Objective function to be minimized in finding root at f_init
    real(kind = 8) :: NPR_sum_recd, NPR_sum_unrecd
! Initialize
    temp_mature = 0.d0    
    temp_1plus = 0.d0     
    rec_init = 0.d0       
    pred_rec = 0.d0       
    initial_oneplus = 0.
    objf_f_init = -99.d0
! Do NPR calculations given f_init and N_recruits = NPR(age=0)
    Call calc_NPR_age(f_rate = f_init, N_recruits = b_sex_ratio, &  
        N_age_tmp = NPR_age_tmp, sum_1plus = temp_1plus, &
        sum_mature = temp_mature, NPR_sum_recd = NPR_sum_recd, NPR_sum_unrecd = NPR_sum_unrecd)        
!====== +++ === === +++ === === +++ === ! NPR_oneplus is global variable, based on NPR calculated at f_init = 0 (set by main program)
    rec_init = init_depl_i * NPR_oneplus / temp_1plus 	 ! Rescale to get the desired 1+ depletion conditioned on Finit 
    rec_init = rec_init / 2                              ! Divide by two to get female recruits conditioned on Finit      
!====== +++ === === +++ === === +++ === ! Calculate actual numbers per recruit conditioned on f_init       
    NPR_age_tmp(0) = rec_init
    Call calc_NPR_age(f_rate = f_init, N_recruits = rec_init, &
        N_age_tmp = NPR_age_tmp, sum_1plus = initial_oneplus, &
        sum_mature = foo, NPR_sum_recd = NPR_sum_recd, NPR_sum_unrecd = NPR_sum_unrecd) 
    pred_rec = b_init * temp_mature 
    pred_rec = pred_rec - (init_depl_i * NPR_oneplus / initial_oneplus)
    pred_rec = b_sex_ratio * pred_rec                               ! In terms of female component (probably not necessary for root finding)
    objf_f_init = pred_rec * pred_rec 
    return
  end function Initial_F    
! *    
!-----------------------------------------------------------------------------------------    
  subroutine rescale_NPR(k_1plus_tmp, initial_oneplus_tmp, N_age_unscaled, N_age_scaled)      
!-----------------------------------------------------------------------------------------  
! Rescale numbers at age per-recruit to initial numbers at age vector, given initial depletion 
!-----------------------------------------------------------------------------------------
    real(kind = 8), intent(in) :: k_1plus_tmp             ! Carrying capacity (in terms of age 1+) 
    real(kind = 8), intent(in) :: initial_oneplus_tmp     ! Initial one plus for stock i        
    real(kind = 8), intent(in) :: N_age_unscaled(0:age_x) ! Unscaled numbers at age in initial year. Note: age_x is global variable
    real(kind = 8), intent(out) :: N_age_scaled(0:age_x)  ! Unscaled numbers at age in initial year. Note: age_x is global variable              
    real(kind = 8) :: scale_pop                           ! Scalar to map numbers per recruit to initial numbers at age    

    scale_pop = b_sex_ratio * k_1plus_tmp / initial_oneplus_tmp ! Assume 50:50 birth rate (so same re-scaling factor used for each sex)  
    N_age_scaled = N_age_unscaled * scale_pop   ! This now in terms of females - so, just assign another vector equal to this one to initialize males

    return     
  end subroutine rescale_NPR

!-----------------------------------------------------------------------------------------
  function assign_area_stock_prop(p_a1_s1, p_a2_s1, p_a2_s2, p_a3_s2, p_a4_s2) result(area_stock_prop)
!-----------------------------------------------------------------------------------------  
! Assign percentage of each stock (columns) to areas (rows) of area x stock matrix
! TODO: All of this hard-coding is no good. Need to revise in future versions. -JRB    
!-----------------------------------------------------------------------------------------
    real(kind = 8) :: p_a1_s1       ! Percentage of stock 1 in area 1
    real(kind = 8) :: p_a2_s1       ! Percentage of stock 1 in area 1        
    real(kind = 8) :: p_a2_s2       ! Percentage of stock 1 in area 1
    real(kind = 8) :: p_a3_s2       ! Percentage of stock 1 in area 1
    real(kind = 8) :: p_a4_s2       ! Percentage of stock 1 in area 1        
    real(kind = 8), dimension(1:4, 1:n_stocks) :: area_stock_prop ! Hard-coded, could be improved: TODO

    if (n_stocks == 1) then
        area_stock_prop(1, 1) = p_a1_s1 ! Assign percentage of each stock (columns) to areas (rows)
        area_stock_prop(2, 1) = p_a2_s1
        area_stock_prop(3, 1) = 0.0d0
        area_stock_prop(4, 1) = 0.0d0
    else if (n_stocks == 2) then
        area_stock_prop(1, 1) = p_a1_s1 ! Assign percentage of each stock (columns) to areas (rows)
        area_stock_prop(2, 1) = p_a2_s1
        area_stock_prop(3, 1) = 0.0d0
        area_stock_prop(4, 1) = 0.0d0
        area_stock_prop(1, 2) = 0.0d0
        area_stock_prop(2, 2) = p_a2_s2
        area_stock_prop(3, 2) = p_a3_s2
        area_stock_prop(4, 2) = p_a4_s2
    else
        print *, "Error from assign_area_stock_prop(): "
        print *, "n_stocks > 2: Not supported currently."
        stop
    end if

    return
  end function assign_area_stock_prop ! End function

!-----------------------------------------------------------------------------------------
  function calc_lambda(s_juv_tmp) result(objf_lambda)
!-----------------------------------------------------------------------------------------  
! Given a value for juvenile survival calculate lambda_max_tmp, and return the squared difference between that and target lambda_max
!  Notes: 
!   - Holding other life history parameters (max birth rate, adult survival etc.) constant from the values in input.par file    
!   - Lambda_max is the maximum population growth rate (i.e. calculated assuming birth rate = b_max)
!   - PBR base-case scenarios for cetaceans (whales and dolphins) assume lambda_max = 1.04 (i.e. R_max = 0.04)
!   - Life history parameters (survival etc.) and R_max are read from input.par, and declared with global scope in Declare_variables_module.f90
!   - This function is called by Brent's method, to find the minimum of objf_lambda = (lambda_max_tmp - lambda_max)^2 
!   - Brent's method iterates over values for juvenile survival, until a solution is found corresponding with the target lambda_max (the root of objf_lambda)
!   - Eigen() found in Eigen_module.f90. 
!       Eigen() is a wrapper for DGEEV(), the LAPACK (Linear Algebra PACKage) procedure for eigenvalues/vectors from a non-symmetric matrix.
!       See also, www.netlib.org            
!-----------------------------------------------------------------------------------------
    real(kind = 8), intent(in) :: s_juv_tmp         ! Juvenile survival rate
    real(kind = 8) :: objf_lambda  ! Objective function for lambda. The squared difference between target and calculated
    real(kind = 8) :: transition_matrix_tmp(0 : age_x , 0 : age_x) ! Temp transition matrix
    real(kind = 8) :: S_age_tmp(0 : age_x)          ! Temporary vector of survival rates at age
    real(kind = 8) :: lambda_tmp                    ! Temporary value for dominant real eigenvalue of transition matrix with s_juv_tmp
    integer(kind = 4) :: aa, jj                     ! Counters

    objf_lambda = 0.0d0
    transition_matrix_tmp = 0.0d0                   ! Need to initialize matrix to zero (otherwise DGEEV is unreliable)

    Call assign_par_vectors(a_r, a_m, a_t, age_x, s_adult, s_juv_tmp, &  ! Slightly repetitious, but need to re-assign temporary juvenile survival rate
                S_age_tmp, selectivity, prop_mat_age)   

! fecundity_max is in units of female calves per female, e.g. 0.50 * birth rate (b_max) for eigen analysis of projection matrix
! Need to use fecundity instead of birth rate to get eigenvalue (lambda) in terms of population (males & female) growth rate                     
    transition_matrix_tmp = assign_transition_matrix(a_m, a_t, age_x, fecundity_max, S_age_tmp, prop_mat_age)                    

    call eigen(transition_matrix_tmp, (age_x + 1), lambda_tmp)    ! Calculate dominant real eigenvalue   

    objf_lambda = (lambda_tmp - (R_max + 1.0d0)) * (lambda_tmp - (R_max + 1.0d0)) ! Square of difference between target and proposed

    return
  end function calc_lambda

!-----------------------------------------------------------------------------------------
! THIS FUNCTION IS CURRENTLY BROKEN -- NOT GETTING SAME RESULTS AS LAPACK CALLS (ALSO TESTED IN R)   
!    real(kind = 8) function characteristic_eq(lambda_tmp)
!! Calculate the characteristic (polynomial) equation of the projection matrix
!! This is an attempt at the algebraic solution to the det(A-lambda*I) = 0
!! The equation used here is taken from the right hand side of Eqn 23 of Punt (1999)
!!   Punt, A.E. 1999. A Full Description of the Standard BALEEN II Model and Some Variants Thereof. J Cet Res Manage 1(Suppl): 267-276
!-----------------------------------------------------------------------------------------
!        real(kind = 8), intent(in) :: lambda_tmp    ! Proposed value for lambda
!        integer(kind = 4) :: ii,jj                  ! Counters
!        real(kind = 8) :: tmp1, tmp2, tmp3
!
!        tmp1 = 0.0d0    ! Initialize temporary variables
!        tmp2 = 0.0d0
!        tmp3 = 0.0d0
!    
!        print *, 
!        print *, "Hello from characteristic_eq: "
!        print *, "S_age(0:age_x-1): ", S_age(0:age_x-1)   
!        do ii = a_m, (age_x -1)
!            tmp1 = PRODUCT(S_age(0:ii-1)) ! SUM ( X, MASK = X .GT. 0.0)
!            print *, "ii: ", ii
!            print *, "tmp1: ", tmp1
!        end do
!   
!        characteristic_eq = tmp1    ! DEBUGGING
!        print *, "characteristic_eq: ", characteristic_eq
!        print *, 
!    !    characteristic_eq = b_max * (tmp2 + tmp3) ! Should equal 1.0 if lambda_tmp is the solution (root) of the equation
!
!        return
!    end function characteristic_eq
    
END MODULE initialize_pop

