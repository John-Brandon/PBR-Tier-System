MODULE initialize_pop    
    use Declare_variables_module ! Access global variables      
    implicit none   ! Declared at this point in the module, this rule now applies to each procedure below    
    contains    
!###### +++ ### ### +++ ### ### +++ ###     
!====== +++ === === +++ === === +++ ===     
!###### +++ ### ### +++ ### ### +++ ### 
    subroutine calc_NPR_age(f_rate, N_recruits, N_age, sum_1plus, sum_mature)
! Does numbers-per-recruit calculations (NPR). 
! Used to calculate the unexploited and exploited age-structure in the first year of the projections        
!====== +++ === === +++ === === +++ === 
        real(kind = 8), intent(in) :: f_rate                ! Human caused mortality rate
        real(kind = 8), intent(in) :: N_recruits            ! Female recruits that are age zero (female calves) for NPR calcs and Initial_F()
        real(kind = 8), intent(out) :: N_age(0:age_x)       ! Numbers at age vector   
        real(kind = 8), intent(out) :: sum_1plus            ! Size of 1+ component 
        real(kind = 8), intent(out) :: sum_mature           ! Size of mature component 
        real(kind = 8) :: prop_N_age(0:age_x)               ! Proportion at age : sum(prop_NPR) = 1.0
        real(kind = 8) :: N_age_mat(0:age_x)                ! Numbers at age per that are mature
        real(kind = 8) :: N_age_immat(0:age_x)              ! Numbers at age per that are immature
        real(kind = 8) :: sum_N_age                         ! The sum of numbers across ages        
        integer(kind = 4) :: kk                             ! Counter for looping over ages
! If this is the first year, initialize to numbers per recruit
        N_age(0) = N_recruits                               ! Fraction of female offspring per recruit, assuming 50:50 sex ratio at birth
! Loop over ages 
        do kk = 1, (age_x - 1)                                              ! NPR calculations from Age 1 to Age x-1
            N_age(kk) = N_age(kk - 1) * S_age(kk - 1)                       ! Calculate numbers-at-age per recruit
            N_age(kk) = N_age(kk) * (1 - selectivity(kk - 1) * f_rate)      ! Subtract any human caused mortality 
        end do
! Calculate numbers in the plus-group        
        N_age(age_x) = N_age(age_x - 1) * S_age(age_x - 1)                  ! Natural mortality into plus-group
        N_age(age_x) = N_age(age_x) * (1 - selectivity(age_x - 1) * f_rate) ! Human caused mortality into plus-group         
        N_age(age_x) = N_age(age_x) / (1 - S_age(age_x) * (1 - selectivity(age_x) * f_rate)) ! Plus-group mortality
! Note these are vectorized operations below 
        N_age_mat = N_age * prop_mat_age          ! Mature females
        N_age_immat = N_age * (1 - prop_mat_age)  ! Numbers immature at age per recruit
        sum_1plus = sum(N_age(1 : age_x))         ! Keep track of one_plus component (females)
        sum_mature = sum(N_age_mat)               ! Sum mature females per female recruit
        sum_N_age = sum(N_age)                    ! Sum across ages in NPR vector
        prop_N_age = N_age / sum_N_age            ! Proportions at each age : sum(prop_N_age) = 1.0
        return
    end subroutine calc_NPR_age   
!###### +++ ### ### +++ ### ### +++ ###     
!====== +++ === === +++ === === +++ ===     
!###### +++ ### ### +++ ### ### +++ ###     
    real(kind = 8) function Initial_F(f_init) result(objf_f_init)
!====== +++ === === +++ === === +++ === 
! Given numbers per-recruit at equilibrium - apply initial human caused mortality rate (F_init)
! Global Vars    
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
!====== +++ === === +++ === === +++ === 
        real(kind = 8), intent(inout) :: f_init ! Initial human caused mortality rate 
        integer(kind = 4) :: ii, jj             ! Index counters for loops
        real(kind = 8) :: temp_mature           ! Size of mature component per female recruit given f_init = x_x
        real(kind = 8) :: temp_1plus            ! Size of 1+ component per female recruit given f_init = x_x
        real(kind = 8) :: foo                   ! calc_NPR_age returns sum_mature per recruit, which isn't necessary when finding f_init at the root of objf_f_init()
        real(kind = 8) :: rec_init              ! Initial recruitment conditioned on Finit (in terms of females)
        real(kind = 8) :: pred_rec              ! Predicted recruitment, used when finding f_init at the root of objf_f_init()
        real(kind = 8) :: initial_oneplus       ! Initial size of the 1+ component per female recruit (used in scaling)
        real(kind = 8) :: objf_f_init	! Objective function to be minimized in finding root at f_init
!====== +++ === === +++ === === +++ === ! Initialize
        temp_mature = 0.    
        temp_1plus = 0.     
        rec_init = 0.       
        pred_rec = 0.       
        initial_oneplus = 0.
        objf_f_init = -99.d0
!====== +++ === === +++ === === +++ === ! Initialize        
        Call calc_NPR_age(f_rate = f_init, N_recruits = b_sex_ratio, &  ! Do NPR calculations given f_init and N_recruits = NPR(age=0)
            N_age = NPR_age_tmp, sum_1plus = temp_1plus, &
            sum_mature = temp_mature)        
!====== +++ === === +++ === === +++ === ! Note: NPR_oneplus is global variable, based on NPR calculated at f_init = 0 (set by main program)
        rec_init = init_depl_i * NPR_oneplus / temp_1plus 	 ! Rescale to get the desired 1+ depletion conditioned on Finit 
        rec_init = rec_init / 2                                  ! Divide by two to get female recruits conditioned on Finit      
!====== +++ === === +++ === === +++ === ! Calculate actual numbers per recruit conditioned on f_init       
        NPR_age_tmp(0) = rec_init
        Call calc_NPR_age(f_rate = f_init, N_recruits = rec_init, &
            N_age = NPR_age_tmp, sum_1plus = initial_oneplus, &
            sum_mature = foo) 
        pred_rec = b_init * temp_mature 
        pred_rec = pred_rec - (init_depl_i * NPR_oneplus / initial_oneplus)
        pred_rec = b_sex_ratio * pred_rec                               ! In terms of female component
        objf_f_init = pred_rec * pred_rec 
        return
    end function Initial_F    
!###### +++ ### ### +++ ### ### +++ ###     
!====== +++ === === +++ === === +++ ===     
!###### +++ ### ### +++ ### ### +++ ### 
    subroutine assign_par_vectors(a_r, a_m, age_x, S_adult, S_juv, S_age, selectivity_age, prop_mat_age)
! Assign values to survival, proportion mature and selectivity at age vectors
!====== +++ === === +++ === === +++ === 
        integer(kind = 4), intent(in) :: a_r
        integer(kind = 4), intent(in) :: a_m
        integer(kind = 4), intent(in) :: age_x        
        real(kind = 8), intent(in) :: S_adult
        real(kind = 8), intent(in) :: S_juv
        real(kind = 8), intent(out) :: S_age(0:age_x)        
        real(kind = 8), intent(out) :: selectivity_age(0:age_x)
        real(kind = 8), intent(out) :: prop_mat_age(0:age_x)
! Note vectorized style for assignments below
        S_age(0 : (a_m - 1)) = S_juv            ! Assign juvenile survival rates to ages < a_m
        S_age(a_m : age_x) = S_adult            ! Assign adult survival rates to a_m <= ages = age_x                  
        prop_mat_age(0 : (a_m - 1)) = 0.0       ! Assign maturity-at-age vector (0 for immature ages, 1 for mature ages)
        prop_mat_age(a_m : age_x) = 1.0
        selectivity_age(0 : (a_r - 1)) = 0.0    ! Assign selectivity at age
        selectivity_age(a_r : age_x) = 1.0        
        return
    end subroutine assign_par_vectors    
!###### +++ ### ### +++ ### ### +++ ###     
!====== +++ === === +++ === === +++ ===     
!###### +++ ### ### +++ ### ### +++ ### 
!    subroutine rescale_NPR(k_1plus_i, init_depl_i, initial_oneplus_i) 
!====== +++ === === +++ === === +++ === 
!! Rescale numbers at age per-recruit to initial numbers at age vector, given initial depletion 
!        implicit none
!        
!        real(kind = 8) :: k_1plus_i         ! Carrying capacity (in terms of age 1+) for stock i
!        real(kind = 8) :: init_depl_i       ! Initial depletion (in terms of age 1+) for stock i
!        real(kind = 8) :: initial_oneplus_i ! Initial one plus for stock i        
!        real(kind = 8) :: scale_pop         ! Scalar to map numbers per recruit to initial numbers at age
!        integer(kind = 4) :: aa             ! Index for age in numbers at age array 
!
!!        print *, "first_yr, yr_max"
!!        print *, first_yr, yr_max
!
!! Time series of population components below (e.g. number of calving females each year)
!!  vector Nplus(first_yr , yr_max)          ! Vector of 1+ population size for all years
!!  vector Ntot(first_yr , yr_max)           ! Total (0+) population size each year, indexed to last year +1 - final pop size is at start of last year (2005?) + 1
!!  vector N_calf(first_yr , yr_max)         ! Vector of calf production for all years
!
!!  3darray NAll(1,4,First_yr,Last_yr+1,0,age_x);
!
!!  vector N_imm_yr_1(First_yr,Last_yr+1);	! Number of immature females by year
!!  vector N_recpt_yr_1(First_yr,Last_yr+1);	! Number of receptive females by year
!!  vector N_calvn_yr_1(First_yr,Last_yr+1);	! Number of calving females by year
!!  vector N_m_yr_1(First_yr,Last_yr+1);		! Numbers of males by year
!        
!!  vector N_dead(First_yr,Last_yr+1);	! Vector of natural dead each year        
!
!!  N_m_yr_1.initialize();
!!  N_calvn_yr_1.initialize();
!!  N_recpt_yr_1.initialize();
!!  N_imm_yr_1.initialize();
!!  
!          
!        scale_pop = 0.5 * k_1plus_i * init_depl_i / initial_oneplus_i ! Assume 50:50 birth rate (so same re-scaling factor used for each sex)  
!        Nage = Nage * scale_pop	! This now in terms of females - so, just assign another vector equal to this one to initialize males
!
!        print *, "Nage: ", Nage ! DEBUGGING    
!!
!!   do aa = 0, age_x
!!   
!!   
!!   end do
!   
!!  for(aa=0;aa<=age_x;aa++) {
!!	  NAll(1,First_yr,aa)=Nage(aa)*(1-prop_mat_a(aa)); 	! immature numbers at age for females
!!	  NAll(2,First_yr,aa)=Nage(aa)*prop_mat_a(aa)*(1-b_init);  ! receptive numbers at age for females
!!	  NAll(3,First_yr,aa)=Nage(aa)*prop_mat_a(aa)*b_init; 	! calving numbers at age for females
!!	  NAll(4,First_yr,aa)=Nage(aa);   			! intial numbers at age for males  	
!!     N_imm_yr_1(First_yr)   += NAll(1,First_yr,aa);		! total numbers at stage
!!     N_recpt_yr_1(First_yr) += NAll(2,First_yr,aa);
!!     N_calvn_yr_1(First_yr) += NAll(3,First_yr,aa);
!!     N_m_yr_1(First_yr)     += NAll(4,First_yr,aa);
!!     Ntot(First_yr)         += NAll(1,First_yr,aa)+NAll(2,First_yr,aa)+NAll(3,First_yr,aa)+NAll(4,First_yr,aa);
!!  }
!
!        do ii = 1, n_stocks
!            depl_i_t(ii , first_yr) = init_depl(ii)
!        end do
!
!        print *, "depl_i_t(): ", depl_i_t(1:2, 0)
!
!!  N_calf(First_yr) = NAll(1,First_yr,0)*2;           // Initial number of calves (for output)
!!  Nplus(First_yr) = Ntot(First_yr)-N_calf(First_yr);
!        
!    end subroutine rescale_NPR
!
!! Function below is R code for assigning life history values to elements of a transition matrix    
!!    real function init_matrix(s_j, s_a, a_mat, age_xx, b_maxx)
!!        implicit none
!!        !S_age(0)
!!    
!!    end function init_matrix
!!    init_matrix = function(s_j, s_a, a_mat, age_xx, b_maxx){
!!  A_matrix = matrix(data = 0, nrow = age_xx, ncol = age_xx) # dimension matrix and fill with zeros
!!  A_matrix[1, (a_mat:age_xx)] = b_maxx # assign birth rates to mature ages in first row of matrix
!!  for(ii in 1:(a_mat-1)) A_matrix[ii+1, ii] = s_j # assign juvenile survival rates
!!  A_matrix[age_xx,a_mat] = s_a # adult survival assumed for maturing animals transitioning into plus-group
!!  A_matrix[age_xx, age_xx] = s_a # adult survival assumed for plus-group
!!  return(A_matrix)
!!}
!!A = init_matrix(s_juv, s_adult, a_m, age_x, b_max) # Call function to initialize matrix A

END MODULE initialize_pop
