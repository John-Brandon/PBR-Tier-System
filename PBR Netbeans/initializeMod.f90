!#######################################################
! File:   initializeMod.f90
! Author: johnbrandon
!
! Created on April 24, 2015, 12:19 PM
!
! SUBROUTINES:
!  initialize_age_struc()
!#######################################################

MODULE initialize_pop
    
    use Declare_variables_module ! Access global variables  
    
    implicit none
    
    contains 
!!#######################################################    
!    subroutine initialize_age_struc(init_depl, age_x, a_m, &
!        prop_mat_a, S_juv, S_adult, S_age, theta, b_max, &
!        NPR_age, NPR_oneplus) 
!! Calculate the stable age structure and equilibrium birth rate based on numbers-per-recruit calculations  
!! Note: Vectors-at-age are indexed starting at (age) zero, to be consistent with documented methods
!! Author: John R. Brandon               
!!#######################################################
!        real(kind = 8), intent(in) :: init_depl             ! Initial depletion
!        integer(kind = 4), intent(in) :: age_x              ! Plus-group age        
!        integer(kind = 4), intent(in) :: a_m                ! Age at sexual maturity (first estrous)  
!        real(kind = 8), intent(out) :: prop_mat_a(0:age_x)  ! Proportion mature at age
!        real(kind = 8), intent(in) :: S_juv                 ! Juvenile survival rate
!        real(kind = 8), intent(in) :: S_adult               ! Adult survival rate
!        real(kind = 8), intent(inout) :: S_age(0:age_x)     ! Survival rate at age vector
!        real(kind = 8), intent(in) :: theta                 ! Density dependence shape parameter
!        real(kind = 8), intent(in) :: b_max                 ! Maximum birth rate        
!        real(kind = 8), intent(inout) :: NPR_age(0:age_x)   ! Numbers at age per female recruit
!        real(kind = 8), intent(out) :: NPR_oneplus          ! Numbers aged 1+ per female recruit 
!        real(kind = 8) :: Nage_imm_0(0:age_x)
!        real(kind = 8) :: Nage_mat_0(0:age_x)
!        real(kind = 8) :: temp_1plus
!        real(kind = 8) :: temp_mat
!        real(kind = 8) :: sum_NPR
!        real(kind = 8) :: prop_NPR(0:age_x)
!        real(kind = 8) :: b_eq
!        real(kind = 8) :: b_init
!        integer(kind = 4) :: ii, jj, kk                     ! Counters for loops        
!! INITIALIZE ##############################                        
!        print *, "Hello from initialize_age_struc"
!        temp_1plus = 0.0                       ! Temp variable to keep track of numbers aged 1+ years 
!        temp_mat = 0.0                         ! Temp variable to keep track of numbers mature
!        b_eq = 0.0                             ! Equilibrium birth rate at carrying capacity
!        b_init = 0.0                              ! Birth rate in first year of projection (function of initial depletion level)
!        NPR_oneplus = 0.0                      ! Size of the one-plus component scaled to per-recruit at carrying capacity - used to rescale initial recruitment with fishing                      !  
!! ###################################
!!        print *, "S_age: ", S_age
!        S_age(0 : (a_m - 1)) = S_juv              ! Assign juvenile survival rates to ages < a_m
!        S_age(a_m : age_x) = S_adult            ! Assign adult survival                   
!!        print *, "S_age: ", S_age
!        prop_mat_a(0 : (a_m - 1)) = 0.0           ! Assign maturity-at-age vector (0 for immature ages, 1 for mature ages)
!        prop_mat_a(a_m : age_x) = 1.0
!! ###################################              
!        NPR_age(0) = 0.5                    ! Numbers of females per recruit, assuming 50:50 sex ratio at birth
!        Nage_imm_0(0) = NPR_age(0)          ! All calves assumed immature
!
!        do kk = 1, (age_x - 1)                                   ! NPR calculations from Age 1 to Age x-1
!            NPR_age(kk) = NPR_age(kk - 1) * S_age(kk - 1)        ! Calculate numbers-at-age per recruit
!            Nage_imm_0(kk) = NPR_age(kk) * (1 - prop_mat_a(kk))  ! Numbers immature at age per recruit
!            Nage_mat_0(kk) = NPR_age(kk) * prop_mat_a(kk)        ! Numbers mature at age per recruit
!            temp_1plus = temp_1plus + NPR_age(kk)                ! Keep track of total age one-plus per recruit
!            temp_mat = temp_mat + Nage_mat_0(kk)                 ! "" mature females per recruit
!        end do
!
!        NPR_age(age_x) = NPR_age(age_x - 1) * S_age(age_x - 1) / (1 - S_age(age_x)) ! PLUS GROUP numbers at age
!!            Nage_imm_0(age_x) = NPR_age(age_x) * (1 - prop_mat_a(age_x)) 	 ! By definition = 0.0 given assumption that age at transition to plus group equals age at maturity
!        Nage_mat_0(age_x) = prop_mat_a(age_x) * NPR_age(age_x)          ! Mature females
!
!        temp_1plus = temp_1plus + NPR_age(age_x)                        ! Keep track of one_plus component (females)
!        temp_mat = temp_mat + Nage_mat_0(age_x)                         ! Keep track of mature female component
!
!        sum_NPR = sum(NPR_age)                                          ! Sum across ages in NPR vector
!        prop_NPR = NPR_age / sum_NPR                                    ! Note this is a vectorized operation in Fortran 90/95
!
!        print *, "Hello from 'initialize_age_struc'"
!        print *, "init_depl: ", init_depl
!
!        b_eq = 1.0 / temp_mat                                           ! Equilibrium birth rate on a per recruit basis
!        print *, "b_eq: ", b_eq        
!!        print *, "b_max: ", b_max                
!!        print *, "theta: ", theta                        
!        
!        b_init = b_eq + (b_max - b_eq) * (1 - (init_depl ** theta))        ! Birth rate in first year of projection, given initial depletion
!        print *, "b_init: ", b_init
!        
!        NPR_oneplus = temp_1plus
!        print *, "NPR_oneplus: ", NPR_oneplus
!            
!        return
!    end subroutine initialize_age_struc
    
!! ####################### 
!    real(kind = 8) function calc_NPR(f_init, selectivity, b_max, theta, a_r, a_m, age_x, &
!        S_adult, S_juv, S_age, prop_mat_age)
!! This function calculates an initial age vector based on a numbers per recruit (i.e. female calf) approach.     
!! It assumes that the sex ratio at birth is 50:50, i.e. there are 0.50 age zero females per calf
!        real(kind = 8), intent(in) :: f_init
!        real(kind = 8), intent(inout) :: selectivity
!        real(kind = 8), intent(in) :: b_max
!        real(kind = 8), intent(in) :: theta
!        integer(kind = 4), intent(inout) :: a_r
!        integer(kind = 4), intent(inout) :: a_m
!        integer(kind = 4), intent(inout) :: age_x
!        real(kind = 8), intent(inout) :: S_adult
!        real(kind = 8), intent(inout) :: S_juv        
!        real(kind = 8), intent(out) :: S_age(0:age_x)
!        real(kind = 8), intent(out) :: prop_mat_age(0:age_x)        
!
!! Assign input values to vectors with selectivity, proportion mature and survival at age        
!!        Call assign_par_vectors(a_r, a_m, age_x, S_adult, S_juv, S_age, selectivity, prop_mat_age)
!        
!!        Call calc_age_vector(age_x, f_rate, selectivity, S_age, prop_mat_age, rec, N_age)
!        
!        calc_NPR = 0.
!        
!        return    
!    end function calc_NPR 
    
! ####################### ++ ####################### ++ #######################
    subroutine calc_NPR_age(f_rate, N_recruits, N_age, sum_1plus, sum_mature)
! TODO : Add description of subroutine here        
! ####################### ++ ####################### ++ #######################        
!        integer(kind = 4), intent(in) :: age_x              ! Plus group age (from input.par file)
        real(kind = 8), intent(in) :: f_rate                ! Human caused mortality rate
!        real(kind = 8), intent(in) :: selectivity(0:age_x)  ! Selectivity at age vector        
!        real(kind = 8), intent(in) :: S_age(0:age_x)        ! Survival rate at age vector
!        real(kind = 8), intent(in) :: prop_mat_age(0:age_x) ! Proportion mature at age vector
        real(kind = 8), intent(in) :: N_recruits            ! Numbers at age zero (female calves)
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
        
!        print *, "From calc_age_vector: "
!        print *,"b_eq: ", b_eq
!        print *,"b_init: ", b_init
!        print *,"b_max: ", b_max
!        print *,"init_depl_i: ", init_depl_i
!        print *,"theta: ", theta

!        print *, "S_age: "
!        print *, S_age
!        print *, "selectivity: "
!        print *, selectivity        
!        print *, "f_rate: "
!        print *, f_rate        
!        print *, "N_age: "
!        print *, N_age

        return
    end subroutine calc_NPR_age   
    
! #######################  
    real(kind = 8) function Initial_F(f_init) result(objf_f_init)
! #######################    
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
! ####################### 
        implicit none
        real(kind = 8), intent(inout) :: f_init
        integer(kind = 4) :: ii, jj				! Index counters for loops
        real(kind = 8) :: temp_mature
        real(kind = 8) :: temp_1plus
        real(kind = 8) :: foo
        real(kind = 8) :: rec_init ! initial recruitment conditioned on Finit (in terms of females)
        real(kind = 8) :: pred_rec
        real(kind = 8) :: initial_oneplus
        real(kind = 8) :: objf_f_init	! Objective function to be minimized in finding root at f_init

! INITIALIZE variables in Initial_F() ##############################        
        temp_mature = 0.    
        temp_1plus = 0.     
        rec_init = 0.       
        pred_rec = 0.       
        initial_oneplus = 0.
        objf_f_init = -99.d0
              
        Call calc_NPR_age(f_rate = f_init, N_recruits = b_sex_ratio, &
            N_age = NPR_age_tmp, sum_1plus = temp_1plus, &
            sum_mature = temp_mature)        
! ##############################                      
!  NPR_oneplus is global variable calculated at f_init = 0 
        rec_init = init_depl_i * NPR_oneplus / temp_1plus 	 ! Rescale to get the desired 1+ depletion conditioned on Finit 
        rec_init = rec_init / 2                                  ! Divide by two to get female recruits conditioned on Finit
!        print *, "Hello from Initial_F(): "
!        print *, "f_init: ", f_init
!        print *, "NPR_oneplus: ", NPR_oneplus
!        print *, "temp_1plus: ", temp_1plus  		  		
!        print *, "temp_1plus: ", temp_1plus        
!        print *, "temp_mature: ", temp_mature  		
!        print *, "prop_mat_age: "
!        print *, prop_mat_age
!        print *, "selectivity: "
!        print *, selectivity
        
! Calculate actual numbers per recruit conditioned on Finit
        NPR_age_tmp(0) = rec_init
        Call calc_NPR_age(f_rate = f_init, N_recruits = rec_init, &
            N_age = NPR_age_tmp, sum_1plus = initial_oneplus, &
            sum_mature = foo) 

        pred_rec = b_init * temp_mature 
        pred_rec = pred_rec - (init_depl_i * NPR_oneplus / initial_oneplus)
        pred_rec = 0.5 * pred_rec

        objf_f_init = pred_rec * pred_rec 
        
!        print *, "NPR_age_tmp: "
!        print *, NPR_age_tmp
!        print *, "initial_oneplus: ", initial_oneplus  
!        print *, "b_init: ", b_init  
!        print *, "pred_rec: ", pred_rec          
        
!        Initial_F = objf_f_init
        
        return
    end function Initial_F    
! #######################  
    
! #######################      
    real(kind = 8) function initial_F_oldcode(f_init)
! #######################    
! Given numbers per-recruit at equilibrium - apply initial human caused mortality rate (F_init)
! Global Vars    
!        real(kind = 8) :: NPR_age                  ! Numbers per female recruit, see initialize_age_struc()      
!        real(kind = 8) :: Nage                     ! Numbers at age conditioned on F_init
!        real(kind = 8) :: selectivity(age_x)       ! Selectivity at age
!        integer(kind = 4) :: a_r                   ! Age at recruitment (vulnerability) to human caused mortality
!        integer(kind = 4) :: age_x                 ! Plus group age (user input)
!        real(kind = 8) :: temp_mature          
!        real(kind = 8) :: temp_1plus
!        real(kind = 8) :: rec_init ! initial recruitment conditioned on Finit (in terms of females)
!        real(kind = 8) :: pred_rec
!        real(kind = 8) :: b_init
!        real(kind = 8) :: initial_oneplus
! ####################### 
        implicit none
        real(kind = 8), intent(inout) :: f_init  
        real(kind = 8) :: rec_init      ! Initial recruitment conditioned on Finit (in terms of females)
        real(kind = 8) :: pred_rec      ! Predicted recruitment given some value for F_init (used to solve for equilibrium F_init)
        real(kind = 8) :: objf_f_init     ! Objective function value for solving for initial human caused mortality rate
        real(kind = 8) :: temp_1plus    !
        real(kind = 8) :: temp_mat      !
        real(kind = 8) :: temp_mature   !
        real(kind = 8) :: initial_oneplus        
        real(kind = 8) :: prop_mat_a(0:age_x)
! INITIALIZE variables in Initial_F() ##############################        
        temp_mature = 0.    ! Reset to zero after they were calculated in initialize_age_struc() subroutine 
        temp_1plus = 0.     ! Reset to zero after they were calculated in initialize_age_struc() subroutine 
        rec_init = 0.       ! Reset to zero after they were calculated in initialize_age_struc() subroutine
        pred_rec = 0.       ! Reset to zero after they were calculated in initialize_age_struc() subroutine
        initial_oneplus = 0.! Reset to zero after they were calculated in initialize_age_struc() subroutine
        
        do ii = 0, age_x    ! Implement assumption that selectivity is knife-edge at age 1+
            if(ii < a_r) then
                selectivity(ii) = 0.0
            else 
                selectivity(ii) = 1.0
            end if
            if(ii < a_m) then
                prop_mat_a(ii) = 0.0
            else 
                prop_mat_a(ii) = 1.0
            end if
        end do
              
        Nage = NPR_age      ! Initialize numbers at age given number per recruit -- using a dummy vector for this routine
        
! ##############################                
        do ii = 1, (age_x -1)
            Nage(ii) = Nage(ii - 1) * S_age(ii - 1) ! Natural mortality - Sage vector initialized in initialize_age_struc() subroutine
            Nage(ii) = Nage(ii) * (1 - selectivity(ii-1) * f_init) ! Apply fishing mortality at start of time series
            temp_mature = temp_mature + prop_mat_a(ii) * Nage(ii) ! Track numbers of mature per recruit
            temp_1plus = temp_1plus + Nage(ii)                    ! Track ages 1+ years per recruit
        end do
!
        Nage(age_x) = Nage(age_x - 1) * S_age(age_x - 1)           ! plus group per recruit
        Nage(age_x) = Nage(age_x) * (1 - selectivity(age_x - 1) * f_init)
        Nage(age_x) = Nage(age_x) / (1 - S_age(age_x) * (1 - selectivity(age_x) * f_init))
        temp_mature = temp_mature + prop_mat_a(age_x) * Nage(age_x)
        temp_1plus = temp_1plus + Nage(age_x)
!  
        rec_init = init_depl_i * NPR_oneplus / temp_1plus 	 ! Rescale to get the desired 1+ depletion conditioned on Finit 
        rec_init = rec_init / 2                                  ! Divide by two to get female recruits conditioned on Finit

!        print *, "Hello from Initial_F(): "
!        print *, "f_init: ", f_init
!        print *, "NPR_oneplus: ", NPR_oneplus
!        print *, "temp_1plus: ", temp_1plus  		  		     
!        print *, "temp_mature: ", temp_mature
!        print *, "prop_mat_a: ", prop_mat_a
!        print *, "selectivity: ", selectivity
! Calculate actual numbers per recruit conditioned on Finit
        Nage(0) = rec_init
        do jj = 1, (age_x - 1)
            Nage(jj) = Nage(jj - 1) * S_age(jj - 1)                      ! Natural mortality
            Nage(jj) = Nage(jj) * (1 - selectivity(jj - 1) * f_init) ! Human caused mortality at start of time series
            initial_oneplus = initial_oneplus + Nage(jj)
        end do
!   
        Nage(age_x) = Nage(age_x - 1) * S_age(age_x - 1)                 ! Plus group per recruit
        Nage(age_x) = Nage(age_x) * (1 - selectivity(age_x - 1) * f_init)
        Nage(age_x) = Nage(age_x) / (1 - S_age(age_x) * (1 - selectivity(age_x) * f_init))
        initial_oneplus = initial_oneplus + Nage(age_x)

        pred_rec = b_init * temp_mature 
        pred_rec = pred_rec - (init_depl_i * NPR_oneplus / initial_oneplus)
        pred_rec = 0.5 * pred_rec

        objf_f_init = pred_rec * pred_rec 

!        print *, "Nage: "
!        print *, Nage
!        print *, "initial_oneplus: ", initial_oneplus  
!        print *, "b_init: ", b_init  
!        print *, "pred_rec: ", pred_rec  
        
        initial_F_oldcode = objf_f_init
        
        return
        
    end function initial_F_oldcode
!    
! ####################### ++ ####################### ++ #######################
    subroutine assign_par_vectors(a_r, a_m, age_x, S_adult, S_juv, S_age, selectivity_age, prop_mat_age)
! Assign values to survival, proportion mature and selectivity at age vectors
! ####################### ++ ####################### ++ #######################        
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
!!#######################     
!    subroutine rescale_NPR(k_1plus_i, init_depl_i, initial_oneplus_i) 
!!####################### 
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
