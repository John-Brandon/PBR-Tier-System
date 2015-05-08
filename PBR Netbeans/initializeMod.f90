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
    
    use Declare_variables_module
    
    implicit none
       
    contains 
    
    subroutine initialize_age_struc() !(a_m, npr, S_tmp1, delt_s)
!#######################################################
! Calculate the stable age structure and equilibrium birth rate based on numbers-per-recruit calculations  
! Note: Vectors-at-age are indexed starting at (age) zero, to be consistent with documented methods
! Author: John R. Brandon        
!#######################################################
        implicit none ! Not necessary if this routine is located inside a module (e.g. initialize_pop) that has "implicit none" declared across module processes (routines and functions).
            
!        integer(kind = 4) :: a_m, npr
!        integer(kind = 4) :: ii, jj, kk, ll, mm
!
!        real(kind = 8), dimension(0:(a_m+1)) :: Female_age, Male_age    ! Vector with numbers-at-age for females & males
!        real(kind = 8), dimension(0:(a_m+1)) :: S_age                ! Vector of survival-at-age
!        real(kind = 8), dimension(0:(a_m+1)) :: P_age                ! Vector of rate of maturation-at-age
!        real(kind = 8), dimension(0:(a_m+1)) :: prop_mat_a          ! Proportion mature at age (all zeros until a_m)
!        real(kind = 8), dimension(0:(a_m+1)) :: NPR_age             ! Numbers-at-age-per-recruit
!        real(kind = 8), dimension(0:(a_m+1)) :: Nage_imm_0          ! numbers-at-age that are immature
!        real(kind = 8), dimension(0:(a_m+1)) :: Nage_mat_0          ! numbers-at-age that are mature
!        real(kind = 8), dimension(0:(a_m+1)) :: prop_NPR            ! Proportions in each age of the NPR vector
!
!        real(kind = 8) :: S_tmp1, delt_s       ! delt_s is the difference between juvenile and adult survival
!        real(kind = 8) :: temp_1plus, temp_mat, NPR_oneplus
!        real(kind = 8) :: b_eq, b_max, b_1, depl_init, theta
!        real(kind = 8) :: sum_NPR     ! The sum across ages in the NPR vector

! INITIALIZE ##############################                        
        temp_1plus = 0.0                       ! Temp variable to keep track of numbers aged 1+ years 
        temp_mat = 0.0                         ! Temp variable to keep track of numbers mature
        b_eq = 0.0                             ! Equilibrium birth rate at carrying capacity
        b_1 = 0.0                              ! Birth rate in first year of projection (function of initial depletion level)
        NPR_oneplus = 0.0                      ! Size of the one-plus component scaled to per-recruit at carrying capacity - used to rescale initial recruitment with fishing
        age_x = a_m + 1
! #########################################
        do ii = 0, age_x               ! assign maturity-at-age vector (0 for immature ages, 1 for mature ages)
            if(ii < a_m) then       
                prop_mat_a(ii) = 0.0
            else 
                prop_mat_a(ii) = 1.0
            end if
        enddo
        print *, "prop_mat_a: ", prop_mat_a
! #########################################        
        print *, "Assigning survival vector: ", S_juv
        print *, "S_juv: ", S_juv
        print *, "S_adult: ", S_adult   
        print *, "delt_s: ", delt_s        

        do jj = 0, (a_m - 1)
            S_age(jj) = S_juv         ! assign juvenile survival       
        enddo

        do kk = a_m, age_x            ! assign adult survival       
            S_age(kk) = S_adult 
        enddo  

        NPR_age(0) = 0.5             ! Numbers of females per recruit, assuming 50:50 sex ratio at birth
        Nage_imm_0(0) = NPR_age(0)          ! All calves assumed immature

        do kk = 1, (age_x - 1)       ! NPR calculations from Age 1 to Age x-1
            NPR_age(kk) = NPR_age(kk - 1) * S_age(kk - 1)            ! Calculate numbers-at-age per recruit
            Nage_imm_0(kk) = NPR_age(kk) * (1 - prop_mat_a(kk))    ! Numbers immature at age per recruit
            Nage_mat_0(kk) = NPR_age(kk) * prop_mat_a(kk)        ! Numbers mature at age per recruit
            temp_1plus = temp_1plus + NPR_age(kk)                ! Keep track of total age one-plus per recruit
            temp_mat = temp_mat + Nage_mat_0(kk)                 ! "" mature females per recruit
        end  do

        NPR_age(age_x) = NPR_age(age_x - 1) * S_age(age_x - 1) / (1 - S_age(age_x)) ! PLUS GROUP numbers at age
!            Nage_imm_0(age_x) = NPR_age(age_x) * (1 - prop_mat_a(age_x)) 	 ! By definition = 0.0 given assumption that age at transition to plus group equals age at maturity
        Nage_mat_0(age_x) = prop_mat_a(age_x) * NPR_age(age_x)          ! Mature females

        temp_1plus = temp_1plus + NPR_age(age_x)                        ! Keep track of one_plus component (females)
        temp_mat = temp_mat + Nage_mat_0(age_x)                         ! Keep track of mature female component

        sum_NPR = sum(NPR_age)                                          ! Sum across ages in NPR vector
        print *, "sum_NPR", sum_NPR
        prop_NPR = NPR_age / sum_NPR                                    ! Note this is a vectorized operation in Fortran 90/95
        print *, "prop_NPR", prop_NPR           ! DEBUGGING

        b_eq = 1.0 / temp_mat                                     ! Equilibrium birth rate on a per recruit basis
        b_1 = b_eq + (b_max - b_eq) * (1 - (init_depl ** theta))        ! Birth rate in first year of projection, given initial depletion

        NPR_oneplus = temp_1plus

        print *, "prop_mat_a = : ", prop_mat_a  ! DEBUGGING
        print *, "NPR_age = : ", NPR_age        ! DEBUGGING
        print *, "Nage_imm_0 = : ", Nage_imm_0  ! DEBUGGING
        print *, "Nage_mat_0 = : ", Nage_mat_0  ! DEBUGGING
        print *, "temp_1plus = : ", temp_1plus  ! DEBUGGING
        print *, "temp_mat = : ", temp_mat      ! DEBUGGING


        print *, "S_age = : ", S_age            ! DEBUGGING
        print *, "S_tmp1 = : ", S_tmp1          ! DEBUGGING  
        print *, "b_eq = : ", b_eq              ! DEBUGGING
        print *, "b_max = : ", b_max              ! DEBUGGING        
        print *, "init_depl = : ", init_depl      ! DEBUGGING                
        print *, "b_1 = : ", b_1              ! DEBUGGING        
        print *, "NPR_oneplus = : ", NPR_oneplus ! DEBUGGING
        print *, "prop_NPR = : ", prop_NPR      ! DEBUGGING
            
        return
    end subroutine initialize_age_struc
!//####################### 

!//#######################  
    real(kind = 8) function Initial_F(f_init)
!// given numbers per-recruit at equilibrium - apply initial human caused mortality rate (F_init)
!//####################### 
        ! use
        implicit none

!        real(kind = 8) :: temp_mature, temp_1plus
!        real(kind = 8) :: rec_init ! initial recruitment conditioned on Finit (in terms of females)
!        real(kind = 8) :: pred_rec, initial_oneplus
         real(kind = 8) :: f_init
! INITIALIZE variables in Initial_F() ##############################        
        temp_mature = 0.    ! Reset to zero after they were calculated in initialize_age_struc() subroutine 
        temp_1plus = 0.     ! Reset to zero after they were calculated in initialize_age_struc() subroutine 
        rec_init = 0.
        pred_rec = 0.
        initial_oneplus = 0.
        
        do ii = 0, age_x    ! Implement assumption that selectivity is knife-edge at age 1+
            if(ii < 1) then
                selectivity(ii) = 0.0
            else 
                selectivity(ii) = 1.0
            end if
        end do
        
        !Initial_F = 0.20    ! Starting guess 
        Nage = NPR_age      ! Initialize numbers at age given number per recruit -- using a dummy vector for this routine
!        print *, "Hello from Initial_F"
!        print *, "Nage: ", Nage  
!        print *, "Input for f_init: ", f_init          
        
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
        rec_init = init_depl * NPR_oneplus / temp_1plus 	 ! now, rescale to get the desired 1+ depletion conditioned on Finit 
        rec_init = rec_init / 2                                  ! divide by two to get female recruits conditioned on Finit
!  
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

        print *, ""
        print *, "Hello from initial_F()"
        print *, "After applying f_init = ", f_init
!        print *, "Nage: ", Nage        
        
!         pred_rec = (0.50 * b_1 * temp_mature) 
!         print *, "pred_rec (temp1): ", pred_rec
!         print *, "b_1: ", b_1
!         print *, "temp_mature: ", temp_mature
!         pred_rec = 0.50 * init_depl * NPR_oneplus / initial_oneplus
!         print *, "pred_rec (temp2): ", pred_rec
!         print *, "init_depl: ", init_depl
!         print *, "NPR_oneplus: ", NPR_oneplus
!         print *, "initial_oneplus: ", initial_oneplus
         pred_rec = b_1 * temp_mature 
         pred_rec = pred_rec - (init_depl * NPR_oneplus / initial_oneplus)
         pred_rec = 0.5 * pred_rec
!         print *, "pred_rec (temp3): ", pred_rec         
         !        pred_rec = (b_1 * temp_mature) - (init_depl * NPR_oneplus / initial_oneplus)        
!        pred_rec = 0.50 * pred_rec
        objf_newt = pred_rec * pred_rec

        ! DEBUGGING (Subtracting one from objective function in hopes of helping ZBRENT root finding)
        !objf_newt = objf_newt ! - 0.00000001

        print *, "objf_newt: ", objf_newt
        print *, ""
        Initial_F = objf_newt
        return
    end function Initial_F
   
!#######################  
    real(kind = 8) function newtons_root()
!#######################
! Using the Newton-Raphson method, find the root of a function known to lie in the interval [x1,x2].
!  The root (rtnewt) will be refined until its accuracty is known within +/- (xacc). (funcd) is a user
!  supplied subroutine that returns both the function value and the first derivative of the function at the point (x).'
!  - description from Numerical Recipes book
!  In this case going to be solving for human caused mortality rate which would leave the population at initial depletion level
!#######################
        implicit none
! Declare local variables        
        integer :: jcount  ! counter specific to this function
        integer :: jmax    ! maximum number of iterations
        real(kind = 8) :: delta_Finit
        real(kind = 8) :: rtnewt
        real(kind = 8) :: f1_newt
        real(kind = 8) :: df_newt
        real(kind = 8) :: dx_newt
        real(kind = 8) :: f_init
 !       real(kind = 8), external :: Initial_F
!
        f_init = 0.
        f1_newt = 0.
        objf_newt = 0.
!
        jmax = 20
        delta_Finit = 0.000001
        rtnewt = 0.50             ! Initial guess for human caused mortality rate in year zero

        do jcount = 1, jmax
              print *, "jcount: ", jcount
              f_init = rtnewt
              objf_newt = Initial_F(f_init)              ! Subroutine call
              f1_newt = objf_newt	
              f_init = rtnewt + delta_Finit
              objf_newt = Initial_F(f_init)              ! Subroutine call
              df_newt = (objf_newt - f1_newt) / delta_Finit
              dx_newt = f1_newt / df_newt
              rtnewt = rtnewt - dx_newt
              print *, "rtnewt: ", rtnewt              
        end do

        print *, "Hello from newtons_root"
        print *, "Max iterations: ", jmax
        print *, "Iterations completed: ", jcount - 1        
        print *, "newtons_root: ", newtons_root
        newtons_root = rtnewt
        return
        
    end function newtons_root
!#######################     
    subroutine rescale_NPR() 
!####################### 
! Rescale numbers at age per-recruit to initial numbers at age vector, given initial depletion 
        implicit none
        
        real(kind = 8) :: scale_pop ! scalar to map numbers per recruit to initial numbers at age
        integer :: aa ! index for age in numbers at age array 
        print *, "first_yr, yr_max"
        print *, first_yr, yr_max
! time series of population components below (e.g. number of calving females each year)
!  vector Nplus(First_yr,Last_yr+1);	// Vector of 1+ population size for all years
!  vector Ntot(First_yr,Last_yr+1);	// Total (0+) population size each year, indexed to last year +1 - final pop size is at start of last year (2005?) + 1
!  vector N_calf(First_yr,Last_yr+1);	// Vector of calf production for all years

!  3darray NAll(1,4,First_yr,Last_yr+1,0,age_x);

!  vector N_imm_yr_1(First_yr,Last_yr+1);	// Number of immature females by year
!  vector N_recpt_yr_1(First_yr,Last_yr+1);	// Number of receptive females by year
!  vector N_calvn_yr_1(First_yr,Last_yr+1);	// Number of calving females by year
!  vector N_m_yr_1(First_yr,Last_yr+1);		// Numbers of males by year
        
!  vector N_dead(First_yr,Last_yr+1);	// Vector of natural dead each year        


!  N_m_yr_1.initialize();
!  N_calvn_yr_1.initialize();
!  N_recpt_yr_1.initialize();
!  N_imm_yr_1.initialize();
!  
!  scale_pop = 0.5 * k_1plus * init_depl / initial_oneplus
!  Nage = Nage*scale_pop;	// this now in terms of females - so, just assign another vector equal to this one to initialize males
!
!  for(aa=0;aa<=age_x;aa++) {
!	  NAll(1,First_yr,aa)=Nage(aa)*(1-prop_mat_a(aa)); 	// immature numbers at age for females
!	  NAll(2,First_yr,aa)=Nage(aa)*prop_mat_a(aa)*(1-b_1); // receptive numbers at age for females
!	  NAll(3,First_yr,aa)=Nage(aa)*prop_mat_a(aa)*b_1; 	// calving numbers at age for females
!	  NAll(4,First_yr,aa)=Nage(aa);   							// intial numbers at age for males  	
!     N_imm_yr_1(First_yr)   += NAll(1,First_yr,aa);		// total numbers at stage
!     N_recpt_yr_1(First_yr) += NAll(2,First_yr,aa);
!     N_calvn_yr_1(First_yr) += NAll(3,First_yr,aa);
!     N_m_yr_1(First_yr)     += NAll(4,First_yr,aa);
!     Ntot(First_yr)         += NAll(1,First_yr,aa)+NAll(2,First_yr,aa)+NAll(3,First_yr,aa)+NAll(4,First_yr,aa);
!  }
!  depl(First_yr)=depl30;
!  N30=K*depl30;		//treating depletion in 1930 as esitmated parameter, so multiply by K to get 1+ population size in 1930
!  N_calf(First_yr) = NAll(1,First_yr,0)*2;           // Initial number of calves (for output)
!  Nplus(First_yr) = Ntot(First_yr)-N_calf(First_yr);
        
    end subroutine  

! Function below is R code for assigning life history values to elements of a transition matrix    
!    real function init_matrix(s_j, s_a, a_mat, age_xx, b_maxx)
!        implicit none
!        !S_age(0)
!    
!    end function init_matrix
!    init_matrix = function(s_j, s_a, a_mat, age_xx, b_maxx){
!  A_matrix = matrix(data = 0, nrow = age_xx, ncol = age_xx) # dimension matrix and fill with zeros
!  A_matrix[1, (a_mat:age_xx)] = b_maxx # assign birth rates to mature ages in first row of matrix
!  for(ii in 1:(a_mat-1)) A_matrix[ii+1, ii] = s_j # assign juvenile survival rates
!  A_matrix[age_xx,a_mat] = s_a # adult survival assumed for maturing animals transitioning into plus-group
!  A_matrix[age_xx, age_xx] = s_a # adult survival assumed for plus-group
!  return(A_matrix)
!}
!A = init_matrix(s_juv, s_adult, a_m, age_x, b_max) # Call function to initialize matrix A

END MODULE initialize_pop
