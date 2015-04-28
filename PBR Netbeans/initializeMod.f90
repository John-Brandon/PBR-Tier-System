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
    
    contains 
    
    subroutine initialize_age_struc(a_m, npr, S_tmp1, delt_s)
!#######################################################
! Calculate the stable age structure and equilibrium birth rate based on numbers-per-recruit calculations  
! Note: Vectors-at-age are indexed starting at (age) zero, to be consistent with documented methods
! Author: John R. Brandon        
!#######################################################
            implicit none
            
            integer(kind = 4) :: a_m, npr
            integer(kind = 4) :: ii, jj, kk, ll, mm

            real(kind = 8), dimension(0:(a_m+1)) :: Female_age, Male_age    ! Vector with numbers-at-age for females & males
            real(kind = 8), dimension(0:(a_m+1)) :: S_age                ! Vector of survival-at-age
            real(kind = 8), dimension(0:(a_m+1)) :: P_age                ! Vector of rate of maturation-at-age
            real(kind = 8), dimension(0:(a_m+1)) :: prop_mat_a          ! Proportion mature at age (all zeros until a_m)
            real(kind = 8), dimension(0:(a_m+1)) :: NPR_age             ! Numbers-at-age-per-recruit
            real(kind = 8), dimension(0:(a_m+1)) :: Nage_imm_0          ! numbers-at-age that are immature
            real(kind = 8), dimension(0:(a_m+1)) :: Nage_mat_0          ! numbers-at-age that are mature
            real(kind = 8), dimension(0:(a_m+1)) :: prop_NPR            ! Proportions in each age of the NPR vector
            
            real(kind = 8) :: S_tmp1, delt_s       ! delt_s is the difference between juvenile and adult survival
            real(kind = 8) :: temp_1plus, temp_mat, NPR_oneplus
            real(kind = 8) :: b_eq, b_max, b_1, depl_init, theta
            real(kind = 8) :: sum_NPR     ! The sum across ages in the NPR vector
! INITIALIZE ##############################                        
            temp_1plus = 0.0                       ! Temp variable to keep track of numbers aged 1+ years 
            temp_mat = 0.0                         ! Temp variable to keep track of numbers mature
            b_eq = 0.0                             ! Equilibrium birth rate at carrying capacity
            b_1 = 0.0                              ! Birth rate in first year of projection (function of initial depletion level)
            NPR_oneplus = 0.0                      ! Size of the one-plus component scaled to per-recruit at carrying capacity - used to rescale initial recruitment with fishing
            age_x = a_m + 1
! #######################################################
            do ii = 0, age_x               ! assign maturity-at-age vector (0 for immature ages, 1 for mature ages)
                if(ii < a_m) then       
                    prop_mat_a(ii) = 0.0
                else 
                    prop_mat_a(ii) = 1.0
                end if
!                print *, "Hello again from init_age_distribution" ! DEBUGGING
            enddo
                  
            do jj = 0, (a_m - 1)
                S_age(jj) = S_tmp1 - delt_s         ! assign juvenile survival       
            enddo
            
            do kk = a_m, age_x                      ! assign adult survival       
                S_age(kk) = S_tmp1 
            enddo  

            NPR_age(0) = 0.5             ! Numbers of females per recruit, assuming 50:50 sex ratio at birth
            Nage_imm_0(0) = 0.5          ! All calves assumed immature

            do kk = 1, (age_x - 1)       ! NPR calculations from Age 1 to Age x-1
                 NPR_age(kk) = NPR_age(kk-1) * S_age(kk-1)            ! Calculate numbers-at-age per recruit
                 Nage_imm_0(kk) = NPR_age(kk) * (1-prop_mat_a(kk))    ! Numbers immature at age per recruit
                 Nage_mat_0(kk) = NPR_age(kk) * prop_mat_a(kk)        ! Numbers mature at age per recruit
                 temp_1plus = temp_1plus + NPR_age(kk)                ! Keep track of total age one-plus per recruit
                 temp_mat = temp_mat + Nage_mat_0(kk)                 ! "" mature females per recruit
            end  do

            NPR_age(age_x) = NPR_age(age_x - 1) * S_age(age_x - 1) / (1 - S_age(age_x))   ! PLUS GROUP numbers at age
!            Nage_imm_0(age_x) = NPR_age(age_x) * (1 - prop_mat_a(age_x)) 	 ! By definition = 0.0 given assumption that age at transition to plus group equals age at maturity
            Nage_mat_0(age_x) = prop_mat_a(age_x) * NPR_age(age_x)		 ! Mature females

            temp_1plus = temp_1plus + NPR_age(age_x)			 ! Keep track of one_plus component (females)
            temp_mat = temp_mat + Nage_mat_0(age_x)			 ! Keep track of mature female component

            sum_NPR = sum(NPR_age)                                       ! Sum across ages in NPR vector
            print *, "sum_NPR", sum_NPR
            prop_NPR = NPR_age / sum_NPR                                 ! Note this is a vectorized operation in Fortran 90/95
            print *, "prop_NPR", prop_NPR
!            do ll = 0, age_x
!                prop_NPR(ll) = NPR_age(ll) / sum_NPR
!            enddo
!            
            b_eq = 1.0 / (temp_mat - 1)                        ! Equilibrium birth rate on a per recruit basis
            b_1 = b_eq + (b_max - b_eq) * (1 - (depl_init ** theta))   ! Birth rate in first year of projection, given initial depletion

            NPR_oneplus = temp_1plus

            print *, "prop_mat_a = : ", prop_mat_a ! DEBUGGING
            print *, "NPR_age = : ", NPR_age ! DEBUGGING
            print *, "Nage_imm_0 = : ", Nage_imm_0 ! DEBUGGING
            print *, "Nage_mat_0 = : ", Nage_mat_0 ! DEBUGGING
            print *, "temp_1plus = : ", temp_1plus ! DEBUGGING
            print *, "temp_mat = : ", temp_mat ! DEBUGGING
            

            print *, "S_age = : ", S_age ! DEBUGGING
            print *, "S_tmp1 = : ", S_tmp1 ! DEBUGGING  
            print *, "b_eq = : ", b_eq ! DEBUGGING
            print *, "NPR_oneplus = : ", NPR_oneplus ! DEBUGGING
            
        return
    end subroutine initialize_age_struc

    subroutine rescale_NPR() ! Rescale numbers at age per-recruit to initial numbers at age vector
!// Rescale initial numbers at age to 1+ population size given depletion_30
!//####################### 
!  dvariable scale_pop;
!  int aa;
!  
!  scale_pop.initialize();
!  Nplus.initialize();
!  N_calf.initialize();
!  NAll.initialize();
!  Ntot.initialize();
!  N_m_yr_1.initialize();
!  N_calvn_yr_1.initialize();
!  N_recpt_yr_1.initialize();
!  N_imm_yr_1.initialize();
!  
!  scale_pop=0.5*K*depl30/initial_oneplus;
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
