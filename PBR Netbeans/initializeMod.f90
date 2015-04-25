!     
! File:   initializeMod.f90
! Author: johnbrandon
!
! Created on April 24, 2015, 12:19 PM
!

MODULE initialize_pop
    
    contains 
    
    subroutine initialize_age_struc(a_m, npr, S_tmp1, delt_s)
!#######################################################
! Calculate the stable age structure and equilibrium birth rate based on numbers-per-recruit calculations  
! Note: Vectors-at-age are indexed starting at (age) zero, to be consistent with documented methods
! Author: John R. Brandon        
!#######################################################
            implicit none
            
            integer(kind = 4) :: a_m, npr
            integer(kind = 4) :: ii, jj, kk

            real(kind = 8), dimension(0:(a_m+1)) :: Female_age, Male_age    ! Vector with numbers-at-age for females & males
            real(kind = 8), dimension(0:(a_m+1)) :: S_age                ! Vector of survival-at-age
            real(kind = 8), dimension(0:(a_m+1)) :: P_age                ! Vector of rate of maturation-at-age
            real(kind = 8), dimension(0:(a_m+1)) :: prop_mat_a          ! Proportion mature at age (all zeros until a_m)
            real(kind = 8), dimension(0:(a_m+1)) :: NPR_age             ! Numbers-at-age-per-recruit
            real(kind = 8), dimension(0:(a_m+1)) :: Nage_imm_0          ! numbers-at-age that are immature
            real(kind = 8), dimension(0:(a_m+1)) :: Nage_mat_0          ! numbers-at-age that are mature
            
            real(kind = 8) :: S_tmp1, delt_s       ! delt_s is the difference between juvenile and adult survival
            real(kind = 8) :: temp_1plus, temp_mat, NPR_oneplus
            real(kind = 8) :: b_eq, b_max, b_1, depl_init, theta
! INITIALIZE ##############################                        
            temp_1plus = 0.0                       ! Temp variable to keep track of numbers aged 1+ years 
            temp_mat = 0.0                         ! Temp variable to keep track of numbers mature
            b_eq = 0.0                             ! Equilibrium birth rate at carrying capacity
            b_1 = 0.0                              ! Birth rate in first year of projection (function of initial depletion level)
            NPR_oneplus = 0.0                      ! Size of the one-plus component scaled to per-recruit at carrying capacity - used to rescale initial recruitment with fishing

! #######################################################
            do 10 ii = 1, a_m               ! assign maturity-at-age vector (0 for immature ages, 1 for mature ages)
                if(ii <= a_m) then
                    prop_mat_a(ii) = 0.0
                else 
                    prop_mat_a(ii) = 1.0
                end if
                print *, "Hello again from init_age_distribution" ! DEBUGGING
10          continue
! #######################################################
! This rate of maturation was used in ADMB code to have a differentiable function (logistic maturity ogive) 
!  which made it possible to estimate age at 50% maturity as a parameter
!  because we're not estimating anything for the PBR code, this approach is not necessary
! #######################################################
            do 11 jj = 1, a_m               ! calculate rate of maturation, given proportion mature at age
              if(1-prop_mat_a(jj)<0.001) then
               P_age(jj-1) = 1.0
              else
               P_age(jj-1) = (prop_mat_a(jj)-prop_mat_a(jj-1))
               P_age(jj-1) = P_age(jj-1) / (1-prop_mat_a(jj-1))
              end if
11          continue              
         
            S_age(0) = S_tmp1 - delt_s      ! assign calf survival

            do kk = 1,(a_m + 1)             ! fill survival-at-age vector
                S_age(kk) = S_tmp1 
                print *, "age = : ", kk           
            enddo   ! syntax for looping, available in f90 

            NPR_age(0) = 0.5             ! Numbers of females per recruit, assuming 50:50 sex ratio at birth
            Nage_imm_0(0) = 0.5          ! All calves assumed immature
            
            do kk = 1, (a_m - 1)
                 NPR_age(kk) = NPR_age(kk-1) * S_age(kk-1)            ! Calculate numbers-at-age per recruit
                 Nage_imm_0(kk) = NPR_age(kk) * (1-prop_mat_a(kk))    ! Initialize Numbers-at-Stage, based on maturity ogive
                 Nage_mat_0(kk) = prop_mat_a(kk) * NPR_age(kk)
                 temp_1plus = temp_1plus + NPR_age(kk)                ! Keep track of one-plus females per recruit
                 temp_mat = temp_mat + Nage_mat_0(kk)                 ! "" mature females per recruit
            end do

            NPR_age(a_m) = NPR_age(a_m - 1) * S_age(a_m - 1) / (1 - S_age(a_m - 1))   ! PLUS GROUP numbers at age
            Nage_imm_0(a_m) = NPR_age(a_m) * (1 - prop_mat_a(a_m)) 	 ! By definition = 0.0 given assumption that age at transition to plus group equals age at maturity
            Nage_mat_0(a_m) = prop_mat_a(a_m) * NPR_age(a_m)		 ! Mature females

            temp_1plus = temp_1plus + NPR_age(a_m)			 ! Keep track of one_plus component (females)
            temp_mat = temp_mat + Nage_mat_0(a_m)			 ! Keep track of mature female component

            b_eq = 1.0 / (S_tmp1 * (temp_mat - 1))                        ! Equilibrium birth rate on a per recruit basis
            b_1 = b_eq 
            b_1 = b_1 + (b_max - b_eq) * (1 - (depl_init ** theta))   ! Birth rate in first year of projection, given initial depletion
            NPR_oneplus = temp_1plus

            print *, "S_age = : ", S_age ! DEBUGGING
            print *, "S_age = : ", S_age ! DEBUGGING
            
        return
    end subroutine initialize_age_struc

END MODULE initialize_pop
