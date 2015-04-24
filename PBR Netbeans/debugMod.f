!     
! File:   debugMod.f90
! Author: johnbrandon
!
! Created on April 23, 2015, 6:38 PM
!


      module debug

         contains 
    
         real function debugz(a_m)
            real(kind = 4) :: a_m
            debugz = a_m / 4
            print *, "Hello from debugz"
            return
         end function debugz

         subroutine init_age_distribution(a_m, npr)
!#######################################################
! Calculate the stable age structure and equilibrium birth rate based on numbers-per-recruit calculations        
!#######################################################
            integer(kind = 4) :: a_m, npr
            integer(kind = 4) :: ii, jj, kk
! INITIALIZE ##############################
            real(kind = 4) Male_age(a_m+1)        ! Vector with numbers-at-age for males
            real(kind = 4) Female_age(a_m+1)          ! "" females

            real(kind = 4) Sage(a_m+1)                ! Vector of survival-at-age
            real(kind = 4) Page(a_m+1)                ! Vector of rate of maturation-at-age
            real(kind = 4) prop_mat_a(a_m+1)           ! Proportion mature at age (all zeros until a_m)
            real(kind = 4) NPR_age(a_m+1)             ! Numbers-at-age-per-recruit
            real(kind = 4) Nage_imm_0(a_m+1)          !
            real(kind = 4) Nage_mat_0(a_m+1)          ! numbers-at-age that are mature

            temp_1plus = 0.0;
            temp_mat = 0.0;
            b_eq = 0.0;                             ! equilibrium birth rate at carrying capacity
            NPR_oneplus = 0.0;                      ! size of the one-plus component scaled to per-recruit at carrying capacity - used to rescale initial recruitment with fishing


! //////////////////////////////////////////////////
            do 10 ii = 1, a_m
                if(ii <= a_m) then
                    prop_mat_a(ii) = 0.0
                else 
                    prop_mat_a(ii) = 1.0
                end if
                print *, "Hello from init_age_distribution"
10          continue

!       do 11 jj = 1, a_m               ! calculate rate of maturation, given proportion mature at age
!         if(1-prop_mat_a(jj)<0.001) then
!          Page(jj-1) = 1.0
!         else
!          Page(jj-1) = (prop_mat_a(jj)-prop_mat_a(jj-1))/(1-prop_mat_a(jj-1))
!         end if
!11     continue              

        return
        end subroutine init_age_distribution
         
      end module debug
!
!      
!
!          public void Do_NPR(double S_tmp1){
!
!
!        for (jj=1;jj<=a_m;jj++)               // calculate rate of maturation, given proportion mature at age
!         {
!          if(1-prop_mat_a[jj]<0.001)
!           Page[jj-1] = 1.0;
!          else
!           Page[jj-1] = (prop_mat_a[jj]-prop_mat_a[jj-1])/(1-prop_mat_a[jj-1]);
!         }
!
!        Sage[0] = S_tmp1 - delt_s;             // assign calf survival
!
!        for(ii = 1; ii <= a_m; ii++)
!            Sage[ii] = S_tmp1;                 // fill survival-at-age vector
!
!        NPR_age[0]=0.5;                        // Numbers of females per recruit, assuming 50:50 sex ratio at birth
!        Nage_imm_0[0]=0.5;                     // All calves assumed immature
!
!       for(kk=1;kk<=a_m-1;kk++){
!         NPR_age[kk] = NPR_age[kk-1]*Sage[kk-1];             //calculate numbers-at-age per recruit
!         Nage_imm_0[kk] = NPR_age[kk]*(1-prop_mat_a[kk]);    //Initialize Numbers-at-Stage, based on maturity ogive
!         Nage_mat_0[kk] = prop_mat_a[kk]*NPR_age[kk];
!         temp_1plus += NPR_age[kk];                        // keep track of one-plus females per recruit
!         temp_mat += Nage_mat_0[kk];                       // "" mature females per recruit
!       }
!
!      NPR_age[a_m] = NPR_age[a_m-1]*Sage[a_m-1]/(1-Sage[a_m-1]);   // PLUS GROUP numbers at age
!      Nage_imm_0[a_m] = NPR_age[a_m]*(1-prop_mat_a[a_m]); 		 // by definition = 0.0 given assumption that age at transition to plus group equals age at maturity
!      Nage_mat_0[a_m] = prop_mat_a[a_m]*NPR_age[a_m];		 // mature females
!
!      temp_1plus = temp_1plus + NPR_age[a_m];									// keep track of one_plus component (females)
!      temp_mat = temp_mat + Nage_mat_0[a_m];									// keep track of mature female component
!
!      b_eq = 1.0 / (S_tmp1*(temp_mat - 1));                                  // equilibrium birth rate on a per recruit basis
!      b_1 = b_eq + (b_max - b_eq) * (1 - Math.pow(depl_init,z));    // Birth rate in first year of projection, given intial depletion
!      NPR_oneplus = temp_1plus;
!
!    }      