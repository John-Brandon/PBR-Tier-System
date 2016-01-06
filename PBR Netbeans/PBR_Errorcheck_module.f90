!====== +++ === === +++ === === +++ === ===     
! File:   PBR_Errorcheck_module.f90
! Author: John R. Brandon
! eMail:  jbrandon at gmail
! Date :  Apr 2015
! Developed under OS:  Mac OS 10.9.5 (x86_64-apple-darwin10.8.0 (64-bit))
! Language : Fortran 90/95
! Originally compiled using: GCC GNU gfortran v 4.2 (free open source Fortran 90/95 compiler) -- https://gcc.gnu.org/
! IDE: Netbeans 8.0.2
!====== +++ === === +++ === === +++ === ===
MODULE PBR_Errorcheck_module
    use Declare_variables_module ! Global variable declarations -- mostly contains the variables read from the input file. 
    implicit none
    contains
!====== +++ === === +++ === === +++ === ===    
    integer(kind = 4) function error_check_input()
!====== +++ === === +++ === === +++ === ===    
! Do some preliminary error checking on the input file 
! For example, the number of stocks entered should be either 1 or 2 (no less and no more)    
!  This error checking code is a work in progress. It is certainly not yet as comprehensive as it could / should be.     
!====== +++ === === +++ === === +++ === ===    
        error_check_input = 0
        
        if (a_d .ne. 1) then      ! Check density dependent component
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Currently code only supports density dependence on age 1+ (i.e. a_d = 1). Check input file."
          error_check_input = -9
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        
        if (n_stocks > 2 .or. n_stocks < 1) then      ! Check that number of stocks is supported
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Code only currently supports either one or two stocks. Check input file."
          error_check_input = -9
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if

        if (p_a1_s1 + p_a2_s1 .ne. 1.0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Percentages for stock 1 in each area do not add to 1.0. Check input file."
          error_check_input = -8            
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if

        if (p_a2_s2 + p_a3_s2 + p_a4_s2 .ne. 1.0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Percentages for stock 2 in each area do not add to 1.0. Check input file."
          error_check_input = -8            
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        
        if (a_m > age_x .or. a_t > age_x) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : age_x must be greater than or equal to a_m and a_t. Check input file."
          error_check_input = -7            
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        
        if (a_t < 1 .or. a_t > age_x) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Condition (1 < a_t <= age_x) not met. Check input file."
          error_check_input = -7            
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        
        if (tier < 1 .or. tier > 4) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Condition (0 < tier < 3) not met. Check input file."
          error_check_input = -6            
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        
        if (tier .eq. 3 .and. n_yrs_avg > 8) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Tier 3 can not average over more than 8 years. Check input file. Change (n_yrs_avg)"
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if
        
        if (prop_catch_fem < 0.d0 .or. prop_catch_fem > 1.d0) then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : Condition (0 < prop_catch_fem < 1) not met. Check input file."
          error_check_input = -6            
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if        
        
        if (determ_pbr /= "Y" .and. determ_pbr /= "N") then
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print *, "ERROR : determ_pbr must be 'Y' or 'N'. Check input file."
          error_check_input = -6            
          print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        end if            
        
        
        return
    end function error_check_input

    ! TODO : Add functionality to check that solution for S_juv < S_adult
    
END MODULE PBR_Errorcheck_module
