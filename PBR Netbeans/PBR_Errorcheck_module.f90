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
    use Declare_variables_module ! Global variable declarations
    implicit none
    contains
    
    integer(kind = 4) function error_check_input()
! Do some preliminary error checking on the input file 
! For example, the number of stocks entered should be either 1 or 2 (no less and no more)    
        error_check_input = 0

        if (n_stocks > 2 .or. n_stocks < 1) then      ! Check that number of stocks is supported
            print *, "ERROR : Code only currently supports either one or two stocks. Check input file, input.par"
            error_check_input = -9
        end if

        if (p_a1_s1 + p_a2_s1 .ne. 1.0) then
            print *, "ERROR : Percentages for stock 1 in each area do not add to 1.0. Check input file, input.par"
            error_check_input = -8            
        end if

        if (p_a2_s2 + p_a3_s2 + p_a4_s2 .ne. 1.0) then
            print *, "ERROR : Percentages for stock 2 in each area do not add to 1.0. Check input file, input.par"
            error_check_input = -8            
        end if

        return
    end function error_check_input

    ! TODO : Add functionality to check that solution for S_juv < S_adult
    
END MODULE PBR_Errorcheck_module
