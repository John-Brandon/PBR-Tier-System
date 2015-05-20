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
        print *, "error_check_input: ", error_check_input
        if (n_stocks > 2 .or. n_stocks < 1) then      ! Check that number of stocks is supported
            print *, "ERROR : Code only currently supports either one or two stocks. Check input file, input.par"
            error_check_input = -9
        end if
            
        return
    end function error_check_input

END MODULE PBR_Errorcheck_module
