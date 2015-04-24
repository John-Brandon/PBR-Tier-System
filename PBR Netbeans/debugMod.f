!#######################################################    
! File:   debugMod.f90
! Author: johnbrandon
!
! Created on April 23, 2015, 6:38 PM
!#######################################################
! Temporary module for developing and debugging code
!#######################################################


      module debug

         contains 
    
         real(kind = 8) function debugz(a_m)
            implicit none
            integer(kind = 4), intent(in) :: a_m ! intent(in) specifies for compiler that 'a_m' is being passed into function
            debugz = real(a_m) / 4 ! note the use of real(a_m), otherwise FORTRAN will round the quotient of two integers
            print *, "a_m: ", a_m 
            print *, "a_m / 4: ", a_m / 4
            print *, "Hello from debugz"
            return
         end function debugz
        
      end module debug
