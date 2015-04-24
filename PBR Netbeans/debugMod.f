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

      end module debug
     