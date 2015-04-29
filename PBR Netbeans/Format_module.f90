!     
! File:   Format_module.f90
! Author: johnbrandon
!
! Created on April 25, 2015, 6:52 PM
!

module Format_module
    implicit none

    character(len=20), parameter :: format1 = "(a10 f5.3 a10 f8.3)"
    character(len=20), parameter :: format2 = "(a10 a10)" ! example header for output file
    character(len=20), parameter :: format3 = "(f10.3 f10.3)" ! example format for two column output file
    character(len=20), parameter :: format4 = "(I10 f10.3)" ! example format for two column output file    

!10	 format(a10 f5.3 a10 f8.3)
!20	 format(a10 a10) ! header for CV_N output file
!30	 format(f10.3 f10.3) ! header for CV_N output file    

end module Format_module
