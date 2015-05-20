!=== === +++ === === +++ === === +++ === ===
! File:   Format_module.f90
! Author: John R. Brandon
! eMail:  jbrandon at gmail
! Date :  Apr 2015
! Developed under OS:  Mac OS 10.9.5 (x86_64-apple-darwin10.8.0 (64-bit))
! Language : Fortran 90/95
! Originally compiled using: GCC GNU gfortran v 4.2 (free open source Fortran 90/95 compiler) -- https://gcc.gnu.org/
! IDE: Netbeans 8.0.2
!====== +++ === === +++ === === +++ === ===
! Purpose : Declare format statements as character strings
!
!  This is something of a hack, with the intent to organize a globally accessible list of format statements in Fortran.
!  As an example, outside code can make use of elements of this list in a read (or write) statement as follows:
!   READ(UNIT, format2) first_character, second_character
!    where: format2 is declared in this module as:
!     character(len=20), parameter :: format2 = "(a10 a10)" 
!====== +++ === === +++ === === +++ === ===

module Format_module
    implicit none
!
!    character(len=20), parameter :: format1 = "(a10 f5.3 a10 f8.3)" ! perhaps give these more expressive names? formatX is non-informative
!    character(len=20), parameter :: format2 = "(a10 a10)" ! example header for output file
!    character(len=20), parameter :: format3 = "(f10.3 f10.3)" ! example format for two column output file
!    character(len=20), parameter :: format4 = "(I10 f10.3)" ! example format for two column output file    

!10	 format(a10 f5.3 a10 f8.3)
!20	 format(a10 a10) ! header for CV_N output file
!30	 format(f10.3 f10.3) ! header for CV_N output file    

end module Format_module
