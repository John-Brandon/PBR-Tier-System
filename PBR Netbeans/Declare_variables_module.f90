!=== === +++ === === +++ === === +++ === ===
! Author: John R. Brandon
! eMail:  jbrandon at gmail
! Date :  Apr 2015
! OS   :  Mac OS 10.9.5 (x86_64-apple-darwin10.8.0 (64-bit))
! Language : Fortran 90/95
!====== +++ === === +++ === === +++ === ===
! Purpose : Declare variables for PBR Tier program
!  These variables will be available to any routine that "uses" this module
!====== +++ === === +++ === === +++ === ===

MODULE Declare_variables_module
    implicit none
    
    real(kind=8) :: CV_N, CV_MORTALITY, THETA, R_MAX, F_R 
    real(kind=8) :: INIT_DEPL, LOWER_TAIL, b_max, S_adult, S_juv
    integer(kind = 4) :: YR_MAX, SURV_FREQ, a_m

    integer(kind = 4) :: ii, jj, kk, ll ! counters
    integer(kind = 4) :: xx, nn
    integer(kind = 4) :: npr, maxage, n_stocks 

    real(kind = 8) :: aa, bb
    real(kind = 8) :: start, finish
    real(kind = 8) :: N_best, z_variate
    real(kind = 8) :: N_min
    real(kind = 8) :: x_min, x_max, step_CV
    real(kind = 8) :: yy
    real(kind = 8) :: delt_s ! delt_s is the difference between adult and juvenile survival rates
         
END MODULE Declare_variables_module
