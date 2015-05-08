!=== === +++ === === +++ === === +++ === ===
! File:   Declare_variables_module.f90
! Author: John R. Brandon
! eMail:  jbrandon at gmail
! Date :  Apr 2015
! Developed under OS:  Mac OS 10.9.5 (x86_64-apple-darwin10.8.0 (64-bit))
! Language : Fortran 90/95
! Originally compiled using: GCC GNU gfortran v 4.2 (free open source Fortran 90/95 compiler) -- https://gcc.gnu.org/
! IDE: Netbeans 8.0.2
!====== +++ === === +++ === === +++ === ===
! Purpose : Declare variables for PBR Tier program
!  These variables will be available to any routine that "uses" this module
!  This approach is analagous with declaring global variables (or common block in Fortran 77). 
!  It is not necessarily great programming practice, but while lazy, it does make the program code easier to read.
!====== +++ === === +++ === === +++ === ===

MODULE Declare_variables_module
    implicit none
      
    CHARACTER REF*10,DESC*50 ! Reference tag for scenario and brief description: Naming adopted from AEP's code for File IO
    character cseed*1 ! Should the seed for the RNG be user specified (Y) or generated from CPU clock (N)
    INTEGER IOUT, IN2 ! File IO
    DATA IOUT/8/,IN2/10/    ! File unit numbers for file IO
         
    real(kind=8) :: cv_n            ! The coefficient of variation (CV) associated with sampling error in abundance estimates
    real(kind=8) :: cv_mortality    ! The CV associated with sampling error for human caused mortality of marine mammals
    real(kind=8) :: lower_tail      ! The percentile of the (log-normal) abundance estimates = N_min for PBR
    real(kind=8) :: f_r             ! Recovery factor parameter for PBR
    real(kind=8) :: theta           ! The shape parameter for density dependence (assumed to act through birth rates)
    real(kind=8) :: r_max           ! The maximum annual growth rate of the population
    real(kind=8) :: b_max           ! Maximum birth rate (in the absence of carrying capacity, i.e. in the limit of zero abundance)
    real(kind=8) :: b_eq            ! Equilibrium birth rate (at carrying capacity)
    real(kind=8) :: b_1             ! Birth rate in first year of projections (a function of initial depletion in abundance)
    real(kind=8) :: s_adult         ! Adult survival rate
    real(kind=8) :: s_juv           ! Juvenile survival rate
    real(kind=8) :: init_depl       ! Initial depletion (fraction of carrying capacity) in first year of projections
    real(kind=8) :: initial_oneplus ! Size of  female 1+ component on per-recruit basis given Finit    
    integer(kind = 4) :: first_yr
    integer(kind = 4) :: yr_max
    integer(kind = 4) :: surv_freq 
    integer(kind = 4) :: a_m    ! age at maturity
    integer(kind = 4) :: ii, jj, kk, ll, mm ! counters
    integer(kind = 4) :: xx, nn
    integer(kind = 4) :: npr, maxage, age_x, n_stocks 
    integer(kind = 4) :: iseed
    integer(kind = 4) :: k_1plus ! carrying capacity (in terms of the age 1+ component of the population)    

    real(kind = 8) :: aa, bb
    real(kind = 8) :: start, finish
    real(kind = 8) :: N_best, z_variate
    real(kind = 8) :: N_min
    real(kind = 8) :: x_min, x_max, step_CV
    real(kind = 8) :: yy
    real(kind = 8) :: delt_s ! delt_s is the difference between adult and juvenile survival rates
!    real(kind = 8) :: f_init ! Initial human caused mortality rate
    real(kind = 8) :: objf_newt ! Objective function value for solving for initial human caused mortality rate using Newton's method  
    real(kind = 8) :: S_tmp1       ! Temporary vector for survival rates -- don't think this is needed here. 
    real(kind = 8) :: temp_1plus, temp_mat, NPR_oneplus
    real(kind = 8) :: sum_NPR     ! The sum across ages in the NPR vector
    
    real(kind = 8) :: temp_mature ! Not sure if this is used ***
    real(kind = 8) :: rec_init ! initial recruitment conditioned on Finit (in terms of females)
    real(kind = 8) :: pred_rec
        
    real(kind = 8), allocatable :: Female_age(:) ! Vector with numbers-at-age for females 
    real(kind = 8), allocatable :: Male_age(:) ! Vector with numbers-at-age for males 
    real(kind = 8), allocatable :: S_age(:)                ! Vector of survival-at-age
!   real(kind = 8), allocatable :: P_age(:)                ! Vector of rate of maturation-at-age
    real(kind = 8), allocatable :: prop_mat_a(:)          ! Proportion mature at age (all zeros until a_m)
    real(kind = 8), allocatable :: NPR_age(:)             ! Numbers-at-age-per-recruit
    real(kind = 8), allocatable :: Nage_imm_0(:)          ! numbers-at-age that are immature
    real(kind = 8), allocatable :: Nage_mat_0(:)          ! numbers-at-age that are mature
    real(kind = 8), allocatable :: prop_NPR(:)            ! Proportions in each age of the NPR vector
    real(kind = 8), allocatable :: Nage(:)                ! Vector of numbers at age, used to solve for initial human caused mortality rate
    real(kind = 8), allocatable :: selectivity(:)         ! Selectivity at age
    
END MODULE Declare_variables_module
