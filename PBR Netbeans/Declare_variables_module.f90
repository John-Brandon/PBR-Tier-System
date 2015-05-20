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
!  It is not necessarily great programming practice, but while lazy, it does make the program code easier to read in some cases.
!====== +++ === === +++ === === +++ === ===

MODULE Declare_variables_module
    implicit none
!    INTEGER IOUT, IN1       ! File IO 
!    DATA IOUT/8/,IN1/10/    ! File unit numbers for file IO
!====== +++ === === +++ === === +++ === ! Variables read from input.par file
    CHARACTER REF*10,DESC*50            ! Reference tag for scenario and brief description: Naming adopted from AEP's code for File IO
    character cseed*1                   ! Should the seed for the RNG be user specified (Y) or generated from CPU clock (N)
    integer(kind = 4) :: iseed          ! Seed for RNG (if applicable)        
    integer(kind = 4) :: n_stocks       ! Number of stocks to model (at present, either 1 or 2 stocks) 
    integer(kind = 4) :: yr_max         ! Number of years to project over   
    integer(kind = 4) :: surv_freq(2)   ! Interval (yrs) between abundance surveys (first abundance survey in year 1)    
    real(kind = 8) :: k_1plus(2)        ! Carrying capacity (in terms of the age 1+ component of the population)    
    real(kind=8) :: cv_n(2)             ! The coefficient of variation (CV) associated with sampling error in abundance estimates
    real(kind=8) :: cv_mortality(2)     ! The CV associated with sampling error for human caused mortality of marine mammals
    real(kind=8) :: theta               ! The shape parameter for density dependence (assumed to act through birth rates)
    real(kind=8) :: r_max               ! The maximum annual growth rate of the population
    real(kind=8) :: f_r(2)              ! Recovery factor parameter for PBR    
    real(kind=8) :: init_depl(2)        ! Initial depletion (fraction of carrying capacity) in first year of projections
    real(kind=8) :: lower_tail          ! The percentile of the (log-normal) abundance estimates = N_min for PBR
    real(kind=8) :: b_max               ! Maximum birth rate (in the absence of carrying capacity, i.e. in the limit of zero abundance)
    real(kind=8) :: b_sex_ratio         ! Sex ratio at birth
    real(kind=8) :: s_adult             ! Adult survival rate
    real(kind=8) :: s_juv               ! Juvenile survival rate
    integer(kind = 4) :: a_m            ! Age at sexual maturity  
    integer(kind = 4) :: a_r            ! Age at recruitment (vulnerability) to human caused mortality - knife edge selectivity assumed
!====== +++ === === +++ === === +++ === ! End variables read from input.par
!====== +++ === === +++ === === +++ === ! These are probably useful also as global variables     
    integer(kind = 4) :: age_x                      ! Age at which individuals enter the 'plus-group' (identical for each stock at present)
    real(kind = 8), allocatable :: selectivity(:)   ! Selectivity at age (identical for each stock at present)
    real(kind = 8), allocatable :: S_age(:)         ! Vector of survival-at-age (identical for each stock at present)
    real(kind = 8), allocatable :: prop_mat_age(:)  ! Proportion mature at age (all zeros until a_m) (identical for each stock at present)
                                                    ! Variables below can be stock specific given different initial depletion conditions for each stock
    real(kind = 8) :: init_depl_i                   ! Initial depletion for stock i     
    real(kind = 8), allocatable :: NPR_age(:)       ! Numbers-at-age-per-recruit
    real(kind = 8), allocatable :: NPR_age_tmp(:)   ! Vector of numbers at age, used to solve for initial human caused mortality rate 
    real(kind = 8) :: NPR_oneplus                   ! Total numbers of ages 1+ per female recruit (zero human caused mortality)    
!    real(kind=8) :: initial_oneplus                 ! Size of female 1+ component on per-recruit basis given F_init    
!====== +++ === === +++ === === +++ === ! These maybe NOT so much     
    real(kind=8) :: b_eq                            ! Equilibrium birth rate (at carrying capacity)
    real(kind=8) :: b_init                          ! Birth rate in first year of projections (a function of initial depletion in abundance)
!    integer(kind = 4) :: first_yr       ! First year of projection

    integer(kind = 4) :: ii, jj, kk, ll, mm ! Counters
    
    integer(kind = 4) :: stock_i            ! Stock ID number   

    real(kind = 8) :: N_best                ! Mean of the abundance estimate 
    real(kind = 8) :: N_min                 ! 
    
    real(kind = 8) :: z_variate             ! Standard normal random variate, i.e. ~ N(mu = 0, sigma = 1)

!    real(kind = 8) :: f_init ! Initial human caused mortality rate
!    real(kind = 8) :: objf_f_init     ! Objective function value for solving for initial human caused mortality rate

!    real(kind = 8) :: temp_1plus    !
!    real(kind = 8) :: temp_mat      !
!    real(kind = 8) :: temp_mature   ! 
 
    real(kind = 8) :: sum_NPR       ! The sum across ages in the NPR vector (including age zero calves)
    
!    real(kind = 8) :: rec_init      ! Initial recruitment conditioned on Finit (in terms of females)
!    real(kind = 8) :: pred_rec      ! Predicted recruitment given some value for F_init (used to solve for equilibrium F_init)
        
    real(kind = 8), allocatable :: Female_age(:)        ! Vector with numbers-at-age for females 
    real(kind = 8), allocatable :: Male_age(:)          ! Vector with numbers-at-age for males 
    real(kind = 8), allocatable :: Nage_imm_0(:)        ! numbers-at-age that are immature
    real(kind = 8), allocatable :: Nage_mat_0(:)        ! numbers-at-age that are mature
    real(kind = 8), allocatable :: prop_NPR(:)          ! Proportions in each age of the NPR vector
    
    real(kind = 8), allocatable :: Nplus(:)    ! Vector of 1+ population size for stock 'i' across all years 
    real(kind = 8), allocatable :: Ntot(:)        ! Total (0+) population size each year for stock 'i' 
    real(kind = 8), allocatable :: N_calf(:)       ! Vector of calf production for stock 'i' across all years
    real(kind = 8), allocatable :: depl_i_t(:, :)         ! Depletion for each stock each year
    real(kind = 8), allocatable :: Nage(:)       ! Numbers-at-age-per-recruit    
!    real(kind = 8), allocatable :: selectivity_i(:, :)    ! Selectivity at age for stock 'i'
!    
!    real(kind = 8), allocatable :: N_oneplus_i_t(:, :, :)    ! Vector of 1+ population size for stock 'i' across all years 
!    real(kind = 8), allocatable :: N_tot_i_t(:, :, :)        ! Total (0+) population size each year for stock 'i' 
!    real(kind = 8), allocatable :: N_calf_i_t(:, :, :)       ! Vector of calf production for stock 'i' across all years
!    real(kind = 8), allocatable :: depl_i_t(:, :, :)         ! Depletion for each stock each year
    
END MODULE Declare_variables_module
