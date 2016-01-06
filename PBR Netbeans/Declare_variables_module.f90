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
! Purpose : Declare variables for PBR Tier System code
!  These variables will be available to any routine that "uses" this module
!====== +++ === === +++ === === +++ === ===
MODULE Declare_variables_module
    implicit none
!====== +++ === === +++ === === +++ === ! Variables read from input.par file
    CHARACTER REF*10,DESC*50            ! Reference tag for scenario and brief description: Naming adopted from AEP's code for File IO
    character cseed*1                   ! Should the seed for the RNG be user specified (Y) or generated from CPU clock (N)
    integer(kind = 4) :: iseed          ! Seed for RNG (if applicable)        
    integer(kind = 4) :: n_sims         ! Number of simulated population projections
    integer(kind = 4) :: n_stocks       ! Number of stocks to model (at present, either 1 or 2 stocks) 
    integer(kind = 4) :: yr_max         ! Number of years to project over   
    integer(kind = 4) :: surv_freq      ! Interval (yrs) between abundance surveys (first abundance survey in year 1) 
! TODO ? : Look into making the vectors in this list allocatable    
    real(kind = 8) :: k_1plus(2)        ! Carrying capacity (in terms of the age 1+ component of the population)    
    real(kind = 8) :: cv_n              ! The coefficient of variation (CV) associated with sampling error in abundance estimates
    real(kind = 8) :: cv_mortality(2)   ! The CV associated with sampling error for human caused mortality of marine mammals
    real(kind = 8) :: theta             ! The shape parameter for density dependence (assumed to act through birth rates)
    real(kind = 8) :: r_max             ! The maximum annual growth rate of the population
    real(kind = 8) :: f_r(2)            ! Recovery factor parameter for PBR    
    real(kind = 8) :: init_depl(2)      ! Initial depletion (fraction of carrying capacity) in first year of projections
    real(kind = 8) :: lower_tail        ! The percentile of the (log-normal) abundance estimates = N_min for PBR
    real(kind = 8) :: b_max             ! Maximum birth rate (in the absence of carrying capacity, i.e. in the limit of zero abundance)
    real(kind = 8) :: b_sex_ratio       ! Sex ratio at birth, in terms of female calves per female (i.e. 0.75 -> 3/4 calves are female)
    real(kind = 8) :: s_adult           ! Adult survival rate
    real(kind = 8) :: s_juv             ! Juvenile survival rate
    integer(kind = 4) :: a_d            ! Density dependence acts on (a_d:age_x) component, e.g. if a_d = 1, then DD on 1+ component
    integer(kind = 4) :: a_t            ! Age at transition to adult survival 
    integer(kind = 4) :: a_m            ! Age at sexual maturity  
    integer(kind = 4) :: age_x          ! Age at which individuals enter the 'plus-group' (identical for each stock at present)    
    integer(kind = 4) :: a_r            ! Age at recruitment (vulnerability) to human caused mortality - knife edge selectivity assumed    
    real(kind = 8) :: p_a1_s1           ! Percentage of stock_1 in subarea_1
    real(kind = 8) :: p_a2_s1           ! Percentage of stock_1 in subarea_2
    real(kind = 8) :: p_a2_s2           ! Percentage of stock_2 in subarea_2
    real(kind = 8) :: p_a3_s2           ! Percentage of stock_2 in subarea_3
    real(kind = 8) :: p_a4_s2           ! Percentage of stock_2 in subarea_4  
    real(kind = 8) :: omega(1:4)          ! Relative vulnerability of animals by area (regardless of stock ID)
    real(kind = 8) :: m_bias            ! Bias in mortality estimates (1.0 = no bias; > 1.0 = estimated mortality less than actual mortality) : Wade (1998) trial 1
    real(kind = 8) :: n_bias            ! Bias in abundance estimates (1.0 = no bias; < 1.0 = estimated abundance more than actual abundance) : Wade (1998) trial 2
    real(kind = 8) :: r_bias            ! Bias in R_max (1.0 = no bias; < 1.0 = true R_max is less than assumed by default PBR value) : Wade (1998) trial 3
    real(kind = 8) :: cv_n_true         ! Actual CV of abundance (Eqn 3) Wade (1998) used to generate N_hat : Wade (1998) trial 4
!    real(kind = 8) :: cv_m_true         ! True CV of mortality : Not sure this is needed, could just change cv_mortality and have seperate input files: Wade (1998) trial 5    
    integer(kind = 4) :: tier           ! Data tier
    integer(kind = 4) :: n_yrs_avg      ! Window of years to average across
    character determ_pbr*1              ! N = Interannual stochastic human caused mortality; Y = Deterministic (e.g. equal exactly to PBR or some fraction therein) 
    real(kind = 8) :: prop_catch_fem    ! Proportion of catch that is female
    real(kind = 8) :: F_rate            ! Constant fishing mortality rate (DEBUGGING)
!    character dd_component*1
!====== +++ === === +++ === === +++ === ! End variables read from input.par
!====== +++ === === +++ === === +++ === ! These are probably useful also as global variables     
    real(kind = 8), allocatable :: selectivity(:)   ! Selectivity at age (identical for each stock at present)  
    real(kind = 8), allocatable :: S_age(:)         ! Vector of survival-at-age (identical for each stock at present)
    real(kind = 8), allocatable :: prop_mat_age(:)  ! Proportion mature at age (all zeros until a_m) (identical for each stock at present)
    real(kind = 8), allocatable :: transition_matrix(:, :) ! Transition matrix (with survival and birth rates)    
                                                    ! Variables below can be stock specific given different initial depletion conditions for each stock
    real(kind = 8) :: init_depl_i                   ! Initial depletion for stock i     
    real(kind = 8), allocatable :: NPR_age(:)       ! Numbers-at-age-per-recruit
    real(kind = 8), allocatable :: NPR_age_tmp(:)   ! Vector of numbers at age, used to solve for initial human caused mortality rate 
    real(kind = 8) :: NPR_oneplus                   ! Total numbers of ages 1+ per female recruit (zero human caused mortality)    
!    real(kind=8) :: initial_oneplus                ! Size of female 1+ component on per-recruit basis given F_init    
!====== +++ === === +++ === === +++ === ! Variables below might not need to be in this module
    real(kind=8) :: b_eq                            ! Equilibrium birth rate (at carrying capacity)
    real(kind=8) :: b_init                          ! Birth rate in first year of projections (a function of initial depletion in abundance)
    real(kind=8) :: fecundity_max                   ! Female calves per female, for eigen analysis of projection matrix 
   
    integer(kind = 4) :: stock_i                    ! Stock ID number   

    real(kind = 8) :: N_best                          ! Mean of the abundance estimate, E[N], given sampling error (CV_N)
    real(kind = 8), allocatable :: N_min_yr_sim(:, :) ! N_min estimates by year and simulation number
    
    real(kind = 8) :: sum_NPR       ! The sum across ages in the NPR vector (including age zero calves)
           
    real(kind = 8), allocatable :: Female_age(:)        ! Vector with numbers-at-age for females 
    real(kind = 8), allocatable :: Male_age(:)          ! Vector with numbers-at-age for males 
    real(kind = 8), allocatable :: Nage_imm_0(:)        ! numbers-at-age that are immature
    real(kind = 8), allocatable :: Nage_mat_0(:)        ! numbers-at-age that are mature
    real(kind = 8), allocatable :: prop_NPR(:)          ! Proportions in each age of the NPR vector
    
    real(kind = 8), allocatable :: N_plus(:)            ! Vector of 1+ population size for stock 'i' across all years 
    real(kind = 8), allocatable :: N_tot(:)             ! Total (0+) population size each year for stock 'i' 
    real(kind = 8), allocatable :: N_calf(:)            ! Vector of calf production for stock 'i' across all years
    real(kind = 8), allocatable :: N_age(:)             ! Numbers-at-age

!====== +++ === === +++ === === +++ === 
    contains
!###### +++ ### ### +++ ### ### +++ ###     
    subroutine initialize_global_vars()
!====== +++ === === +++ === === +++ ===   
!###### +++ ### ### +++ ### ### +++ ###         
! Initialize variables declared in this module. Set values zero. 
! Only do this for those variables with values not read by input file.            
!====== +++ === === +++ === === +++ ===
        b_eq = 0.0d0                 ! Equilibrium birth rate (at carrying capacity)
        b_init = 0.0d0               ! Birth rate in first year of projections (a function of initial depletion in abundance)
        fecundity_max = 0.0d0        ! Female calves per female, for eigen analysis of projection matrix 
        
        stock_i = 0                  ! Stock ID number (integer)
        
        N_best = 0.0d0               ! Expectation of an abundance estimate
        N_min_yr_sim = 0.0d0                ! See above for remaining variable definitions

        sum_NPR = 0.0d0
        Female_age = 0.0d0
        Male_age = 0.0d0     
        Nage_imm_0 = 0.0d0        
        Nage_mat_0 = 0.0d0        
        prop_NPR = 0.0d0        
        N_plus = 0.0d0 
        N_tot = 0.0d0
        N_calf = 0.0d0
        N_age = 0.0d0
       
        return
    end subroutine initialize_global_vars
! *
END MODULE Declare_variables_module
