!=== === +++ === === +++ === === +++ === ===
! Author: John R. Brandon
! eMail:  jbrandon at gmail
! Date :  Apr 2015
! OS   :  Mac OS 10.9.5 (x86_64-apple-darwin10.8.0 (64-bit))
! Language : Fortran 90/95
!====== +++ === === +++ === === +++ === ===
! Purpose : This module contains routines for file input/output (IO)
!====== +++ === === +++ === === +++ === ===
module PBR_FileIO_Module
    
    use Declare_variables_module 
    
    implicit none
    
    contains

    subroutine read_inits()! (CV_N, CV_MORTALITY, THETA, R_MAX, F_R, INIT_DEPL, &
         !LOWER_TAIL, YR_MAX, SURV_FREQ, KK, IPAR)

!     Read parameter file and initializes variables
       integer :: INPUT_FILE
       DATA INPUT_FILE /7/
       
       OPEN (INPUT_FILE, FILE='INPUT.PAR', STATUS='OLD')

801    FORMAT(16X, I6)
802    FORMAT(16X, F6.3)
803    FORMAT(16X, A6)
804    FORMAT(16X, A1)

!     Read in model parameters/options, checking they are within allowed range
       READ (INPUT_FILE,'(T37, A, /A)') REF,DESC !,PARFIL,MATFIL
       print *, "REF : ", REF
       print *, "DESC : ", DESC
       
       READ (INPUT_FILE, *) ! skip lines with comments / blanks in header of input file
       READ (INPUT_FILE, *)
       READ (INPUT_FILE, 804) cseed    ! Will the user be supplying a seed for the random number generator (RNG)?
       READ (INPUT_FILE, 801) iseed    ! Seed for RNG, ignored if cseed = N
       READ (INPUT_FILE, 801) n_stocks
       READ (INPUT_FILE, 801) YR_MAX 
       READ (INPUT_FILE, 801) SURV_FREQ 
       READ (INPUT_FILE, 801) KK 
       
       READ (INPUT_FILE, 802) CV_N 
       READ (INPUT_FILE, 802) CV_MORTALITY 
       READ (INPUT_FILE, 802) THETA 
       READ (INPUT_FILE, 802) R_MAX 
       READ (INPUT_FILE, 802) F_R 
       READ (INPUT_FILE, 802) INIT_DEPL 
       READ (INPUT_FILE, 802) LOWER_TAIL 
       READ (INPUT_FILE, 802) b_max 
       READ (INPUT_FILE, 802) S_adult
       READ (INPUT_FILE, 802) S_juv
       READ (INPUT_FILE, 801) a_m       
       

       WRITE(*,*) "cseed: ", cseed ! DEBUGGING
       WRITE(*,*) "iseed: ", iseed ! DEBUGGING
       WRITE(*,*) "n_stocks", n_stocks
       WRITE(*,*) "YR_MAX", YR_MAX
       WRITE(*,*) "SURV_FREQ", SURV_FREQ
       WRITE(*,*) "KK", KK
       WRITE(*,*) "CV_N", CV_N
       WRITE(*,*) "CV_MORTALITY", CV_MORTALITY
       WRITE(*,*) "THETA", THETA
       WRITE(*,*) "R_MAX", R_MAX
       WRITE(*,*) "F_R", F_R
       WRITE(*,*) "INIT_DEPL", INIT_DEPL
       WRITE(*,*) "LOWER_TAIL", LOWER_TAIL
       WRITE(*,*) "b_max", b_max
       WRITE(*,*) "S_adult", S_adult
       WRITE(*,*) "S_juv", S_juv       
       WRITE(*,*) "a_m", a_m       
       
       CLOSE(INPUT_FILE)
       
      RETURN
      END subroutine read_inits

END MODULE PBR_FileIO_Module
!# Input parameters for PBR Tier system simulations -- base case 
!yr_max          100     # Number of years to project over
!surv_freq       4       # interval (yrs) between abundance surveys (first abundance survey in year 1)
!KK              10000   # Carrying capacity
!cv_n            0.20    # CV(N) , will be transformed to the standard deviation in log-space
!cv_mortality    0.30    # Variability in mortality, used to generate random normal deviate with mean = PBR
!theta           1.0     # Density dependence shape parameter
!R_max           0.04    # Default for cetaceans 
!F_r             0.50    # Base case recovery factor
!init_depl       0.30    # Initial depletion (Abundance as a fraction of carrying capacity)
!lower_tail      0.20    # Percentile of log-normal distribution, used to calculate N.min, given CV(N)
!surv.yr.tmp     1       # Counter used in assigning survey years
!B_max           0.5     # Maximum birth rate (in the absence of density dependence)
!S_adult         0.95    # Adult (mature) survival rate
!S_juv           0.80    # Juvenile (immature) survival rate
!a_m             10      # Age at sexual maturity  