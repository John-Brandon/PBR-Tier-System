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
    
    use Declare_variables_module ! Access global variables
    
    implicit none
    
    contains

    subroutine read_inits()! (cseed, iseed, CV_N, CV_MORTALITY, THETA, R_MAX, F_R, INIT_DEPL, &
         !LOWER_TAIL, YR_MAX, SURV_FREQ, KK, IPAR)

!     Read parameter file and initializes variables
       integer :: input_file
       data input_file /7/   ! move this into Declare_variables_module?
       
       open (input_file, file = 'INPUT.PAR', status = 'OLD')

801    FORMAT(16X, I6) ! move these into Format_module?
802    FORMAT(16X, F8.0)
!803    FORMAT(16X, A6)
804    FORMAT(16X, A1)

!     Read in model parameters/options, checking they are within allowed range
       read (input_file,'(T37, A, /A)') REF,DESC !,PARFIL,MATFIL
       print *, "REF : ", REF
       print *, "DESC : ", DESC
       
       read (input_file, *) ! Skip lines with comments / blanks in header of input file
       read (input_file, *) ! Each read statement in Fortran moves to next line before reading
       read (input_file, 804) cseed    ! Will the user be supplying a seed for the random number generator (RNG)?
       read (input_file, 801) iseed    ! Seed for RNG, ignored if cseed = N
       read (input_file, 801) n_sims
       read (input_file, 801) n_stocks
       read (input_file, 801) yr_max 
       read (input_file, 801) surv_freq(1)
       read (input_file, 801) surv_freq(2)       
       read (input_file, 802) k_1plus(1)  ! Carrying capacity for stock 1
       read (input_file, 802) k_1plus(2)  ! Carrying capacity for stock 2       
       
       read (input_file, 802) cv_N(1) 
       read (input_file, 802) cv_N(2)       
       read (input_file, 802) cv_mortality(1)
       read (input_file, 802) cv_mortality(2)       
       read (input_file, 802) theta 
       read (input_file, 802) r_max 
       read (input_file, 802) f_r(1) 
       read (input_file, 802) f_r(2)       
       read (input_file, 802) init_depl(1)
       read (input_file, 802) init_depl(2)
       read (input_file, 802) lower_tail
       read (input_file, 802) b_max 
       read (input_file, 802) b_sex_ratio       
       read (input_file, 802) S_adult
       read (input_file, 802) S_juv
       read (input_file, 801) a_t       
       read (input_file, 801) a_m
       read (input_file, 801) age_x
       read (input_file, 801) a_r
       read (input_file, 802) p_a1_s1
       read (input_file, 802) p_a2_s1       
       read (input_file, 802) p_a2_s2
       read (input_file, 802) p_a3_s2
       read (input_file, 802) p_a4_s2       

! Output to screen for checking
       write(*,*) "cseed: ", cseed ! DEBUGGING
       write(*,*) "iseed: ", iseed ! DEBUGGING
       write(*,*) "n_sims", n_sims
       write(*,*) "n_stocks", n_stocks
       write(*,*) "YR_MAX", YR_MAX
       write(*,*) "SURV_FREQ", SURV_FREQ
       write(*,*) "k_1plus(1)", k_1plus(1)
       write(*,*) "k_1plus(2)", k_1plus(2)       
       write(*,*) "CV_N", CV_N
       write(*,*) "CV_MORTALITY", CV_MORTALITY
       write(*,*) "THETA", THETA
       write(*,*) "R_MAX", R_MAX
       write(*,*) "F_R", F_R
       write(*,*) "INIT_DEPL", INIT_DEPL
       write(*,*) "LOWER_TAIL", LOWER_TAIL
       write(*,*) "b_max", b_max
       write(*,*) "S_adult", S_adult
       write(*,*) "S_juv", S_juv       
       write(*,*) "a_t", a_t
       write(*,*) "a_m", a_m
       write(*,*) "age_x", age_x
       write(*,*) "a_r", a_r        
       write(*,*) "p_a1_s1", p_a1_s1               
       write(*,*) "p_a2_s1", p_a2_s1               
       write(*,*) "p_a2_s2", p_a2_s2               
       write(*,*) "p_a3_s2", p_a3_s2               
       write(*,*) "p_a4_s2", p_a4_s2                                    

       close(input_file)
       
       return
      end subroutine read_inits

end module PBR_FileIO_Module
