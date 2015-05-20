!#######################################################    
! File:   debugMod.f90
! Author: johnbrandon
!
! Created on April 23, 2015, 6:38 PM
!#######################################################
! Temporary module for developing and debugging code
!#######################################################
      module debug
    
         implicit none
         
         contains 
!#######################################################    
         real(kind = 8) function debugz(a_m_tmp)
            implicit none
            integer(kind = 4), intent(in) :: a_m_tmp ! intent(in) specifies for compiler that 'a_m' is being passed into function
            debugz = real(a_m_tmp) / 4 ! note the use of real(a_m), otherwise FORTRAN will round the quotient of two integers
            print *, "a_m_tmp: ", a_m_tmp 
            print *, "a_m_tmp / 4: ", real(a_m_tmp) / 4
            print *, "Hello from debugz"
            return
         end function debugz
         
!#######################################################     
	subroutine test_cpu_time(start, finish, aa, bb)
	
	 real start, finish
	 real aa, bb
	 call cpu_time(start)
! put code to test here
! Want to benchmark matrix multiplication in FORTRAN vs. that in R (for age-structured model)
	  
	 call cpu_time(finish)
	 print '("Time = ",f6.3," seconds.")',finish-start
	return
	end subroutine test_cpu_time

!#######################  
        real(kind = 8) function newtons_root()
!#######################
! Using the Newton-Raphson method, find the root of a function known to lie in the interval [x1,x2].
!  The root (rtnewt) will be refined until its accuracty is known within +/- (xacc). (funcd) is a user
!  supplied subroutine that returns both the function value and the first derivative of the function at the point (x).'
!  - description from Numerical Recipes book
!  In this case going to be solving for human caused mortality rate which would leave the population at initial depletion level
!#######################
        implicit none
! Declare local variables        
        integer :: jcount  ! counter specific to this function
        integer :: jmax    ! maximum number of iterations
        real(kind = 8) :: delta_Finit
        real(kind = 8) :: rtnewt
        real(kind = 8) :: f1_newt
        real(kind = 8) :: df_newt
        real(kind = 8) :: dx_newt
        real(kind = 8) :: f_init
!       real(kind = 8), external :: Initial_F
!
        f_init = 0.
        f1_newt = 0.
        objf_newt = 0.
!
        jmax = 20
        delta_Finit = 0.000001
        rtnewt = 0.50             ! Initial guess for human caused mortality rate in year zero

        do jcount = 1, jmax
              print *, "jcount: ", jcount
              f_init = rtnewt
              objf_newt = Initial_F(f_init)              ! Subroutine call
              f1_newt = objf_newt	
              f_init = rtnewt + delta_Finit
              objf_newt = Initial_F(f_init)              ! Subroutine call
              df_newt = (objf_newt - f1_newt) / delta_Finit
              dx_newt = f1_newt / df_newt
              rtnewt = rtnewt - dx_newt
              print *, "rtnewt: ", rtnewt              
        end do

        print *, "Hello from newtons_root"
        print *, "Max iterations: ", jmax
        print *, "Iterations completed: ", jcount - 1        
        print *, "newtons_root: ", newtons_root
        newtons_root = rtnewt
        return
        
        end function newtons_root
      
      end module debug
!end end end end end end end end end end end end end end end end
