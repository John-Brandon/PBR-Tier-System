MODULE initialize_pop    
  use Declare_variables_module ! Access global variables      
  implicit none                ! This rule now applies to each procedure below    
!====== +++ === === +++ === === +++ ===
  contains        
!   * subroutine assign_par_vectors() *    
! Assign values to survival, proportion mature and selectivity at age vectors
!    
!   * subroutine calc_NPR_age() *    
! Does numbers-per-recruit calculations (NPR). 
! Used to calculate the unexploited and exploited age-structure in the first year of the projections                 
!    
!   * function Initial_F() *
! Given numbers per-recruit at equilibrium - apply initial human caused mortality rate (F_init)
! This function is called by Brent's method (file: BRENT.f90) to find f_init that results in a stable age structure, given 
!  the initial depletion level of the population relative to carrying capacity, i.e. the root of the objective function value
!  returned by Initial_F()    
!====== +++ === === +++ === === +++ === 

!====== +++ === === +++ === === +++ ===   
!###### +++ ### ### +++ ### ### +++ ###         
  subroutine calc_NPR_age(f_rate, N_recruits, N_age, sum_1plus, sum_mature)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===     
! Does numbers-per-recruit calculations (NPR). 
! Used to calculate the unexploited and exploited age-structure in the first year of the projections                 
!====== +++ === === +++ === === +++ === 
    real(kind = 8), intent(in) :: f_rate                ! Human caused mortality rate
    real(kind = 8), intent(in) :: N_recruits            ! Female recruits that are age zero (female calves) for NPR calcs and Initial_F()
    real(kind = 8), intent(out) :: N_age(0:age_x)       ! Numbers at age vector   
    real(kind = 8), intent(out) :: sum_1plus            ! Size of 1+ component 
    real(kind = 8), intent(out) :: sum_mature           ! Size of mature component 
    real(kind = 8) :: prop_N_age(0:age_x)               ! Proportion at age : sum(prop_NPR) = 1.0
    real(kind = 8) :: N_age_mat(0:age_x)                ! Numbers at age per that are mature
    real(kind = 8) :: N_age_immat(0:age_x)              ! Numbers at age per that are immature
    real(kind = 8) :: sum_N_age                         ! The sum of numbers across ages        
    integer(kind = 4) :: kk                             ! Counter for looping over ages
! Initialize to number of age zero recruits
    N_age(0) = N_recruits                               ! Fraction of female offspring per recruit, assuming 50:50 sex ratio at birth
! Loop over ages 
    do kk = 1, (age_x - 1)                                              ! NPR calculations from Age 1 to Age x-1
        N_age(kk) = N_age(kk - 1) * S_age(kk - 1)                       ! Subtract natural mortality 
        N_age(kk) = N_age(kk) * (1 - selectivity(kk - 1) * f_rate)      ! Subtract any human caused mortality 
    end do
! Calculate numbers in the plus-group        
    N_age(age_x) = N_age(age_x - 1) * S_age(age_x - 1)                  ! Natural mortality into plus-group
    N_age(age_x) = N_age(age_x) * (1 - selectivity(age_x - 1) * f_rate) ! Human caused mortality into plus-group         
    N_age(age_x) = N_age(age_x) / (1 - S_age(age_x) * (1 - selectivity(age_x) * f_rate)) ! Plus-group mortality
! Note these are vectorized operations below 
    N_age_mat = N_age * prop_mat_age          ! Mature females
    N_age_immat = N_age * (1 - prop_mat_age)  ! Numbers immature at age per recruit
    sum_1plus = sum(N_age(1 : age_x))         ! Keep track of one_plus component (females)
    sum_mature = sum(N_age_mat)               ! Sum mature females per female recruit
    sum_N_age = sum(N_age)                    ! Sum across ages in NPR vector
    prop_N_age = N_age / sum_N_age            ! Proportions at each age : sum(prop_N_age) = 1.0
    return
  end subroutine calc_NPR_age   
  
!====== +++ === === +++ === === +++ === 
!###### +++ ### ### +++ ### ### +++ ###     
  real(kind = 8) function Initial_F(f_init) result(objf_f_init)
!###### +++ ### ### +++ ### ### +++ ###       
!====== +++ === === +++ === === +++ === 
! Given numbers per-recruit at equilibrium - apply initial human caused mortality rate (F_init)
! Global Vars accessed:    
!        real(kind = 8) :: NPR_age                  ! Numbers per female recruit, see initialize_age_struc()      
!        real(kind = 8) :: NPR_age_tmp              ! Numbers at age conditioned on F_init
!        real(kind = 8) :: selectivity(age_x)       ! Selectivity at age
!        integer(kind = 4) :: a_r                   ! Age at recruitment (vulnerability) to human caused mortality
!        integer(kind = 4) :: age_x                 ! Plus group age (user input)
!        real(kind = 8) :: temp_mature          
!        real(kind = 8) :: temp_1plus
!        real(kind = 8) :: rec_init ! initial recruitment conditioned on Finit (in terms of females)
!        real(kind = 8) :: pred_rec
!        real(kind = 8) :: b_init
!        real(kind = 8) :: initial_oneplus
! b_sex_ratio : sex ratio at birth
!====== +++ === === +++ === === +++ === 
    real(kind = 8), intent(inout) :: f_init ! Initial human caused mortality rate 
    integer(kind = 4) :: ii, jj             ! Index counters for loops
    real(kind = 8) :: temp_mature           ! Size of mature component per female recruit given f_init = x_x
    real(kind = 8) :: temp_1plus            ! Size of 1+ component per female recruit given f_init = x_x
    real(kind = 8) :: foo                   ! calc_NPR_age returns sum_mature per recruit, which isn't necessary when finding f_init at the root of objf_f_init()
    real(kind = 8) :: rec_init              ! Initial recruitment conditioned on Finit (in terms of females)
    real(kind = 8) :: pred_rec              ! Predicted recruitment, used when finding f_init at the root of objf_f_init()
    real(kind = 8) :: initial_oneplus       ! Initial size of the 1+ component per female recruit (used in scaling)
    real(kind = 8) :: objf_f_init	! Objective function to be minimized in finding root at f_init
!====== +++ === === +++ === === +++ === ! Initialize
    temp_mature = 0.    
    temp_1plus = 0.     
    rec_init = 0.       
    pred_rec = 0.       
    initial_oneplus = 0.
    objf_f_init = -99.d0
!====== +++ === === +++ === === +++ === ! Initialize        
    Call calc_NPR_age(f_rate = f_init, N_recruits = b_sex_ratio, &  ! Do NPR calculations given f_init and N_recruits = NPR(age=0)
        N_age = NPR_age_tmp, sum_1plus = temp_1plus, &
        sum_mature = temp_mature)        
!====== +++ === === +++ === === +++ === ! Note: NPR_oneplus is global variable, based on NPR calculated at f_init = 0 (set by main program)
    rec_init = init_depl_i * NPR_oneplus / temp_1plus 	 ! Rescale to get the desired 1+ depletion conditioned on Finit 
    rec_init = rec_init / 2                                  ! Divide by two to get female recruits conditioned on Finit      
!====== +++ === === +++ === === +++ === ! Calculate actual numbers per recruit conditioned on f_init       
    NPR_age_tmp(0) = rec_init
    Call calc_NPR_age(f_rate = f_init, N_recruits = rec_init, &
        N_age = NPR_age_tmp, sum_1plus = initial_oneplus, &
        sum_mature = foo) 
    pred_rec = b_init * temp_mature 
    pred_rec = pred_rec - (init_depl_i * NPR_oneplus / initial_oneplus)
    pred_rec = b_sex_ratio * pred_rec                               ! In terms of female component (probably not necessary for root finding)
    objf_f_init = pred_rec * pred_rec 
    return
  end function Initial_F    
  
!====== +++ === === +++ === === +++ ===  
!###### +++ ### ### +++ ### ### +++ ###     
  subroutine assign_par_vectors(a_r, a_m, age_x, S_adult, S_juv, S_age, selectivity_age, prop_mat_age)
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===         
! Assign values to survival, proportion mature and selectivity at age vectors
!====== +++ === === +++ === === +++ === 
    integer(kind = 4), intent(in) :: a_r
    integer(kind = 4), intent(in) :: a_m
    integer(kind = 4), intent(in) :: age_x        
    real(kind = 8), intent(in) :: S_adult
    real(kind = 8), intent(in) :: S_juv
    real(kind = 8), intent(out) :: S_age(0:age_x)        
    real(kind = 8), intent(out) :: selectivity_age(0:age_x)
    real(kind = 8), intent(out) :: prop_mat_age(0:age_x)
! Note vectorized style for assignments below
    S_age(0 : (a_m - 1)) = S_juv            ! Assign juvenile survival rates to ages < a_m
    S_age(a_m : age_x) = S_adult            ! Assign adult survival rates to a_m <= ages = age_x                  
    prop_mat_age(0 : (a_m - 1)) = 0.0       ! Assign maturity-at-age vector (0 for immature ages, 1 for mature ages)
    prop_mat_age(a_m : age_x) = 1.0
    selectivity_age(0 : (a_r - 1)) = 0.0    ! Assign selectivity at age
    selectivity_age(a_r : age_x) = 1.0        
    return
  end subroutine assign_par_vectors  
  
!====== +++ === === +++ === === +++ ===       
!###### +++ ### ### +++ ### ### +++ ###     
    subroutine rescale_NPR(k_1plus_tmp, initial_oneplus_tmp, N_age_unscaled, N_age_scaled) ! , init_depl_i, initial_oneplus_i
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===     
!! Rescale numbers at age per-recruit to initial numbers at age vector, given initial depletion 
!====== +++ === === +++ === === +++ ===             
        real(kind = 8), intent(in) :: k_1plus_tmp                           ! Carrying capacity (in terms of age 1+) 
!        real(kind = 8) :: init_depl_i                                        ! Initial depletion (in terms of age 1+) for stock i
        real(kind = 8), intent(in) :: initial_oneplus_tmp                   ! Initial one plus for stock i        
        real(kind = 8) :: scale_pop                            ! Scalar to map numbers per recruit to initial numbers at age
        real(kind = 8), intent(in) :: N_age_unscaled(0:age_x) ! Unscaled numbers at age in initial year. Note: age_x is global variable
        real(kind = 8), intent(out) :: N_age_scaled(0:age_x) ! Unscaled numbers at age in initial year. Note: age_x is global variable              
        integer(kind = 4) :: aa             ! Index for age in numbers at age array 
!====== +++ === === +++ === === +++ ===                     
        print *, "Hello from rescale_NPR"
!        print *, "size(NPR_age)"
!        aa = size(NPR_age)
!        print *, aa
        print *, "NPR_age_tmp"
        print *, NPR_age_tmp
!        print *, "N_age"
!        print *, N_age
        
             
!!        print *, "first_yr, yr_max"
!!        print *, first_yr, yr_max
!
!! Time series of population components below (e.g. number of calving females each year)
!!  vector Nplus(first_yr , yr_max)          ! Vector of 1+ population size for all years
!!  vector Ntot(first_yr , yr_max)           ! Total (0+) population size each year, indexed to last year +1 - final pop size is at start of last year (2005?) + 1
!!  vector N_calf(first_yr , yr_max)         ! Vector of calf production for all years
!
!!  3darray NAll(1,4,First_yr,Last_yr+1,0,age_x);
!
!!  vector N_imm_yr_1(First_yr,Last_yr+1);	! Number of immature females by year
!!  vector N_recpt_yr_1(First_yr,Last_yr+1);	! Number of receptive females by year
!!  vector N_calvn_yr_1(First_yr,Last_yr+1);	! Number of calving females by year
!!  vector N_m_yr_1(First_yr,Last_yr+1);		! Numbers of males by year
!        
!!  vector N_dead(First_yr,Last_yr+1);	! Vector of natural dead each year        
!
!!  N_m_yr_1.initialize();
!!  N_calvn_yr_1.initialize();
!!  N_recpt_yr_1.initialize();
!!  N_imm_yr_1.initialize();
!!  
!       
        print *, "k_1plus_tmp: "
        print *, k_1plus_tmp
        print *, "initial_oneplus_tmp: "
        print *, initial_oneplus_tmp
        scale_pop = b_sex_ratio * k_1plus_tmp / initial_oneplus_tmp ! Assume 50:50 birth rate (so same re-scaling factor used for each sex)  
        N_age_scaled = N_age_unscaled * scale_pop	! This now in terms of females - so, just assign another vector equal to this one to initialize males
        print *, "scale_pop: "
        print *, scale_pop!
!        print *, "Nage: ", Nage ! DEBUGGING    
!!
!!   do aa = 0, age_x
!!   
!!   
!!   end do
!   
!!  for(aa=0;aa<=age_x;aa++) {
!!	  NAll(1,First_yr,aa)=Nage(aa)*(1-prop_mat_a(aa)); 	! immature numbers at age for females
!!	  NAll(2,First_yr,aa)=Nage(aa)*prop_mat_a(aa)*(1-b_init);  ! receptive numbers at age for females
!!	  NAll(3,First_yr,aa)=Nage(aa)*prop_mat_a(aa)*b_init; 	! calving numbers at age for females
!!	  NAll(4,First_yr,aa)=Nage(aa);   			! intial numbers at age for males  	
!!     N_imm_yr_1(First_yr)   += NAll(1,First_yr,aa);		! total numbers at stage
!!     N_recpt_yr_1(First_yr) += NAll(2,First_yr,aa);
!!     N_calvn_yr_1(First_yr) += NAll(3,First_yr,aa);
!!     N_m_yr_1(First_yr)     += NAll(4,First_yr,aa);
!!     Ntot(First_yr)         += NAll(1,First_yr,aa)+NAll(2,First_yr,aa)+NAll(3,First_yr,aa)+NAll(4,First_yr,aa);
!!  }
!
!        do ii = 1, n_stocks
!            depl_i_t(ii , first_yr) = init_depl(ii)
!        end do
!
!        print *, "depl_i_t(): ", depl_i_t(1:2, 0)
!
!!  N_calf(First_yr) = NAll(1,First_yr,0)*2;           // Initial number of calves (for output)
!!  Nplus(First_yr) = Ntot(First_yr)-N_calf(First_yr);
        
       return     
    end subroutine rescale_NPR
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ === 
    real(kind = 8) function assign_area_stock_prop(p_a1_s1, p_a2_s1, p_a2_s2, p_a3_s2, p_a4_s2) result(area_stock_prop)
! Assign percentage of each stock (columns) to areas (rows)
!====== +++ === === +++ === === +++ ===     
        real(kind = 8) :: p_a1_s1       ! Percentage of stock 1 in area 1
        real(kind = 8) :: p_a2_s1       ! Percentage of stock 1 in area 1        
        real(kind = 8) :: p_a2_s2       ! Percentage of stock 1 in area 1
        real(kind = 8) :: p_a3_s2       ! Percentage of stock 1 in area 1
        real(kind = 8) :: p_a4_s2       ! Percentage of stock 1 in area 1        
        real(kind = 8), dimension(4,2) :: area_stock_prop ! Hard-coded, could be improved
        area_stock_prop(1, 1) = p_a1_s1 ! Assign percentage of each stock (columns) to areas (rows)
        area_stock_prop(2, 1) = p_a2_s1
        area_stock_prop(3, 1) = 0.0d0
        area_stock_prop(4, 1) = 0.0d0
        area_stock_prop(1, 2) = 0.0d0
        area_stock_prop(2, 2) = p_a2_s2
        area_stock_prop(3, 2) = p_a3_s2
        area_stock_prop(4, 2) = p_a4_s2
        return
    end function assign_area_stock_prop
    
    real(kind = 8) function assign_transition_matrix(a_m, age_x, b_rate, S_age, prop_mat_age) result(transition_matrix)
        integer(kind = 4) :: a_m, age_x
        integer(kind = 4) :: ii
        real(kind = 8) :: b_rate
        real(kind = 8) :: S_age(0:age_x)
        real(kind = 8) :: prop_mat_age(0:age_x)
        real(kind = 8) :: transition_matrix(0:age_x, 0:age_x)
    
        transition_matrix(0, :) = b_rate * prop_mat_age            ! Assign birth rates to mature ages in first row of matrix
        
        do ii = 0, (a_m - 1)
            transition_matrix(ii+1, ii) = S_age(ii)   ! Assign juvenile survival rates
        end do
        transition_matrix(age_x, a_m) = S_age(age_x - 1)             ! Adult survival assumed for maturing animals transitioning into plus-group
        transition_matrix(age_x, age_x) = S_age(age_x)               ! Adult survival assumed for plus-group
          
        
    end function assign_transition_matrix
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ === 
    
    subroutine power_method ( n, a, y, it_max, tol, lambda, it_num )
!*****************************************************************************80
!
!! POWER_METHOD applies the power method for a real eigenvalue.
!
!  Discussion:
!
!    For a given NxN matrix A and an N vector Y, the power method produces
!    a series of estimates for LAMBDA, the largest eigenvalue, and Y,
!    the eigenvector corresponding to LAMBDA.
!
!    The iteration repeats the following steps
!
!      AY     = A * Y
!      LAMBDA = || AY ||
!      Y      = AY / LAMBDA
!
!    If the matrix A has a single real eigenvalue of maximum modulus,
!    then this iteration will generally produce a good estimate for that
!    eigenvalue and its corresponding eigenvector.
!
!    If there are multiple distinct eigenvalues of the same modulus,
!    perhaps two values of opposite sign, or complex eigenvalues, then
!    the situation is more complicated.
!
!    Separate issues:
!
!    * when estimating the value of LAMBDA, we use the Rayleigh quotient,
!    LAMBDA = ( y' * A * y ) / ( y' * y ).  Since we normalize Y, the
!    bottom of the fraction is 1.  Using this estimate allows us to
!    easily capture the sign of LAMDBA.  Using the eucldean norm
!    instead, for instance, would always give a positive value.
!
!    * If the dominant eigenvalue is negative, then the iteration
!    as given will produce eigenvector iterates that alternate in sign.
!
!    * It is worth knowing whether the successive eigenvector estimates
!    are tending to some value.  Since an eigenvector is really a direction,
!    we need to normalize the vectors, and we need to somehow treat both
!    a vector and its negative as holding the same information.  This
!    means that the proper way of measuring the difference between two
!    eigenvector estimates is to normalize them both, and then compute
!    the cosine between them as y1'y2, followed by the sine, which is
!    sqrt ( 1 - ( y1'y2)^2 ).  If this sine is small, the vectors y1 and y2
!    are "close" in the sense of direction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Input/output, real ( kind = 8 ) Y(N), the estimate for the eigenvector.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!    1 <= IT_MAX.
!
!    Input, real ( kind = 8 ) TOL, an error tolerance.
!
!    Output, real ( kind = 8 ) LAMBDA, the estimate for the eigenvalue.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
    implicit none

    integer ( kind = 4 ) :: n
    integer ( kind = 4) :: ii, jj ! Counters for debugging added by JRB
    real ( kind = 8 ) :: a(n,n)
    real ( kind = 8 ) :: ay(n)
    real ( kind = 8 ) :: cos_y1y2
    logical, parameter :: debug = .true. ! DEBUGGING -- turn this flag back to .false. when satisfied procedure is working robustly
    integer ( kind = 4 ) :: it_max
    integer ( kind = 4 ) :: it_num
    real ( kind = 8 ) :: lambda
    real ( kind = 8 ) :: lambda_old
    real ( kind = 8 ) :: sin_y1y2
    real ( kind = 8 ) :: tol
    real ( kind = 8 ) :: val_dif
    real ( kind = 8 ) :: y(n)
    real ( kind = 8 ) :: y_old(n)

    if ( debug ) then
      print *, "Hello from Power Method: "
      print *, "Transtion Matrix: "
      do jj = 1, n
        write (*,"(100f4.3)") (a(jj, ii), ii = 1, n)    ! The 100f... is a bit of a hack. Works if <= 100 columns to be printed
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     IT      Lambda          Delta-Lambda    Delta-Y'
      write ( *, '(a)' ) ' '
    end if
    !
    !  Force Y to be a vector of unit norm.
    !
    print *, "y in: "
    print *, y
    y(1:n) = y(1:n) / sqrt ( sum ( y(1:n)**2 ) )
    print *, "y norm: "
    print *, y
    
    it_num = 0

    y_old(1:n) = y(1:n)
    !
    !  Compute AY = A*Y.
    !
    ay(1:n) = matmul ( a(1:n,1:n), y(1:n) )
    !
    !  Estimate LAMBDA = (AY,Y)/(Y,Y).
    !
    lambda = dot_product ( y(1:n), ay(1:n) )
    !
    !  Force AY to have unit norm.
    !  Replace Y by AY.
    !
    y(1:n) = ay(1:n) / sqrt ( sum ( ay(1:n)**2 ) )
    !
    !  The sign of Y is optional.  If LAMBDA is probably negative,
    !  switch sign of new Y to match old one.
    !
    if ( lambda < 0.0D+00 ) then
      y(1:n) = - y(1:n)
    end if

    val_dif = 0.0D+00
    cos_y1y2 = dot_product ( y(1:n), y_old(1:n) )
    sin_y1y2 = sqrt ( ( 1.0D+00 - cos_y1y2 ) * ( 1.0D+00 + cos_y1y2 ) )

    if ( debug ) then
      write ( *, '(2x,i5,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
        it_num, lambda, val_dif, sin_y1y2
    end if
    !
    !  Now repeat these steps in an iteration.
    !
    do it_num = 1, it_max

      lambda_old = lambda
      y_old(1:n) = y(1:n)

      ay(1:n) = matmul ( a(1:n,1:n), y(1:n) )
      lambda = dot_product ( y(1:n), ay(1:n) )
      y(1:n) = ay(1:n) / sqrt ( sum ( ay(1:n)**2 ) )
      if ( lambda < 0.0D+00 ) then
        y(1:n) = - y(1:n)
      end if

      val_dif = abs ( lambda - lambda_old )
      cos_y1y2 = dot_product ( y(1:n), y_old(1:n) )
      sin_y1y2 = sqrt ( ( 1.0D+00 - cos_y1y2 ) * ( 1.0D+00 + cos_y1y2 ) )

      if ( debug ) then
        write ( *, '(2x,i5,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          it_num, lambda, val_dif, sin_y1y2
      end if

      if ( val_dif <= tol ) then
        exit
      end if

    end do

    y(1:n) = ay(1:n) / lambda

    return
    end subroutine power_method

END MODULE initialize_pop

