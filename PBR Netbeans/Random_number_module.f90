!     
! File:   random_numbers_module.f90
! Author: johnbrandon
!
! Created on April 27, 2015, 4:04 PM
!

!    implicit none
MODULE RNG
    
    contains
!
    real(kind=8) function r8_normal_01(seed)
!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
        implicit none

        real ( kind = 8 ) r1
        real ( kind = 8 ) r2
        real ( kind = 8 ) r8_normal_01
        real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
        real ( kind = 8 ) r8_uniform_01
        integer ( kind = 4 ) seed
        real ( kind = 8 ) x

        r1 = r8_uniform_01 ( seed )
        r2 = r8_uniform_01 ( seed )
        x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

        r8_normal_01 = x

        return
    end function r8_normal_01 

!*****************************************************************************80    
    function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
        implicit none

        integer ( kind = 4 ) k
        real ( kind = 8 ) r8_uniform_01
        integer ( kind = 4 ) seed
         
        k = seed / 127773

        print *, "Hi from r8_uniform_01"
        print *, "Seed: ", seed
        print *, "k: ", k
        print *, "seed / 127773: ", seed / 127773
         
        seed = 16807 * ( seed - k * 127773 ) - k * 2836
        
        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if
      !
      !  Although SEED can be represented exactly as a 32 bit integer,
      !  it generally cannot be represented exactly as a 32 bit real number!
      !
        r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

        print *, "Goodbye from r8_uniform_01"
        print *, "r8_uniform_01:", r8_uniform_01
        
        return
    end function r8_uniform_01
    
END MODULE RNG
