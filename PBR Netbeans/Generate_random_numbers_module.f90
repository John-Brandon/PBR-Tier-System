!=== === +++ === === +++ === === +++ === ===
! File:   Generate_random_numbers_module.f90
! Author: John R. Brandon
! eMail:  jbrandon at gmail
! Date :  Apr 2015
! Developed under OS:  Mac OS 10.9.5 (x86_64-apple-darwin10.8.0 (64-bit))
! Language : Fortran 90/95
! Originally compiled using: GCC GNU gfortran v 4.2 (free open source Fortran 90/95 compiler) -- https://gcc.gnu.org/
! IDE: Netbeans 8.0.2
!====== +++ === === +++ === === +++ === ===
! Purpose : Create the seed vector to initialize the random number generator (RNG).
!  This code takes values from an input file (e.g. input.par) and creates the seed vector accordingly.
!  For example, if the input file specifies a certain value for seed (iseed), that will be used. Otherwise,
!   the system clock will be used to create the seed vector. 
!  The RNG is initialized by calling the intrinsic random_seed() routine and passing it the seed vector.
!====== +++ === === +++ === === +++ === ===

MODULE Generate_random_numbers_module
    
    use Declare_variables_module ! contains variables read from input file re: how to seed RNG
    
    implicit none
    
    contains
!#######################################################         
    subroutine set_random_seed()
        implicit none ! this is not necessary to include in procedure, because 'implicit none' is declared at the module level
        integer, allocatable :: seed(:) ! This is a general routine, so leaving variable declarations here for portability 
        integer :: seed_n  ! Length of the seed vector       
        integer :: init_seed ! Value underlying creation of seed vector values  
         
        call random_seed(size = seed_n) ! Call the intrinsic random_seed subroutine and return length of seed vector (might be compiler dependent)
        allocate(seed(seed_n)) ! Dimension the seed vector accordingly (using allocate to dynamically dimension vector)
         
!        print *, "seed(seed_n): ", seed ! in gfortran 4.2 compiler on Mac OS X 10.9, seed is vector of length eight (zeros at this stage)
!        print *, "seed_n: ", seed_n
         
        if(cseed.eq."Y") then      ! User wants to specify the seed for the RNG, e.g. so results are reproducible
            init_seed = iseed      ! User defined seed from input file
        else if(cseed.eq."N") then ! User wants a random seed for RNG instead
            call system_clock(count = init_seed)   ! Seed is set using CPU time
        else
            stop "ERROR: Program stopped running, cseed in input file needs to = Y or N"
        end if     
         
!         print *, "init_seed: ", init_seed
         do ii = 1, seed_n
             seed(ii) = init_seed + 37*(ii-1) ! recommended function for setting seeds of RNG because 37 is prime number and RNG likes primes
!             print *, "ii / seed(ii): ", ii, seed(ii)
         enddo

         call random_seed(put = seed) ! pass the seed vector to initialize the random number generator        
        
    end subroutine set_random_seed

END MODULE Generate_random_numbers_module
