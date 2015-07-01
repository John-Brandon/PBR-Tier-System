module eigen_module
    !use initialize_pop !, only: assign_par_vectors, assign_transition_matrix   ! For re-assigning values during search for s_juv
    use Declare_variables_module
    implicit none
    contains
!====== +++ === === +++ === === +++ ===        
!   * subroutine eigen() *    
! Wrapper for calling LAPACK procedure DGEEV, which calculates the eigenvalues and left and right eigenvectors of a 
!  non-symmetric square matrix
!    
!   * print_eigenvalues() *
! Formatted printing of vector of eigenvalues returned from eigen()
!
!   * print_eigenvectors() *
! Formatted printing of eigenvectors returned from eigen()
!
!   * subroutine extract_lambda() *    
! Extract the dominant (largest absolute value) real (not complex) eigenvalue from a vector of eigenvalues    
!====== +++ === === +++ === === +++ === 

    ! *    
!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===      
    real (kind=8) function extract_lambda(wr, wi, n)
!====== +++ === === +++ === === +++ ===        
! Extract the dominant (largest absolute value) real (not complex) eigenvalue from a vector of eigenvalues
!====== +++ === === +++ === === +++ ===     
        integer(kind = 4), intent(in) :: n         ! Length of vector of eigenvalues
        integer(kind = 4) :: j         ! Counter
        integer(kind = 4) :: j_keep    ! Keep track of index with dominant real eigenvalue
        real(kind = 8) :: lambda       ! Dominant real eigenvalue 
        real(kind = 8) :: wr(n), wi(n) ! The real and imaginary components of the vector of eigenvalues
          
        j_keep = 0  ! Initialize
        lambda = 0.
        
!        print *, "Hello from extract_lambda()" ! DEBUGGING
        
        do j = 1, n
         if( wi( j ) == 0.0 ) THEN   ! Is the imaginary part of complex number greater than zero? If not, this is a real root.
            
!            print *, "Real root at index: ", j
            
            if (ABS(wr(j)) > ABS(lambda))  then 
!                print *, "Previous lambda: ", lambda
!                print *, "Proposed lambda wr(j): ", wr(j)
                lambda = wr(j)
!                print *, "New lambda: ", lambda
            end if
         end if
        end do
        
        extract_lambda = lambda
        return
    end function extract_lambda

!###### +++ ### ### +++ ### ### +++ ### 
!====== +++ === === +++ === === +++ ===       
    subroutine eigen(A, N, lambda)
!====== +++ === === +++ === === +++ ===           
! Wrapper for calling LAPACK procedure DGEEV, which calculates the eigenvalues and left and right eigenvectors of a 
!  non-symetric square matrix   
! Notes:
!  DGEEV changes the contents of the matrix A on return.        
!====== +++ === === +++ === === +++ ===             
        integer(kind = 4), intent(in) :: N
        integer(kind = 4) :: ii, jj         ! Counters
        integer(kind = 4) :: LDA, LDVL, LDVR
        integer(kind = 4), parameter :: LWMAX = 1000
        integer(kind = 4) :: INFO, LWORK

        real(kind = 8), intent(in) :: A( 1:N, 1:N ) ! Square, but not necessarily symmetric matrix 
        real(kind = 8) :: VL( 1:N, 1:N )            ! Matrix with left eigenvectors 
        real(kind = 8) :: VR( 1:N, 1:N )            ! Matrix with right eigenvectors
        real(kind = 8) :: WR( 1:N )                 ! Vector with real part of eigenvalue
        real(kind = 8) :: WI( 1:N )                 ! Vector with complex part of eigenvalue
        real(kind = 8) :: WORK( LWMAX )
        real (kind = 8), intent(out) :: lambda      ! Dominant real eigenvalue

!        EXTERNAL         DGEEV
!        EXTERNAL         PRINT_EIGENVALUES, PRINT_EIGENVECTORS

!*
!*     .. Intrinsic Functions ..
!        INTRINSIC        INT, MIN
!*
!*     .. Executable Statements ..
!        print *, "Hello from testing DGEEV: " 
!        print *, "Transition Matrix: "             ! DEBUGGING
!        do ii = 1, N
!            write (*,"(100f8.3)") (A(ii, jj), jj = 1, N)    ! The 100f... is a bit of a hack. Works if <= 100 columns to be printed
!        end do
      
!        WRITE(*,*)'DGEEV Example Program Results'
!*
!*     Query the optimal workspace.
!*
        LWORK = -1      ! This tells procedure that a workspace query is assumed; the procedure only calculates optimal size of WORK array
        LDVL = N
        LDVR = N
        LDA = N
        !
        CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, &
                   VR, LDVR, WORK, LWORK, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )  ! Set this for next call to DGEEV, which will return eigenvalues / vectors

!        write(*,*) "LWMAX: ", LWMAX                    ! DEBUGGING
!        write(*,*) "INT( WORK(1) ): ", INT( WORK(1) )            
!        write(*,*) "LWORK: ", LWORK
!*
!*     Solve eigenproblem.
!*      Note that matrix A will have been overwritten on return
        CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, &
                 VR, LDVR, WORK, LWORK, INFO )
!*
!*     Check for convergence.
!*
        IF( INFO.GT.0 ) THEN
           WRITE(*,*)'The DGEEV algorithm failed to compute eigenvalues.'
           STOP
        END IF

!     Extract dominant real eigenvalue 
        lambda = extract_lambda( WR, WI, N )
!
!     Print eigenvalues.
!
!        CALL print_eigenvalues( 'Eigenvalues', N, WR, WI )
!
!     Print left eigenvectors.
!
!        CALL PRINT_EIGENVECTORS( 'Left eigenvectors', N, WI, VL )
!!
!!     Print right eigenvectors.
!!
!        CALL PRINT_EIGENVECTORS( 'Right eigenvectors', N, WI, VR )
        
        return
      END subroutine eigen
!
!     End of DGEEV Example.
!
!  =============================================================================
!
!     Auxiliary routine: printing eigenvalues.
!
    subroutine print_eigenvalues( DESC, N, WR, WI )
        character*(*), intent(in) :: DESC
        integer(kind = 4), intent(in) :: N
        real(kind = 8), intent(in) :: WR( 1:N )
        real(kind = 8), intent(in) :: WI( 1:N )
        real(kind = 8), parameter :: zero = 0.0d0
        integer(kind = 4) :: jj
!
        WRITE(*,*)
        WRITE(*,*) DESC
        DO jj = 1, N
           IF( WI( jj ) .EQ. ZERO ) THEN   ! Is the imaginary part of complex number greater than zero? If not, is a real root.
              WRITE(*,9998,ADVANCE='NO') WR( jj ) 
           ELSE
              WRITE(*,9999,ADVANCE='NO') WR( jj ), WI( jj )
           END IF
        END DO
        WRITE(*,*)

9998    FORMAT( 11(:,1X,F6.2) )
9999    FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )

        RETURN
    END subroutine print_eigenvalues
!*
!*     Auxiliary routine: printing eigenvectors.
!*
    SUBROUTINE PRINT_EIGENVECTORS( DESC, N, WI, V )
        character*(*), intent(in) :: DESC
        integer(kind = 4), intent(in) :: N
        real(kind = 8), intent(in) :: WI( 1:N ), V( 1:N, 1:N )
        real(kind = 8), parameter :: ZERO = 0.0d0
        INTEGER          I, J
!*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, N
         J = 1
         DO WHILE( J.LE.N )
            IF( WI( J ).EQ.ZERO ) THEN
               WRITE(*,9998,ADVANCE='NO') V( I, J )
               J = J + 1
            ELSE
               WRITE(*,9999,ADVANCE='NO') V( I, J ), V( I, J+1 )
               WRITE(*,9999,ADVANCE='NO') V( I, J ), -V( I, J+1 )
               J = J + 2
            END IF
         END DO
         WRITE(*,*)
      END DO
!*
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END subroutine PRINT_EIGENVECTORS
      
      end module eigen_module      

