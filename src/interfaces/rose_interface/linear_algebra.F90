module linear_algebra

! Module to do linear algebra operations for real, complex and quaternion matrices

  use rose_global

  implicit none


  interface matrix_mul
      module procedure matrix_mul2
      module procedure matrix_mul3
      module procedure matrix_mul4
      module procedure matrix_mul5
  end interface

contains

   subroutine define_algebra_type(spatial_orbitals,restricted,rcq)
        logical, intent(in) :: restricted, spatial_orbitals
        integer, intent(out):: rcq

        write(luout,*) "Select algebra (rcq) based on the use of restricted and spatial orbitals:"
        write(luout,*) "  rcq = 1, 2 or 4  : real, complex or quaternion algebra, respectively"
        write(luout,*) "  restricted       : ",restricted
        write(luout,*) "  spatial_orbitals : ",spatial_orbitals
        if (spatial_orbitals) then
           rcq = 1
        else if (.not. spatial_orbitals .and. .not. restricted) then
           rcq = 2
        else if (.not. spatial_orbitals .and. restricted) then
           rcq = 4
        end if
        write(luout,'(A3,A,I2)') "   ","rcq              : ", rcq

   end subroutine define_algebra_type

   function multiplication(a,b) result(c)
        implicit none
        real(8) :: a(:),b(:)
        real(8) :: c(size(a))
        if (size(a) .ne. size(b)) stop 'a and b not same type, cannot be multiplied'
        if (size(a) .eq. 1) then
          c(1) = a(1)*b(1)
        else if (size(a) .eq. 2) then
          c(1) = a(1)*b(1) - a(2)*b(2)
          c(2) = a(1)*b(2) + a(2)*b(1)
        else if (size(a) .eq. 4) then
          c(1) = a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4)
          c(2) = a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3)
          c(3) = a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2)
          c(4) = a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1)
        endif
   end function multiplication

   function division(a,b) result(c)
        implicit none
        real(8) :: a(:),b(:)
        real(8) :: c(size(a)), norm
        if (size(a) .ne. size(b)) stop 'a and b not same type, cannot be divided'
        if (size(a) .eq. 1) then
          c(1) = a(1)/b(1)
        else if (size(a) .eq. 2) then
          c(1) = (a(1)*b(1) + a(2)*b(2))/norm2(b)**2
          c(2) = (a(2)*b(1) - a(1)*b(2))/norm2(b)**2
        else if (size(a) .eq. 4) then
          c(1) = (b(1)*a(1) + b(2)*a(2) + b(3)*a(3) + b(4)*a(4))/norm2(b)**2
          c(2) = (b(1)*a(2) - b(2)*a(1) + b(3)*a(4) - b(4)*a(3))/norm2(b)**2
          c(3) = (b(1)*a(3) - b(2)*a(4) - b(3)*a(1) + b(4)*a(2))/norm2(b)**2
          c(4) = (b(1)*a(4) + b(2)*a(3) - b(3)*a(2) - b(4)*a(1))/norm2(b)**2
        endif
   end function division

   subroutine transpose_conjg(m,conjgtransm)
        implicit none
        integer             :: i,j
        real(8), intent(in) :: m(:,:,:) ! use assumed shape array
        real(8), intent(out):: conjgtransm(size(m,2),size(m,1),size(m,3))

        conjgtransm = 0.D0
        do i = 1, size(m,1)
         do j = 1, size(m,2)
          if (size(m,3) .eq. 1) then
             conjgtransm(j,i,1) = m(i,j,1)
          else if (size(m,3) .eq. 2) then
             conjgtransm(j,i,1) =   m(i,j,1)
             conjgtransm(j,i,2) = - m(i,j,2)
          else if (size(m,3) .eq. 4) then
             ! q = a + ib + jc + kd, q* = a - ib - jc - kd
             conjgtransm(j,i,1) =   m(i,j,1)
             conjgtransm(j,i,2) = - m(i,j,2)
             conjgtransm(j,i,3) = - m(i,j,3)
             conjgtransm(j,i,4) = - m(i,j,4)
          else
             print*,"Error in the dimension defining the algebra."
             stop
          end if
         end do
        end do
   end subroutine transpose_conjg

   subroutine dot_prod(vec1,vec2,v)
        ! dot_product, i.e. < vec1 | vec2 >. The complex conjugate of vec1 should be considered !
        integer               :: rcq, n
        real(8), intent(in)   :: vec1(:,:), vec2(:,:) ! the second dimension denotes the algebra (real = 1, complex = 2, quaternion = 4)
        real(8), intent(out)  :: v(size(vec1,2)) ! dot_product can be real, complex or quaternion. The dimension of this array corresponds to the algebra. The return value is then an array and not a number !
        real(8), allocatable  :: v11(:),v12(:),v13(:),v14(:)
        real(8), allocatable  :: v21(:),v22(:),v23(:),v24(:)

        if (size(vec1,1) .ne. size(vec2,1)) then
           print*,"Dot product between two array of different size."
           stop
        end if
        if (size(vec1,2) .ne. size(vec2,2)) then
           print*,"Dot product between two array of different algebra."
           stop
        end if
        rcq = size(vec1,2)
        n   = size(vec1,1)

        v = 0.D0
        if (rcq .eq. 1) then
           allocate(v11(n))
           allocate(v21(n))
           v11  = vec1(:,1)
           v21  = vec2(:,1)
           v(1) = dot_product(v11,v21)
           deallocate(v11,v21)
        else if (rcq .eq. 2) then
           allocate(v11(n))
           allocate(v12(n))
           allocate(v21(n))
           allocate(v22(n))
           ! vec1(:,2) has a minus sign because we consider the conjugate (vec1 is in the bra)
           v11  =   vec1(:,1)
           v12  = - vec1(:,2)
           v21  =   vec2(:,1)
           v22  =   vec2(:,2)
           v(1) = dot_product(v11,v21) - dot_product(v12,v22)
           v(2) = dot_product(v11,v22) + dot_product(v12,v21) ! if vec1 = vec2 the imaginary part of the dot product is 0.
           deallocate(v11,v12,v21,v22)
        else if (rcq .eq. 4) then
           allocate(v11(n))
           allocate(v12(n))
           allocate(v13(n))
           allocate(v14(n))
           allocate(v21(n))
           allocate(v22(n))
           allocate(v23(n))
           allocate(v24(n))
           ! Use quaternion multiplication:
           !v1 | 1   i   j   k | --> v2
           !---|---------------|
           ! 1 | 1   i   j   k |
           ! i | i  −1   k  −j |
           ! j | j  −k  −1   i |
           ! k | k   j  −i  −1 |
           !--------------------
           ! vec1(:,2), vec1(:,3) and vec1(:,4) have a minus sign because we consider the conjugate (vec1 is in the bra)
           v11  =   vec1(:,1)
           v12  = - vec1(:,2)
           v13  = - vec1(:,3)
           v14  = - vec1(:,4)
           v21  =   vec2(:,1)
           v22  =   vec2(:,2)
           v23  =   vec2(:,3)
           v24  =   vec2(:,4)
           v(1) = dot_product(v11,v21) - dot_product(v12,v22) - dot_product(v13,v23) - dot_product(v14,v24)
           v(2) = dot_product(v11,v22) + dot_product(v12,v21) + dot_product(v13,v24) - dot_product(v14,v23) ! if v1 = v2 this part is 0
           v(3) = dot_product(v11,v23) - dot_product(v12,v24) + dot_product(v13,v21) + dot_product(v14,v22) ! if v1 = v2 this part is 0
           v(4) = dot_product(v11,v24) + dot_product(v12,v23) - dot_product(v13,v22) + dot_product(v14,v21) ! if v1 = v2 this part is 0
           deallocate(v11,v12,v13,v14,v21,v22,v23,v24)
        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if
   end subroutine dot_prod

   subroutine matrix_mul2(m1,m2,matrix)
        implicit none
        integer              :: i,j,k
        real(8), intent(in)  :: m1(:,:,:), m2(:,:,:) ! use assumed shape array
        real(8), intent(out) :: matrix(size(m1,1),size(m2,2),size(m2,3))

        if (size(m1,3) .ne. size(m2,3)) then
           write(luout,*)"matrix multiplication cannot operate between two different matrix types (real, complex or quaternion)."
           stop
        end if

        matrix = 0.D0
        do i = 1, size(m1,1)
         do j = 1, size(m2,2)
          do k = 1, size(m1,2)
           if (size(m1,3) .eq. 1) then
              matrix(i,j,1) = matrix(i,j,1) + m1(i,k,1)*m2(k,j,1)
           else if (size(m1,3) .eq. 2) then
              matrix(i,j,1) = matrix(i,j,1) + m1(i,k,1)*m2(k,j,1) - m1(i,k,2)*m2(k,j,2)
              matrix(i,j,2) = matrix(i,j,2) + m1(i,k,1)*m2(k,j,2) + m1(i,k,2)*m2(k,j,1)
           else if (size(m1,3).eq. 4) then
              ! Use quaternion multiplication:
              !m1 | 1   i   j   k | --> m2
              !---|---------------|
              ! 1 | 1   i   j   k |
              ! i | i  −1   k  −j |
              ! j | j  −k  −1   i |
              ! k | k   j  −i  −1 |
              !--------------------
              matrix(i,j,1) = matrix(i,j,1) + m1(i,k,1)*m2(k,j,1) - m1(i,k,2)*m2(k,j,2) - m1(i,k,3)*m2(k,j,3) - m1(i,k,4)*m2(k,j,4)
              matrix(i,j,2) = matrix(i,j,2) + m1(i,k,1)*m2(k,j,2) + m1(i,k,2)*m2(k,j,1) + m1(i,k,3)*m2(k,j,4) - m1(i,k,4)*m2(k,j,3)
              matrix(i,j,3) = matrix(i,j,3) + m1(i,k,1)*m2(k,j,3) - m1(i,k,2)*m2(k,j,4) + m1(i,k,3)*m2(k,j,1) + m1(i,k,4)*m2(k,j,2)
              matrix(i,j,4) = matrix(i,j,4) + m1(i,k,1)*m2(k,j,4) + m1(i,k,2)*m2(k,j,3) - m1(i,k,3)*m2(k,j,2) + m1(i,k,4)*m2(k,j,1)
           else
              print*,"Error in the dimension defining the algebra."
              stop
           end if
          end do
         end do
        end do
   end subroutine matrix_mul2

   subroutine matrix_mul3(m1,m2,m3,matrix)
        implicit none
        real(8), intent(in)  :: m1(:,:,:), m2(:,:,:), m3(:,:,:) ! use assumed shape array
        real(8), allocatable :: intermed(:,:,:)
        real(8), intent(out) :: matrix(size(m1,1),size(m3,2),size(m1,3))

        if (size(m1,3) .ne. size(m2,3) .or. size(m1,3) .ne. size(m3,3)) then
           write(luout,*)"matrix multiplication cannot operate between two different matrix types (real, complex or quaternion)."
           stop
        end if

        allocate(intermed(size(m1,1),size(m2,2),size(m1,3)))
        call matrix_mul(m1,m2,intermed)
        call matrix_mul(intermed,m3,matrix)
        deallocate(intermed)

   end subroutine matrix_mul3

   subroutine matrix_mul4(m1,m2,m3,m4,matrix)
        implicit none
        real(8), intent(in)  :: m1(:,:,:), m2(:,:,:), m3(:,:,:), m4(:,:,:) ! use assumed shape array
        real(8), allocatable :: intermed1(:,:,:), intermed2(:,:,:)
        real(8), intent(out) :: matrix(size(m1,1),size(m4,2),size(m1,3))

        if (size(m1,3) .ne. size(m2,3) .or. size(m1,3) .ne. size(m3,3) .or. size(m1,3) .ne. size(m4,3)) then
           write(luout,*)"matrix multiplication cannot operate between two different matrix types (real, complex or quaternion)."
           stop
        end if

        allocate(intermed1(size(m1,1),size(m2,2),size(m1,3)))
        allocate(intermed2(size(m3,1),size(m4,2),size(m1,3)))
        call matrix_mul(m1,m2,intermed1)
        call matrix_mul(m3,m4,intermed2)
        call matrix_mul(intermed1,intermed2,matrix)
        deallocate(intermed1,intermed2)

   end subroutine matrix_mul4

   subroutine matrix_mul5(m1,m2,m3,m4,m5,matrix)
        implicit none
        real(8), intent(in)  :: m1(:,:,:), m2(:,:,:), m3(:,:,:), m4(:,:,:), m5(:,:,:) ! use assumed shape array
        real(8), allocatable :: intermed1(:,:,:), intermed2(:,:,:), intermed3(:,:,:)
        real(8), intent(out) :: matrix(size(m1,1),size(m5,2),size(m1,3))

        if (size(m1,3) .ne. size(m2,3) .or. size(m1,3) .ne. size(m3,3) &
           & .or. size(m1,3) .ne. size(m4,3) .or. size(m1,3) .ne. size(m5,3)) then
           write(luout,*)"matrix multiplication cannot operate between two different matrix types (real, complex or quaternion)."
           stop
        end if

        allocate(intermed1(size(m1,1),size(m2,2),size(m1,3)))
        allocate(intermed2(size(m3,1),size(m4,2),size(m1,3)))
        allocate(intermed3(size(m1,1),size(m4,2),size(m1,3)))
        call matrix_mul(m1,m2,intermed1)
        call matrix_mul(m3,m4,intermed2)
        call matrix_mul(intermed1,intermed2,intermed3)
        call matrix_mul(intermed3,m5,matrix)
        deallocate(intermed1,intermed2,intermed3)

   end subroutine matrix_mul5

   subroutine solve_lineq(A,B,X)
        implicit none
        real(8), intent(in)    :: A(:,:,:),B(:,:,:)
        real(8), intent(inout) :: X(:,:,:)
        integer                :: rcq

        rcq = size(A,3)
        if ((rcq .eq. 1) .or. (rcq .eq. 2)) call lapack_solvelineq(A,B,X)
        if (rcq .eq. 4) call Gauss_Seidel_solvelineq(A,B,X)

   end subroutine solve_lineq

   subroutine lapack_solvelineq(A,B,X)
        ! solve AX = B. A and B are usually replaced by the LAPACK/BLAS subroutine, but it does not matter as they are only given as input.
        implicit none
        integer                :: nmo1,nmo2,rcq,ierr, i ,j, dnmo1, dnmo2
        real(8), intent(in)    :: A(:,:,:),B(:,:,:)
        real(8), intent(inout) :: X(:,:,:)
        real(8),    allocatable:: A1(:,:),X1(:,:)
        complex(8), allocatable:: A2(:,:),X2(:,:)
        integer, allocatable   :: ipiv(:)

        if ((size(A,3) .ne. size(B,3)) .or. (size(A,3) .ne. size(X,3)))&
          & stop 'solve_lineq: A, B or X work in different algebra.'
        if (size(A,1) .ne. size(B,1)) stop 'solve_lineq: A and B have different number of rows.'
        if (size(A,2) .ne. size(X,1)) stop 'solve_lineq: A and X dimensions not consistent'
        if (size(X,2) .ne. size(B,2)) stop 'solve_lineq: B and X dimensions not consistent'

        rcq  = size(A,3)
        nmo1 = size(A,1)
        nmo2 = size(B,2)

        if (rcq .eq. 1) then
           allocate(A1(nmo1,nmo1))
           allocate(X1(nmo1,nmo2))
           A1 = A(:,:,1)
           X1 = B(:,:,1)
           allocate(ipiv(nmo1))
           call dgesv( nmo1, nmo2, A1, nmo1, ipiv, X1, nmo1, ierr )
           deallocate(ipiv)
           if (ierr .gt. 0) then
              write(luout,*)'solve_lineq: The diagonal element of the triangular factor of A,'
              write(luout,*)'U(',ierr,',',ierr,') is zero, so that'
              write(luout,*)'A is singular; the solution could not be computed.'
              stop
           end if
           X(:,:,1) = X1
           deallocate(A1,X1)
        else if (rcq .eq. 2) then
           allocate(A2(nmo1,nmo1))
           allocate(X2(nmo1,nmo2))
           A2 = dcmplx(A(:,:,1),A(:,:,2))
           X2 = dcmplx(B(:,:,1),B(:,:,2))
           call zposv( 'U', nmo1, nmo2, A2, nmo1, X2, nmo1, ierr )
           if (ierr .gt. 0) then
              write(luout,*)'solve_lineq: The leading minor of order ',ierr,' is not positive'
              write(luout,*)'definite; the solution could not be computed.'
              stop
           end if
           X(:,:,1) = real(X2(:,:))
           X(:,:,2) = dimag(X2(:,:))
           deallocate(A2,X2)
        else if (rcq .eq. 4) then
           dnmo1 = 2*nmo1
           dnmo2 = 2*nmo2
           ! A quaternion can be written as a 2x2 complex matrix
           ! q = a + b*i + c*j + d*k
           ! --> | a + b*i , c + d*i |
           !     |-c + d*i , a - b*i |
           allocate (A2(dnmo1,dnmo1))
           allocate (X2(dnmo1,dnmo2))
           do j = 1, nmo1
              do i = 1, nmo1
                 A2(2*i-1,2*j-1) = dcmplx( A(i,j,1), A(i,j,2))
                 A2(2*i-1,2*j  ) = dcmplx( A(i,j,3), A(i,j,4))
                 A2(2*i  ,2*j-1) = dcmplx(-A(i,j,3), A(i,j,4))
                 A2(2*i  ,2*j  ) = dcmplx( A(i,j,1),-A(i,j,2))
              end do
           enddo
           do j = 1, nmo2
              do i = 1, nmo1
                 X2(2*i-1,2*j-1) = dcmplx( B(i,j,1), B(i,j,2))
                 X2(2*i-1,2*j  ) = dcmplx( B(i,j,3), B(i,j,4))
                 X2(2*i  ,2*j-1) = dcmplx(-B(i,j,3), B(i,j,4))
                 X2(2*i  ,2*j  ) = dcmplx( B(i,j,1),-B(i,j,2))
              end do
           end do
           call zposv( 'U', dnmo1, dnmo2, A2, dnmo1, X2, dnmo1, ierr )
           if (ierr .gt. 0) then
              write(luout,*)'solve_lineq: The leading minor of order ',ierr,' is not positive'
              write(luout,*)'definite; the solution could not be computed.'
              stop
           end if
           do j = 1, nmo2
              do i = 1, nmo1
                 ! Just like the diagonalization, it might be that we take the Kramers pair and not the wanted eigenvector because of degeneracy in quaternion algebra.
                 ! I should implement a quaternion version of this subroutine. Using Jacobi's method or Gauss-Seidel method ?
                 X(i,j,1) = real (X2(2*i-1,2*j-1),8)
                 X(i,j,2) = dimag(X2(2*i-1,2*j-1))
                 X(i,j,3) =-real (X2(2*i  ,2*j-1),8)
                 X(i,j,4) = dimag(X2(2*i  ,2*j-1))
              end do
           end do
           deallocate(A2,X2)
        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if

   end subroutine lapack_solvelineq

   subroutine Gauss_Seidel_solvelineq(A,B,X)
        ! solve AX = B
        ! Should converge faster than Jacobi's method, but less appropriate for parallel implementation.
        ! It converges for diagonally dominant matrix A (sufficient but not necessary condition).
        ! This is the case in our project as A is the S11 or S22 matrices. S11 is identity and S22 should be close to diagonally dominant.
        implicit none
        integer                :: nmo1,nmo2,rcq,i,j,k,conv,diag_dominant1,diag_dominant2,iter
        real(8), intent(in)    :: A(:,:,:),B(:,:,:)
        real(8), intent(inout) :: X(:,:,:)
        real(8), allocatable   :: summation(:), X1(:,:), X2(:,:), deviation_vector(:,:)
        real(8)                :: norm, norm_X1, norm_X2

        rcq  = size(A,3)
        nmo1 = size(A,1)
        nmo2 = size(B,2)

        ! Check if the matrix is diagonally dominant:
        diag_dominant1 = 1
        diag_dominant2 = 0
        do i=1,nmo1
         norm = 0
         do j=1,i-1
          norm = norm + euc_norm(A(i,j,:))
         enddo
         do j=i+1,nmo1
          norm = norm + euc_norm(A(i,j,:))
         enddo
         if (norm .gt. euc_norm(A(i,i,:))) diag_dominant1 = 0
         if (norm .lt. euc_norm(A(i,i,:))) diag_dominant2 = 1
        enddo
        if ((diag_dominant1 .eq. 0) .or. ((diag_dominant1 .eq. 1) .and. (diag_dominant2 .eq. 0))) &
          &  print*,"**Warning solvelineq: Matrix is not diagonally dominant. &
                       (a sufficient, but not necessary condition.)"
        
        ! Start the Gauss Seidel algorithm
        allocate(X1(nmo1,rcq),X2(nmo1,rcq),deviation_vector(nmo1,rcq))
        allocate(summation(rcq))
        ! First guess for X. Let's assume A is almost identity, so Xguess is almost B:
        X = B
        
        do i=1,nmo2
         conv = 0
         iter = 0
         do while (conv .ne. 1)
           X1 = X(:,i,:)
           do j=1,nmo1
            summation = 0.D0
            do k=1,nmo1
             if (k .eq. j) cycle
             summation(:) = summation(:) + multiplication(A(j,k,:),X(k,i,:))
            enddo
            X(j,i,:) = (1.D0/A(j,j,1))*(B(j,i,:) - summation(:))
           enddo
           X2 = X(:,i,:)
           deviation_vector(:,:) = X2(:,:) - X1(:,:)
           norm = norm_vector(deviation_vector)
           if (norm .le. 1.D-15) conv = 1
           iter = iter + 1
           if (iter .gt. 2048) then
            print*,"**WARNING: Gauss Seidel algorithm has not converged."
            exit
           endif
         enddo
        enddo

        deallocate(X1,X2,summation,deviation_vector)

   end subroutine Gauss_Seidel_solvelineq

   subroutine diag_matrix(A,eigval)
        implicit none
        real(8), intent(inout)       :: A(:,:,:)
        real(8), intent(  out)       :: eigval(:)
        integer                      :: rcq

        rcq = size(A,3)
        if ((rcq .eq. 1) .or. (rcq .eq. 2)) call lapack_diagmatrix(A,eigval)
        if (rcq .eq. 4) call jacobi_diagmatrix(A,eigval)

   end subroutine diag_matrix

   subroutine lapack_diagmatrix(A,eigval)
        implicit none
        real(8), intent(inout)       :: A(:,:,:)
        real(8), intent(  out)       :: eigval(:)
        integer                      :: n, info, lwork, liwork, lrwork, i, j, dn
        integer, allocatable         :: iwork(:)
        complex(8), allocatable      :: work(:)
        real(8), allocatable         :: rwork(:)
        real(8), allocatable         :: Dreal(:,:)
        complex(8), allocatable      :: Dcomp(:,:)
        real(8), allocatable         :: eigval2(:)

        n = size(A,1)
        if (size(A,3) .eq. 1) then
           allocate(Dreal(n,n))
           Dreal = A(:,:,1)
!
!          Query the optimal workspace and allocate the work array.
!
           lrwork = -1
           liwork = -1
           allocate(rwork(1))
           allocate(iwork(1))
           call dsyevd( 'V', 'L', n, Dreal, n, eigval, rwork, lrwork, &
                       iwork, liwork, info )
           lrwork = int(rwork(1))
           liwork = iwork(1)
           deallocate(rwork)
           deallocate(iwork)
           allocate(rwork(lrwork))
           allocate(iwork(liwork))
!
!          Solve eigenproblem.
!
           call dsyevd( 'V', 'L', n, Dreal, n, eigval, rwork, lrwork, &
                       iwork, liwork, info )
!
!          Check for convergence.
!
           if( info.gt.0 ) then
              write(luout,*)'The algorithm failed to compute eigenvalues.'
              stop
           end if

           A(:,:,1) = Dreal

           deallocate(Dreal)
           deallocate(rwork,iwork)

        else if (size(A,3) .eq. 2) then
           allocate(Dcomp(n,n))
           Dcomp = dcmplx(A(:,:,1),A(:,:,2))
!
!          Query the optimal workspace and allocate the work arrays.
!
           lwork = -1
           liwork = -1
           lrwork = -1
           allocate(work(1))
           allocate(iwork(1))
           allocate(rwork(1))
           call zheevd( 'V', 'L', n, Dcomp, n, eigval, work, lwork, rwork, &
                       lrwork, iwork, liwork, info )
           lwork = int(work(1))
           lrwork = int(rwork(1))
           liwork = iwork(1)
           deallocate(work)
           deallocate(rwork)
           deallocate(iwork)
           allocate(work(lwork))
           allocate(rwork(lrwork))
           allocate(iwork(liwork))
!
!          Solve eigenproblem.
!
           call zheevd( 'V', 'L', n, Dcomp, n, eigval, work, lwork, rwork, &
                       lrwork, iwork, liwork, info )
!
!          Check for convergence.
!
           if( info.gt.0 ) then
              write(luout,*)'The algorithm failed to compute eigenvalues.'
              stop
           end if

           A(:,:,1) = real(Dcomp(:,:),8)
           A(:,:,2) = dimag(Dcomp(:,:))

           deallocate(Dcomp)
           deallocate(work,rwork,iwork)

        else if (size(A,3) .eq. 4) then
           dn = 2*n
           ! A quaternion can be written as a 2x2 complex matrix
           ! q = a + b*i + c*j + d*k
           ! --> | a + b*i , c + d*i |
           !     |-c + d*i , a - b*i |
           allocate(Dcomp(dn,dn))
           do j = 1, n
              do i = 1, n
                 Dcomp(2*i-1,2*j-1) = dcmplx( A(i,j,1), A(i,j,2))
                 Dcomp(2*i-1,2*j  ) = dcmplx( A(i,j,3), A(i,j,4))
                 Dcomp(2*i  ,2*j-1) = dcmplx(-A(i,j,3), A(i,j,4))
                 Dcomp(2*i  ,2*j  ) = dcmplx( A(i,j,1),-A(i,j,2))
              enddo
           enddo

           allocate(eigval2(dn))
!
!          Query the optimal workspace and allocate work arrays.
!
           lwork = -1
           liwork = -1
           lrwork = -1
           allocate(work(1))
           allocate(iwork(1))
           allocate(rwork(1))
           call zheevd( 'V', 'L', dn, Dcomp, dn, eigval2, work, lwork, rwork, &
                       lrwork, iwork, liwork, info )
           lwork = int(work(1))
           lrwork = int(rwork(1))
           liwork = iwork(1)
           deallocate(work)
           deallocate(rwork)
           deallocate(iwork)
           allocate(work(lwork))
           allocate(rwork(lrwork))
           allocate(iwork(liwork))
!
!          Solve eigenproblem.
!
           call zheevd( 'V', 'L', dn, Dcomp, dn, eigval2, work, lwork, rwork, &
                       lrwork, iwork, liwork, info )
!
!          Check for convergence.
!
           if( info.gt.0 ) then
              write(luout,*)'The algorithm failed to compute eigenvalues.'
              stop
           end if

           ! Eigenvalues of the dn x dn matrix should come by pair for quaternions
           ! check:
           do i = 1, n
              if (abs(eigval2(2*i-1) - eigval2(2*i)) .gt. 1D-12) then
                 print*,"** Warning : Eigenvalues of quaternion matrix doesn't come by pair."
              end if
              eigval(i) = eigval2(2*i-1)
           end do

           do j = 1, n
              do i = 1, n
                 ! Assume that the correct eigenvectors are the odd ones. This is wrong, as you might get Kramers pair or not, which have degenerate eigenvalues.
                 ! You should use jacobi_diagmatrix subroutine for the quaternion diagonalization, which avoids the transformation to complex algebra and back.
                 A(i,j,1) = real (Dcomp(2*i-1,2*j-1),8)
                 A(i,j,2) = dimag(Dcomp(2*i-1,2*j-1))
                 A(i,j,3) =-real (Dcomp(2*i  ,2*j-1),8)
                 A(i,j,4) = dimag(Dcomp(2*i  ,2*j-1))
              end do
           end do

           deallocate(Dcomp,eigval2)
           deallocate(work,rwork,iwork)

        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if

   end subroutine lapack_diagmatrix

   subroutine jacobi_diagmatrix(matrix,eigenvalues)
        ! based on http://fourier.eng.hmc.edu/e176/lectures/ch1/node1.html
        ! return matrix as eigenvectors
        implicit none
        real(8), intent(inout) :: matrix(:,:,:)
        real(8), intent(inout) :: eigenvalues(:)
        integer                :: i,j,k,l,m,n,rcq
        real(8)                :: c,s,t,w
        real(8)                :: pivot_norm,matrix_ii,matrix_jj
        real(8), allocatable   :: eigenvectors(:,:,:),off_diag(:)
        real(8), allocatable   :: pivot_phase(:)
        real(8), allocatable   :: matrix_ik(:),matrix_jk(:)
        real(8), allocatable   :: eigenvectors_ki(:),eigenvectors_kj(:)
        integer, allocatable   :: ind(:)
        
        ! Initialization
        if (size(matrix,1) .ne. size(matrix,2)) stop "jacobi diagonalization: matrix has to be square"
        if (size(matrix,1) .ne. size(eigenvalues)) stop "jacobi diagonalization: wrong dimension in matrix or eigenvalues"
        n = size(matrix,1)
        rcq = size(matrix,3)
        allocate(ind(n)) ! denote the maxind of each row from 1 to n
        allocate(eigenvectors(n,n,rcq))
        allocate(pivot_phase(rcq))
        allocate(matrix_ik(rcq),matrix_jk(rcq))
        allocate(eigenvectors_ki(rcq),eigenvectors_kj(rcq))
        allocate(off_diag(n))
        eigenvectors = 0.D0
        forall (i=1:n) eigenvectors(i,i,1) = 1.D0

        ! find the off diagonal element of the greatest absolute value, the pivot
        pivot_norm = 1.D0
        do while (pivot_norm .gt. 1.D-15)
           ! this is very inefficient, doing ind all the time, also for rows that are not changed.
           ! also have the rcq index at the end instead of at the beginning, creating strides in memory access
           do i = 1, n-1
              ind(i) = maxind(matrix(i,:,:),i)
              off_diag(i) = euc_norm(matrix(i,ind(i),:))
           enddo
           m = 1 ! First row. 
           do i = 2, n-1 ! loop on rows (the last one does not contain any off-diagonal elements, obviously)
              if (off_diag(i) .gt. off_diag(m)) m = i
           enddo
           i = m
           j = ind(m)
           pivot_norm = euc_norm(matrix(i,j,:))
           if (pivot_norm .lt. 1.D-15) exit
           pivot_phase(:) = matrix(i,j,:)/pivot_norm
           ! rotate matrix into the real plane by scaling the j^th basis vector
           do k = 1, n
              if (k.eq.j) cycle
              matrix(k,j,:) = division(matrix(k,j,:),pivot_phase)
              matrix(j,k,:) = multiplication(pivot_phase,matrix(j,k,:))
           end do

           ! determine rotation in the real plane (c=cosine, s=sine, t=tangent)
           matrix_ii = matrix(i,i,1)
           matrix_jj = matrix(j,j,1)
           w = (matrix_jj - matrix_ii)/(2.D0*pivot_norm)
           if (w .le. 0) t = -w + dsqrt(w**2 + 1)
           if (w .gt. 0) t = -w - dsqrt(w**2 + 1)
           c = 1.D0/dsqrt(1 + t**2)
           s = t/dsqrt(1 + t**2)

           ! update all elements in the k-th and l-th rows and columns:
           matrix(i,i,1) = matrix_ii - t*pivot_norm
           matrix(j,j,1) = matrix_jj + t*pivot_norm
           matrix(i,j,:) = 0.D0
           matrix(j,i,:) = 0.D0
           do k=1,n
              if ((k.ne.i) .and. (k.ne.j)) then
                matrix_ik(:) = matrix(i,k,:)
                matrix_jk(:) = matrix(j,k,:)
                matrix(i,k,:) = c*matrix_ik(:) - s*matrix_jk(:)
                matrix(j,k,:) = s*matrix_ik(:) + c*matrix_jk(:)
                matrix(k,i,1) = matrix(i,k,1)
                matrix(k,j,1) = matrix(j,k,1)
                matrix(k,i,2:rcq) = -matrix(i,k,2:rcq)
                matrix(k,j,2:rcq) = -matrix(j,k,2:rcq)
              endif
           end do
           ! rotate eigenvectors after first scaling them with the phase factor
           do k=1,n
              eigenvectors(k,j,:) = division(eigenvectors(k,j,:),pivot_phase)
           end do
           do k=1,n
              eigenvectors_ki(:) = eigenvectors(k,i,:)
              eigenvectors_kj(:) = eigenvectors(k,j,:)
              eigenvectors(k,i,:) = c*eigenvectors_ki(:) - s*eigenvectors_kj(:)
              eigenvectors(k,j,:) = s*eigenvectors_ki(:) + c*eigenvectors_kj(:)
           enddo
        enddo

        deallocate(ind,off_diag,pivot_phase,matrix_ik,matrix_jk)
        deallocate(eigenvectors_ki,eigenvectors_kj)

        ! store in eigenvalues:
        forall (i=1:n) eigenvalues(i) = matrix(i,i,1)

        call sort_by_ascending_eigenvalues(eigenvectors,eigenvalues,.true.)
        matrix = eigenvectors
        deallocate(eigenvectors)
          
   end subroutine jacobi_diagmatrix

   integer function maxind(vec,k)  ! index of largest off-diagonal element in row k (largest element in vec)
     implicit none
     real(8), intent(in) :: vec(:,:) ! 2nd dimension is rcq.
     integer, intent(in) :: k
     integer             :: i

     maxind = k + 1 ! working with the upper triangle (no need to read elements below k+1)
     do i = k + 2, size(vec,1)
       if (euc_norm(vec(i,:)) .gt. euc_norm(vec(maxind,:))) maxind = i
     enddo

   end function maxind

   real(8) function euc_norm(vec) ! euclidian (L2) norm
     implicit none
     real(8), intent(in) :: vec(:)
             
     if (size(vec) .eq. 1) then
         euc_norm = abs(vec(1))
     else
         euc_norm = norm2(vec)
     end if

   end function euc_norm

   subroutine SVD(A,eigval,U,VT)
        implicit none
        real(8), intent(inout)       :: A(:,:,:)
        real(8), intent(  out)       :: eigval(:)
        real(8), intent(  out)       :: U(:,:,:)
        real(8), intent(  out)       :: VT(:,:,:)
        integer                      :: rcq

        rcq = size(A,3)
        if ((rcq .eq. 1) .or. (rcq .eq. 2)) call lapack_SVD(A,eigval,U,VT)
        if (rcq .eq. 4) call jacobi_SVD(A,eigval,U,VT)

   end subroutine SVD

   subroutine lapack_SVD(A,eigval,U,VT)
        implicit none
        real(8), intent(in   )       :: A(:,:,:) ! assumed shape array. A is not changed in the output.
        real(8), intent(  out)       :: eigval(:)
        real(8), intent(  out)       :: U(:,:,:)
        real(8), intent(  out)       :: VT(:,:,:)
        integer,    parameter        :: lwmax=100000
        integer                      :: m, n, info, lwork, i, j, dm, dn
        complex(8), dimension(lwmax) :: work
        real(8), allocatable         :: rwork(:)
        real(8), allocatable         :: Dreal(:,:), Unitary(:,:), VnitaryT(:,:)
        complex(8), allocatable      :: Dcomp(:,:), CUnitary(:,:), CVnitaryT(:,:)
        real(8), allocatable         :: eigval2(:)

        m = size(A,1)
        n = size(A,2)

        if (size(A,3) .eq. 1) then
           allocate(Dreal(m,n))
           Dreal = A(:,:,1)
!
!          Query the optimal workspace.
!
           lwork = -1
           allocate(Unitary(m,m))
           allocate(VnitaryT(n,n))

           call dgesvd('All','All',m,n,Dreal,m,eigval,Unitary,m,VnitaryT,n,&
                       work,lwork,info)

           lwork = min(lwmax,int(work(1)))
!
!          Compute SVD
!
           call dgesvd('All','All',m,n,Dreal,m,eigval,Unitary,m,VnitaryT,n,&
                       work,lwork,info)
!
!          Check for convergence.
!
           if( info.gt.0 ) then
              write(luout,*)'The algorithm computing SVD failed to converge. Info:',info
              write(luout,*)'Workspace lwork:',lwork
              stop
           end if

           U(:,:,1)  = Unitary
           VT(:,:,1) = VnitaryT

           deallocate(Dreal,Unitary,VnitaryT)

        else if (size(A,3) .eq. 2) then
           allocate(Dcomp(m,n))
           Dcomp = dcmplx(A(:,:,1),A(:,:,2))
!
!          Query the optimal workspace.
!
           lwork = -1
           allocate(rwork(5*min(m,n)))
           allocate(CUnitary(m,m))
           allocate(CVnitaryT(n,n))

           call zgesvd('All','All',m,n,Dcomp,m,eigval,CUnitary,m,CVnitaryT,n,&
                       work,lwork,rwork,info)

           lwork = min(lwmax,int(work(1)))
!
!          Compute SVD
!
           call zgesvd('All','All',m,n,Dcomp,m,eigval,CUnitary,m,CVnitaryT,n,&
                       work,lwork,rwork,info)
!
!          Check for convergence.
!
           if( info.gt.0 ) then
              write(luout,*)'The algorithm computing SVD failed to converge. Info:',info
              write(luout,*)'Workspace lwork:',lwork
              stop
           end if

           U(:,:,1)  = real(CUnitary(:,:),8)
           U(:,:,2)  = dimag(CUnitary(:,:))
           VT(:,:,1) = real(CVnitaryT(:,:),8)
           VT(:,:,2) = dimag(CVnitaryT(:,:))

           deallocate(Dcomp,CUnitary,CVnitaryT,rwork)

        else if (size(A,3) .eq. 4) then
           dm = 2*m
           dn = 2*n
           ! A quaternion can be written as a 2x2 complex matrix
           ! q = a + b*i + c*j + d*k
           ! --> | a + b*i , c + d*i |
           !     |-c + d*i , a - b*i |
           allocate(Dcomp(dm,dn))
           do j = 1, n
              do i = 1, m
                 Dcomp(2*i-1,2*j-1) = dcmplx( A(i,j,1), A(i,j,2))
                 Dcomp(2*i-1,2*j  ) = dcmplx( A(i,j,3), A(i,j,4))
                 Dcomp(2*i  ,2*j-1) = dcmplx(-A(i,j,3), A(i,j,4))
                 Dcomp(2*i  ,2*j  ) = dcmplx( A(i,j,1),-A(i,j,2))
              enddo
           enddo
!
!          Query the optimal workspace.
!
           lwork = -1
           allocate(rwork(5*min(dm,dn)))
           allocate(eigval2(min(dm,dn)))
           allocate(CUnitary(dm,dm))
           allocate(CVnitaryT(dn,dn))

           call zgesvd('All','All',dm,dn,Dcomp,dm,eigval2,CUnitary,dm,CVnitaryT,dn,&
                       work,lwork,rwork,info)

           lwork = min(lwmax,int(work(1)))
!
!          Compute SVD
!
           call zgesvd('All','All',dm,dn,Dcomp,dm,eigval2,CUnitary,dm,CVnitaryT,dn,&
                       work,lwork,rwork,info)
!
!          Check for convergence.
!
           if( info.gt.0 ) then
              write(luout,*)'The algorithm computing SVD failed to converge. Info:',info
              write(luout,*)'Workspace lwork:',lwork
              stop
           end if

           ! Eigenvalues of the diagonal matrix should come by pair for quaternions
           ! check:
           do i = 1, min(m,n)
              if (abs(eigval2(2*i-1) - eigval2(2*i)) .gt. 1D-12) then
                 print*,"** Warning : Eigenvalues of quaternion matrix do not come in pairs."
              end if
              eigval(i) = eigval2(2*i-1)
           end do

           do j = 1, m
              do i = 1, m
                 ! Just like the diagonalization, it might be that we take the Kramers pair and not the wanted eigenvector because of degeneracy in quaternion algebra.
                 ! I should implement a quaternion version of this subroutine.
                 U(i,j,1)  = real (CUnitary(2*i-1,2*j-1),8)
                 U(i,j,2)  = dimag(CUnitary(2*i-1,2*j-1))
                 U(i,j,3)  =-real (CUnitary(2*i  ,2*j-1),8)
                 U(i,j,4)  = dimag(CUnitary(2*i  ,2*j-1))
              end do
           end do
           do j = 1, n
              do i = 1, n
                 ! This is the transpose, so the eigenvector is on the row and not column.
                 ! Just like the diagonalization, it might be that we take the Kramers pair and not the wanted eigenvector because of degeneracy in quaternion algebra.
                 ! I should implement a quaternion version of this subroutine.
                 VT(i,j,1) = real (CVnitaryT(2*i-1,2*j-1),8)
                 VT(i,j,2) = dimag(CVnitaryT(2*i-1,2*j-1))
                 VT(i,j,3) =-real (CVnitaryT(2*i-1,2*j  ),8)
                 VT(i,j,4) = dimag(CVnitaryT(2*i-1,2*j  ))
              end do
           end do
           deallocate(Dcomp,eigval2,CUnitary,CVnitaryT)

        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if

   end subroutine lapack_SVD

   subroutine jacobi_SVD(A,eigval,U,VT)

        implicit none
        real(8), intent(in   )       :: A(:,:,:) ! assumed shape array. A is not changed in the output.
        real(8), intent(  out)       :: eigval(:)
        real(8), intent(  out)       :: U(:,:,:)
        real(8), intent(  out)       :: VT(:,:,:)
        real(8), allocatable         :: AT(:,:,:),V(:,:,:),B(:,:,:),AV(:,:,:), ATU(:,:,:)
        integer                      :: i,j,k,n,m,nval,rcq

        m   = size(A,1)
        n   = size(A,2)
        rcq = size(A,3)
        nval= min(m,n) ! number of singular values
        if (nval .ne. size(eigval)) stop 'SVD: error in dimension of eigval'

        allocate(AT(n,m,rcq))
        allocate(V(n,n,rcq))
        call transpose_conjg(A,AT)

        if (n.ge.m) then
           allocate(B(m,m,rcq))
           allocate(ATU(n,m,rcq))
           call matrix_mul(A,AT,B)
           call diag_matrix(B,eigval) ! return eigenvectors of B which are equal to U
           eigval = dsqrt(eigval)
           call sort_by_descending_eigenvalues(B,eigval,.true.)
           U = B
           call matrix_mul(AT,U,ATU)
           V = 0.D0
           do i=1,m ! just not considering the m:n eigenvectors...
             if (eigval(i) .ge. 1D-15) V(:,i,:) = ATU(:,i,:)/eigval(i)
             if (eigval(i) .lt. 1D-15) V(:,i,:) = ATU(:,i,:)
           enddo        
           deallocate(B,ATU)
        else
           allocate(B(n,n,rcq))
           allocate(AV(m,n,rcq))
           call matrix_mul(AT,A,B)
           call diag_matrix(B,eigval)
           eigval = dsqrt(eigval)
           call sort_by_descending_eigenvalues(B,eigval,.true.)
           V = B
           call matrix_mul(A,V,AV)
           U = 0.D0
           do i=1,n ! just not considering the n:m eigenvectors...
             if (eigval(i) .ge. 1D-15) U(:,i,:) = AV(:,i,:)/eigval(i)
             if (eigval(i) .lt. 1D-15) U(:,i,:) = AV(:,i,:)
           enddo
           deallocate(B,AV)
        end if
        call transpose_conjg(V,VT)
        deallocate(AT,V)

   end subroutine jacobi_SVD

   subroutine diff_norm2_rcq(m1,m2,norm)
        real(8), intent(in) :: m1(:,:,:), m2(:,:,:)
        real(8)             :: norm

        if (size(m1,3) .ne. size(m2,3)) then
           write(luout,*)"Norm between two different matrix types (real, complex or quaternion)."
           stop
        end if

        if (size(m1,3) .eq. 1) then
           norm = norm2(m1(:,:,1) - m2(:,:,1))
        else if (size(m1,3) .eq. 2) then
           norm = norm2([norm2(m1(:,:,1) - m2(:,:,1)), &
                  norm2(m1(:,:,2) - m2(:,:,2))])
        else if (size(m1,3) .eq. 4) then
           norm = norm2([norm2(m1(:,:,1) - m2(:,:,1)), &
                  norm2(m1(:,:,2) - m2(:,:,2)), &
                  norm2(m1(:,:,3) - m2(:,:,3)), &
                  norm2(m1(:,:,4) - m2(:,:,4))])
        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if
   end subroutine diff_norm2_rcq

   real(8) function norm_matrix(m)
        ! Return norm of a matrix
        real(8), intent(in) :: m(:,:,:)

        if (size(m,3) .eq. 1) then
           norm_matrix = norm2(m(:,:,1))
        else if (size(m,3) .eq. 2) then
           norm_matrix = norm2([norm2(m(:,:,1)), &
                  norm2(m(:,:,2))])
        else if (size(m,3) .eq. 4) then
           norm_matrix = norm2([norm2(m(:,:,1)), &
                  norm2(m(:,:,2)), &
                  norm2(m(:,:,3)), &
                  norm2(m(:,:,4))])
        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if
   end function norm_matrix

   real(8) function norm_vector(v)
        ! Return norm of a vector
        real(8), intent(in) :: v(:,:)

        if (size(v,2) .eq. 1) then
           norm_vector = norm2(v(:,1))
        else if (size(v,2) .eq. 2) then
           norm_vector = norm2([norm2(v(:,1)), &
                  norm2(v(:,2))])
        else if (size(v,2) .eq. 4) then
           norm_vector = norm2([norm2(v(:,1)), &
                  norm2(v(:,2)), &
                  norm2(v(:,3)), &
                  norm2(v(:,4))])
        else
           print*,"Error in the dimension defining the algebra."
           stop
        end if
   end function norm_vector

   subroutine SymOrth(C,S)
        !Symmetrically orthogonalize orbitals C with respect to overlap matrix S (i.e., such that Cnew^T S Cnew == id).

        implicit none
        real(8), parameter           :: ThrDel=1.D-12, Thrsh=1.D-12
        integer                      :: i,j, iDel, rcq, nmo1, nmo2
        real(8),    intent(in)       :: S(:,:,:) ! assumed shape array
        real(8),    intent(inout)    :: C(:,:,:) ! assumed shape array
        real(8),    allocatable      :: Cnew(:,:,:)
        real(8),    allocatable      :: transconjgC(:,:,:)
        real(8),    allocatable      :: transconjgmatrix(:,:,:)
        real(8),    allocatable      :: matrix(:,:,:),matrix2(:,:,:), Id(:,:,:)
        real(8),    allocatable      :: ew(:)
        real(8),    allocatable      :: matrix1(:,:),eigenvalues(:),eigenvectors(:,:)
        real(8)                      :: norm

        if (size(C,3) .ne. size(S,3)) then
           write(luout,*)"SymOrth: two different matrix types (real, complex or quaternion)."
           stop
        end if
        if (size(C,1) .ne. size(S,2)) then
           write(luout,*)"SymOrth: wrong dimensions for C and S, Crow .neq. Scol."
           stop
        end if
        if (size(S,1) .ne. size(S,2)) then
           write(luout,*)"SymOrth: wrong dimensions for S, should be square matrix."
           stop
        end if

        rcq  = size(C,3)
        nmo1 = size(C,1)
        nmo2 = size(C,2)

        ! Compute C^\dagger * S * C
        allocate(matrix(nmo2,nmo2,rcq))
        allocate(ew(nmo2))
        allocate(transconjgC(nmo2,nmo1,rcq))
        call transpose_conjg(C,transconjgC)
        call matrix_mul(transconjgC,S,C,matrix)

        call diag_matrix(matrix,ew)

        deallocate(transconjgC)

        ! Check if some eigenvalues are below threshold. If yes, get rid of them.
        iDel = 0
        do i=1,nmo2
          if (ew(iDel+1) .lt. ThrDel) iDel = iDel + 1
        enddo
        if (iDel == 0) then
          do i=1,nmo2
            ew(i) = 1.0/sqrt(sqrt(ew(i)))
          enddo
          do i=1,nmo2
            matrix(:,i,:)=matrix(:,i,:)*ew(i)
          enddo
          allocate(transconjgmatrix(nmo2,nmo2,rcq))
          allocate(Cnew(nmo1,nmo2,rcq))
          call transpose_conjg(matrix,transconjgmatrix)
          call matrix_mul(C,matrix,transconjgmatrix,Cnew)
          C = Cnew
          deallocate(transconjgmatrix,matrix,ew,Cnew)
        else
          print*,"   ** WARNING: deleted", iDel, " eigenvalues in S^{-1/2} construction."
          allocate(matrix2(nmo2,nmo2-iDel,rcq))
          do i=1,nmo2-iDel
            ew(i+iDel) = 1.0/sqrt(sqrt(ew(i+iDel)))
          enddo
          do i=1,nmo2-iDel
            matrix2(:,i,:)=matrix(:,i+iDel,:)*ew(i+iDel)
          enddo
          allocate(transconjgmatrix(nmo2-iDel,nmo2,rcq))
          allocate(Cnew(nmo1,nmo2,rcq))
          call transpose_conjg(matrix2,transconjgmatrix)
          call matrix_mul(C,matrix2,transconjgmatrix,Cnew)
          C = Cnew
          deallocate(transconjgmatrix,matrix,ew)
          deallocate(matrix2,Cnew)
        endif

        ! Compute C^\dagger * S * C with the new C matrix. (only to check that it is equal to identity))
        allocate(matrix(nmo2,nmo2,rcq))
        allocate(transconjgC(nmo2,nmo1,rcq))
        call transpose_conjg(C,transconjgC)
        call matrix_mul(transconjgC,S,C,matrix)
        deallocate(transconjgC)

        ! Check orthogonalization:
        allocate(Id(nmo2,nmo2,rcq))
        Id = 0.D0
        forall(i=1:nmo2) Id(i,i,1) = 1.D0
        call diff_norm2_rcq(matrix, Id, norm)
        
        if (norm .gt. Thrsh) then
           write(luout,*)
           write (luout,'(A,E12.4)') " **WARNING** Symmetrical orthogonalization has failed (thrsh: 1.D-12), norm:",norm
        else
           !write (luout,'(A,E12.4)') " **SUCCESS** Symmetrical orthogonalization was successful, norm:",norm
        end if

   end subroutine SymOrth

   subroutine matrix_invsqrt (M)
        !Return the inverse square root of the input matrix M
        ! Matrix should be square and positive definite

        implicit none
        real(8), parameter           :: Thrsh=1.D-12
        real(8),    intent(inout)    :: M(:,:,:) ! assumed shape array
        real(8),    allocatable      :: matrix(:,:,:),transconjgmatrix(:,:,:)
        real(8),    allocatable      :: ew(:)
        integer                      :: i, n, rcq

        rcq  = size(M,3)
        n    = size(M,1)
        if (size(M,2) /= n) then
           write(luout,*)"matrix_invsqrt: wrong dimensions, should be square matrix."
           stop
        end if

        allocate(matrix(n,n,rcq))
        allocate(ew(n))
        matrix = M
        call diag_matrix(matrix,ew)

        if (minval(ew) < Thrsh) then
           write(luout,*)"matrix_invsqrt: input matrix is not positive definite."
           write(luout,*)"smallest eigenvalue is ",minval(ew)
           stop
        end if

        forall (i=1:n) matrix(:,i,:) = matrix(:,i,:)/sqrt(sqrt(ew(i)))
        allocate(transconjgmatrix(n,n,rcq))
        call transpose_conjg(matrix,transconjgmatrix)
        call matrix_mul(matrix,transconjgmatrix,M)

   end subroutine matrix_invsqrt

   subroutine matrix_inv (M)
        !Return the inverse of the input matrix M
        ! Matrix should be square and positive definite

        implicit none
        real(8), parameter           :: Thrsh=1.D-12
        real(8),    intent(inout)    :: M(:,:,:) ! assumed shape array
        real(8),    allocatable      :: matrix(:,:,:),transconjgmatrix(:,:,:)
        real(8),    allocatable      :: ew(:)
        integer                      :: i, n, rcq

        rcq  = size(M,3)
        n    = size(M,1)
        if (size(M,2) /= n) then
           write(luout,*)"matrix_inv: wrong dimensions, should be square matrix."
           stop
        end if

        allocate(matrix(n,n,rcq))
        allocate(ew(n))
        matrix = M
        call diag_matrix(matrix,ew)

        if (minval(ew) < Thrsh) then
           write(luout,*)"matrix_inv: input matrix is not positive definite."
           write(luout,*)"smallest eigenvalue is ",minval(ew)
           stop
        end if

        forall (i=1:n) matrix(:,i,:) = matrix(:,i,:)/sqrt(ew(i))
        allocate(transconjgmatrix(n,n,rcq))
        call transpose_conjg(matrix,transconjgmatrix)
        call matrix_mul(matrix,transconjgmatrix,M)

   end subroutine matrix_inv

   subroutine sort_by_ascending_eigenvalues(vec,val,restricted)
        ! Sort the vectors by lowest eigenvalues
        
        implicit none
        integer                :: i
        real(8), intent(inout) :: vec(:,:,:)
        real(8), intent(inout) :: val(:)
        logical, intent(in)    :: restricted
        real(8), allocatable   :: vec_tmp(:,:,:)
        real(8), allocatable   :: vec_alpha(:,:,:)
        real(8), allocatable   :: vec_beta(:,:,:)
        real(8), allocatable   :: tmp(:,:,:)
        real(8), allocatable   :: val_tmp(:)
        real(8), allocatable   :: val_alpha(:)
        real(8), allocatable   :: val_beta(:)
 
        if (size(vec,2) .ne. size(val)) stop 'sort_by_ascending_eigenvalues: error in dimension of vec and val'
 
        allocate(vec_tmp(size(vec,1),size(vec,2),size(vec,3)))
        allocate(val_tmp(size(val)))
        if (restricted) then
           vec_tmp = 0.D0
           val_tmp = 0.D0
           vec_tmp = vec
           val_tmp = val
           do i=1,size(val)
              val_tmp(i)   = val(minloc(val(:),1))
              vec_tmp(:,i,:) = vec(:,minloc(val(:),1),:)
              val(minloc(val(:),1)) = 1.d9
           end do
        else
           allocate(vec_alpha(size(vec,1),size(vec,2)/2,size(vec,3)))
           allocate(vec_beta(size(vec,1),size(vec,2)/2,size(vec,3)))
           allocate(val_alpha(size(val)/2))
           allocate(val_beta(size(val)/2))
           vec_alpha = 0.D0
           vec_beta  = 0.D0
           val_alpha = 0.D0
           val_beta  = 0.D0
           do i=1,size(val_alpha)
              vec_alpha(:,i,:) = vec(:,2*i-1,:)
              vec_beta(:,i,:)  = vec(:,2*i,:)
              val_alpha(i) = val(2*i-1)
              val_beta(i)  = val(2*i)
           enddo
           do i=1,size(val_alpha)
              val_tmp(2*i-1) = val_alpha(minloc(val_alpha(:),1))
              val_tmp(2*i) = val_beta(minloc(val_beta(:),1))
              vec_tmp(:,2*i-1,:) = vec_alpha(:,minloc(val_alpha(:),1),:)
              vec_tmp(:,2*i,:) = vec_beta(:,minloc(val_beta(:),1),:)
              val_alpha(minloc(val_alpha(:),1)) = 1.d9
              val_beta(minloc(val_beta(:),1)) = 1.d9
           enddo
           deallocate(vec_alpha,vec_beta,val_alpha,val_beta)
        endif
 
        vec = vec_tmp
        val = val_tmp
        deallocate(val_tmp,vec_tmp)

   end subroutine sort_by_ascending_eigenvalues

   subroutine sort_by_descending_eigenvalues(vec,val,restricted)
        
        implicit none
        integer                :: i
        real(8), intent(inout) :: vec(:,:,:)
        real(8), intent(inout) :: val(:)
        logical, intent(in)    :: restricted
        real(8), allocatable   :: val_tmp(:)
 
        if (size(vec,2) .ne. size(val)) stop 'sort_by_descending_eigenvalues: error in dimension of vec and val'
 
        allocate(val_tmp(size(val)))
        val_tmp = -val
        call sort_by_ascending_eigenvalues(vec,val_tmp,restricted)
        val = val_tmp
        deallocate(val_tmp)

   end subroutine sort_by_descending_eigenvalues

   subroutine TestLinAlg
        !Driver for unit tests of this module.

!        use rose_utils

        real(8),    allocatable      :: C(:,:,:),S(:,:,:), U(:,:,:), VT(:,:,:), eigval(:)
        integer                      :: i, j, m, n, icq
        real(8)                      :: A(6,5),D(5,6),F(5,5),G(5,3)
        complex(8)                   :: B(3,4),E(4,3),H(4,4),K(4,2)

        DATA             A/&
       &  8.79, 6.11,-9.15, 9.57,-3.49, 9.84,&
       &  9.93, 6.91,-7.93, 1.64, 4.02, 0.15,&
       &  9.83, 5.04, 4.86, 8.83, 9.80,-8.99,&
       &  5.45,-0.27, 4.85, 0.74,10.00,-6.02,&
       &  3.16, 7.98, 3.01, 5.80, 4.27,-5.31&
       &                  /
        DATA             B/&
       & ( 5.91,-5.69),(-3.15,-4.08),(-4.89, 4.20),&
       & ( 7.09, 2.72),(-1.89, 3.27),( 4.10,-6.70),&
       & ( 7.78,-4.06),( 4.57,-2.07),( 3.28,-3.84),&
       & (-0.79,-7.21),(-3.88,-3.30),( 3.84, 1.19)&
       &                  /
        DATA             D/&
       &  8.79, 6.11,-9.15, 9.57,-3.49, 9.84,&
       &  9.93, 6.91,-7.93, 1.64, 4.02, 0.15,&
       &  9.83, 5.04, 4.86, 8.83, 9.80,-8.99,&
       &  5.45,-0.27, 4.85, 0.74,10.00,-6.02,&
       &  3.16, 7.98, 3.01, 5.80, 4.27,-5.31&
       &                  /
        DATA             E/&
       & ( 5.91,-5.69),(-3.15,-4.08),(-4.89, 4.20),&
       & ( 7.09, 2.72),(-1.89, 3.27),( 4.10,-6.70),&
       & ( 7.78,-4.06),( 4.57,-2.07),( 3.28,-3.84),&
       & (-0.79,-7.21),(-3.88,-3.30),( 3.84, 1.19)&
       &                  /
        DATA             F/&
       &  3.14, 0.00, 0.00, 0.00, 0.00,&
       &  0.17, 0.79, 0.00, 0.00, 0.00,&
       & -0.90, 0.83, 4.53, 0.00, 0.00,&
       &  1.65,-0.65,-3.70, 5.32, 0.00,&
       & -0.72, 0.28, 1.60,-1.37, 1.98&
       &                  /
        DATA             G/&
       & -7.29, 9.25, 5.99,-1.94,-8.30,&
       &  6.11, 2.90,-5.05,-3.80, 9.66,&
       &  0.59, 8.88, 7.57, 5.57,-1.67&
       &                  /
        DATA             H/&
       & ( 5.96, 0.00),( 0.40, 1.19),(-0.83, 0.48),(-0.57,-0.40),&
       & ( 0.00, 0.00),( 7.95, 0.00),( 0.33,-0.09),( 0.22,-0.74),&
       & ( 0.00, 0.00),( 0.00, 0.00),( 4.43, 0.00),(-1.09,-0.32),&
       & ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 3.46, 0.00)&
       &                  /
        DATA             K/&
       & (-2.94, 5.79),( 8.12,-9.12),( 9.09,-5.03),( 7.36, 6.77),&
       & ( 8.44, 3.07),( 1.00,-4.62),( 3.64,-2.33),( 8.04, 2.87)&
       &                  /

        print*, " ENTERING UNIT TEST FOR LINEAR ALGEBRA"
        ! test symmetric orthogonalization for matrices of different size
        do icq = 1, 4
        if (icq .eq. 3) cycle ! No algebra for dimension 4, we only have real, complex, quaternion
        print*, " Testing algebra :",icq
        do m = 2, 30 ! loop over the size of the S matrix
           do n = m, m  ! loop over the number of vectors to be orthogonalized
              print*, " Symmetric orthogonalization for m=",m,", n=",n
              allocate (S(m,m,icq))
              allocate (C(m,n,icq))
              ! first basis function overlaps with all others that are otherwise orthogonal
              S = 0.D0
              j = 1
              call random_number(S)
              S = S*0.05D0
              forall(i=1:m) S(i,i,1) = 1.D0
              forall(i=1:m) S(i,i,2:icq) = 0.D0
              !forall(i=2:m) S(i,j,1) = 0.1D0
              !forall(i=2:m) S(j,i,1) = 0.1D0
              do i = 2, m
                 do j = 1, i
                    S(j,i,1) = S(i,j,1)
                    S(j,i,2:icq) = - S(i,j,2:icq)
                 enddo
              enddo
              ! define unit vectors with some kind of phase
              C = 0.D0
              do j = 1, n
              !   C(j,j,mod(j,icq)+1) = 1.D0
                 C(j,j,1) = 1.D0
              end do
              call SymOrth(C,S)
              deallocate (S)
              deallocate (C)
           enddo
        enddo
        enddo

        allocate(S(6,5,1))
        allocate(eigval(5))
        allocate(U(6,6,1))
        allocate(VT(5,5,1))
        S(:,:,1) = A

        call lapack_SVD(S,eigval,U,VT)
        print*
        print*,"(m,n) = (6,5)"
        print*,"eigval SVD:",eigval
        print*,"U SVD:"
        !call print_matrix(U)
        print*,"VT SVD:"
        !call print_matrix(VT)
        call jacobi_SVD(S,eigval,U,VT)
        print*
        print*,"eigvaltest:",eigval
        print*,"U test:"
        !call print_matrix(U)
        print*,"VT test:"
        !call print_matrix(VT)
        deallocate(S,eigval,U,VT)

        allocate(S(5,6,1))
        allocate(eigval(5))
        allocate(U(5,5,1))
        allocate(VT(6,6,1))
        S(:,:,1) = D

        call lapack_SVD(S,eigval,U,VT)
        print*
        print*,"(m,n) = (5,6)"
        print*,"eigval SVD:",eigval
        print*,"U SVD:"
        !call print_matrix(U)
        print*,"VT SVD:"
        !call print_matrix(VT)
        call jacobi_SVD(S,eigval,U,VT)
        print*
        print*,"eigvaltest:",eigval
        print*,"U test:"
        !call print_matrix(U)
        print*,"VT test:"
        !call print_matrix(VT)
        deallocate(S,eigval,U,VT)

        allocate(S(4,3,2))
        allocate(eigval(3))
        allocate(U(4,4,2))
        allocate(VT(3,3,2))
        S(:,:,1) = dreal(E(:,:))
        S(:,:,2) = dimag(E(:,:))

        call lapack_SVD(S,eigval,U,VT)
        print*
        print*,"(m,n) = (4,3)"
        print*,"eigval SVD:",eigval
        print*,"U SVD:"
        !call print_matrix(U)
        print*,"VT SVD:"
        !call print_matrix(VT)
        call jacobi_SVD(S,eigval,U,VT)
        print*
        print*,"eigvaltest:",eigval
        print*,"U test:"
        !call print_matrix(U)
        print*,"VT test:"
        !call print_matrix(VT)
        deallocate(S,eigval,U,VT)

        allocate(S(3,4,2))
        allocate(eigval(3))
        allocate(U(3,3,2))
        allocate(VT(4,4,2))
        S(:,:,1) = dreal(B(:,:))
        S(:,:,2) = dimag(B(:,:))

        call lapack_SVD(S,eigval,U,VT)
        print*
        print*,"(m,n) = (3,4)"
        print*,"eigval SVD:",eigval
        print*,"U SVD:"
        !call print_matrix(U)
        print*,"VT SVD:"
        !call print_matrix(VT)
        call jacobi_SVD(S,eigval,U,VT)
        print*
        print*,"eigvaltest:",eigval
        print*,"U test:"
        !call print_matrix(U)
        print*,"VT test:"
        !call print_matrix(VT)
        deallocate(S,eigval,U,VT)

        allocate(S(5,5,1))
        allocate(C(5,3,1))
        allocate(U(5,3,1))
        S(:,:,1) = F
        C(:,:,1) = G
        call lapack_solvelineq(S,C,U)
        print*,"lapack test solve lineq, X:"
        !call print_matrix(U)
        call Gauss_Seidel_solvelineq(S,C,U) ! SU = C
        print*,"Gauss Seidel test solve lineq, X:"
        !call print_matrix(U)
        deallocate(S,C,U)

        allocate(S(4,4,2))
        allocate(C(4,2,2))
        allocate(U(4,2,2))
        S(:,:,1) = dreal(H)
        S(:,:,2) = dimag(H)
        C(:,:,1) = dreal(K)
        C(:,:,2) = dimag(K)
        call lapack_solvelineq(S,C,U)
        print*,"lapack test solve lineq, X:"
        !call print_matrix(U)
        call Gauss_Seidel_solvelineq(S,C,U) ! SU = C
        print*,"Gauss Seidel test solve lineq, X:"
        !call print_matrix(U)
        deallocate(S,C,U)

        print*, " COMPLETED UNIT TEST FOR LINEAR ALGEBRA"

   end subroutine TestLinAlg

end module linear_algebra
