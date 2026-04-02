program test_davidson
    use iso_c_binding
    use davidson_module
    implicit none

    integer, parameter :: n = 10000
    integer, parameter :: nroots = 20
    real(c_double), allocatable :: A_mat(:,:), diag(:)
    real(c_double), allocatable :: b(:), x(:), x_exact(:)
    real(c_double), allocatable :: eig_guess(:,:), eigenvalues(:)
    real(c_double), allocatable :: eig_ref(:), work(:), A_copy(:,:)
    real(c_double) :: t_start, t_end, err
    logical :: converged
    integer :: i, j, nconv, info, lwork

    print '(A)', '================================================'
    print '(A, I6)', '  Davidson Solver Test, matrix size n =', n
    print '(A)', '================================================'

    ! -------------------------------------------------------------------
    ! Build a diagonally dominant symmetric matrix with known structure
    ! A(i,j) = 0.1/(|i-j|+1) for off-diagonal, A(i,i) = i (diagonal)
    ! -------------------------------------------------------------------
    allocate(A_mat(n,n), diag(n))

    do j = 1, n
        do i = 1, n
            if (i == j) then
                A_mat(i,j) = dble(i)
            else
                A_mat(i,j) = 0.1d0 / dble(abs(i-j) + 1)
            end if
        end do
    end do
    do i = 1, n
        diag(i) = A_mat(i,i)
    end do

    ! -------------------------------------------------------------------
    ! Test 1: Linear solver   A * x = b
    ! -------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '  ====== TEST 1: Linear Solve ======'

    allocate(b(n), x(n), x_exact(n))

    ! Build a known solution and compute b = A * x_exact
    do i = 1, n
        x_exact(i) = 1.0d0 / dble(i)
    end do
    call dgemv('N', n, n, 1.0d0, A_mat, n, x_exact, 1, 0.0d0, b, 1)

    x = 0.0d0
    call cpu_time(t_start)
    call davidson_linear_solve(apply_A, apply_precond, b, x, n, &
                                1.0d-10, 200, 50, converged)
    call cpu_time(t_end)

    ! Compute error vs exact solution
    call daxpy(n, -1.0d0, x_exact, 1, x, 1)
    err = dnrm2(n, x, 1)

    print '(A, L2)',     '  Converged:        ', converged
    print '(A, E12.4)',  '  ||x - x_exact||:  ', err
    print '(A, F10.4, A)', '  Wall time:        ', t_end - t_start, ' s'

    deallocate(b, x, x_exact)

    ! -------------------------------------------------------------------
    ! Test 2: Eigensolver   A * x = lambda * x
    ! -------------------------------------------------------------------
    print '(A)', ''
    print '(A)', '  ====== TEST 2: Eigensolver ======'

    ! Reference: full dsyev
    allocate(A_copy(n,n), eig_ref(n))
    A_copy = A_mat
    lwork = 3*n
    allocate(work(lwork))
    ! call cpu_time(t_start)
    ! call dsyev('V', 'U', n, A_copy, n, eig_ref, work, lwork, info)
    ! call cpu_time(t_end)
    ! print '(A, F10.4, A)', '  dsyev reference:  ', t_end - t_start, ' s'
    ! print '(A)', '  Reference lowest eigenvalues:'
    ! do i = 1, nroots
    !     print '(A, I3, A, F18.10)', '    lambda(', i, ') =', eig_ref(i)
    ! end do
    ! deallocate(work, A_copy)

    ! Davidson eigensolver
    allocate(eig_guess(n, nroots), eigenvalues(nroots))

    ! Initial guesses: unit vectors for smallest diagonal elements
    eig_guess = 0.0d0
    do i = 1, nroots
        eig_guess(i, i) = 1.0d0
    end do

    call cpu_time(t_start)
    call davidson_eigen_solve(apply_A, apply_eigen_precond, eig_guess, eigenvalues, &
                               n, nroots, 1.0d-10, 200, 10*nroots, nconv)
    call cpu_time(t_end)

    print '(A)', '  Davidson eigenvalues:'
    do i = 1, nroots
        print '(A, I3, A, F18.10, A, E12.4)', &
            '    lambda(', i, ') =', eigenvalues(i), &
            '   err =', abs(eigenvalues(i) - eig_ref(i))
    end do
    print '(A, I3, A, I3)',    '  Converged roots:  ', nconv, ' / ', nroots
    print '(A, F10.4, A)',     '  Wall time:        ', t_end - t_start, ' s'

    deallocate(eig_guess, eigenvalues, eig_ref, A_mat, diag)

    print '(A)', ''
    print '(A)', '  Done.'

contains

    ! A * x  (dense matvec using the stored matrix)
    subroutine apply_A(x_in, Ax_out, nn)
        integer, intent(in) :: nn
        real(c_double), intent(in)  :: x_in(nn)
        real(c_double), intent(out) :: Ax_out(nn)
        call dgemv('N', nn, nn, 1.0d0, A_mat, nn, x_in, 1, 0.0d0, Ax_out, 1)
    end subroutine

    ! Diagonal preconditioner for linear solve: M^{-1} = D^{-1}
    subroutine apply_precond(r_in, z_out, nn)
        integer, intent(in) :: nn
        real(c_double), intent(in)  :: r_in(nn)
        real(c_double), intent(out) :: z_out(nn)
        integer :: ii
        do ii = 1, nn
            z_out(ii) = r_in(ii) / diag(ii)
        end do
    end subroutine

    ! Diagonal preconditioner for eigensolver: (D - shift)^{-1}
    subroutine apply_eigen_precond(r_in, z_out, nn, shift)
        integer, intent(in) :: nn
        real(c_double), intent(in)  :: r_in(nn)
        real(c_double), intent(out) :: z_out(nn)
        real(c_double), intent(in)  :: shift
        integer :: ii
        real(c_double) :: denom
        do ii = 1, nn
            denom = diag(ii) - shift
            if (abs(denom) < 1.0d-12) then
                z_out(ii) = 0.0d0
            else
                z_out(ii) = r_in(ii) / denom
            end if
        end do
    end subroutine

end program test_davidson
