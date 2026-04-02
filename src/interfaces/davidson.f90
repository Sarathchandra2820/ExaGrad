module davidson_module

    use math_utils
    use iso_c_binding

    implicit none

    abstract interface
        subroutine matvec_iface(x, Ax, n)
            import :: c_double
            integer, intent(in) :: n
            real(c_double), intent(in)  :: x(n)
            real(c_double), intent(out) :: Ax(n)
        end subroutine

        subroutine precond_eigen_iface(r, z, n, shift)
            import :: c_double
            integer, intent(in) :: n
            real(c_double), intent(in)  :: r(n)
            real(c_double), intent(out) :: z(n)
            real(c_double), intent(in)  :: shift
        end subroutine
    end interface

contains

    !==========================================================================
    ! Orthogonalise vector z against V(:,1:k) using dgemv (two passes).
    ! proj is workspace of length >= k.
    ! Returns the norm of z after orthogonalisation in znorm.
    !==========================================================================
    subroutine orthogonalise(V, n, k, z, proj, znorm)
        integer, intent(in) :: n, k
        real(c_double), intent(in)    :: V(n, *)
        real(c_double), intent(inout) :: z(n)
        real(c_double), intent(out)   :: proj(*), znorm

        ! Pass 1
        call dgemv('T', n, k, 1.0d0, V, n, z, 1, 0.0d0, proj, 1)
        call dgemv('N', n, k, -1.0d0, V, n, proj, 1, 1.0d0, z, 1)
        ! Pass 2 for numerical stability
        call dgemv('T', n, k, 1.0d0, V, n, z, 1, 0.0d0, proj, 1)
        call dgemv('N', n, k, -1.0d0, V, n, proj, 1, 1.0d0, z, 1)

        znorm = dnrm2(n, z, 1)
    end subroutine orthogonalise


    !==========================================================================
    ! Davidson iterative solver for linear equations: A * x = b
    !==========================================================================
    subroutine davidson_linear_solve(sigma_op, precond_op, b, x, n, &
                                      tol, max_iter, max_sub, converged)
        procedure(matvec_iface) :: sigma_op, precond_op
        integer, intent(in)           :: n, max_iter, max_sub
        real(c_double), intent(in)    :: b(n)
        real(c_double), intent(inout) :: x(n)
        real(c_double), intent(in)    :: tol
        logical, intent(out)          :: converged

        real(c_double), allocatable :: V(:,:), AV(:,:)
        real(c_double), allocatable :: H(:,:), Hcopy(:,:)
        real(c_double), allocatable :: rhs(:), coeff(:), proj(:)
        real(c_double), allocatable :: r(:), z(:)
        integer, allocatable        :: ipiv(:)

        real(c_double) :: rnorm, znorm
        integer :: iter, k, j, info

        allocate(V(n, max_sub), AV(n, max_sub))
        allocate(H(max_sub, max_sub), Hcopy(max_sub, max_sub))
        allocate(rhs(max_sub), coeff(max_sub), proj(max_sub))
        allocate(r(n), z(n), ipiv(max_sub))

        converged = .false.

        ! Seed: use preconditioned b if x is zero, else use x
        if (dnrm2(n, x, 1) < 1.0d-15) then
            call precond_op(b, z, n)
        else
            call dcopy(n, x, 1, z, 1)
        end if

        znorm = dnrm2(n, z, 1)
        if (znorm < 1.0d-15) then
            print *, '  Davidson: zero initial vector, aborting'
            deallocate(V, AV, H, Hcopy, rhs, coeff, proj, r, z, ipiv)
            return
        end if

        k = 1
        call dcopy(n, z, 1, V(1,1), 1)
        call dscal(n, 1.0d0 / znorm, V(1,1), 1)
        call sigma_op(V(1,1), AV(1,1), n)
        H(1,1) = ddot(n, V(1,1), 1, AV(1,1), 1)
        rhs(1) = ddot(n, V(1,1), 1, b, 1)

        print '(A)', ''
        print '(A)', '  --- Davidson Linear Solver ---'
        print '(A5, A18)', '  Iter', '||residual||'

        do iter = 1, max_iter

            ! Solve reduced system: H(1:k,1:k) * c = rhs(1:k)
            Hcopy(1:k,1:k) = H(1:k,1:k)
            coeff(1:k) = rhs(1:k)
            call dgesv(k, 1, Hcopy, max_sub, ipiv, coeff, max_sub, info)
            if (info /= 0) then
                print *, '  Davidson: dgesv failed, info =', info
                exit
            end if

            ! Solution: x = V * c
            call dgemv('N', n, k, 1.0d0, V, n, coeff, 1, 0.0d0, x, 1)

            ! Residual: r = AV * c - b
            call dgemv('N', n, k, 1.0d0, AV, n, coeff, 1, 0.0d0, r, 1)
            call daxpy(n, -1.0d0, b, 1, r, 1)

            rnorm = dnrm2(n, r, 1)
            print '(I5, E18.8)', iter, rnorm

            if (rnorm < tol) then
                converged = .true.
                exit
            end if

            ! Precondition the residual
            call precond_op(r, z, n)

            ! Orthogonalise via dgemv (two passes)
            call orthogonalise(V, n, k, z, proj, znorm)

            if (znorm < 1.0d-15) then
                print *, '  Davidson: new direction too small, stopping'
                exit
            end if

            ! Restart if subspace full
            if (k >= max_sub) then
                print *, '  Davidson: restarting'
                call dcopy(n, x, 1, z, 1)
                znorm = dnrm2(n, x, 1)
                k = 1
                call dcopy(n, z, 1, V(1,1), 1)
                call dscal(n, 1.0d0 / znorm, V(1,1), 1)
                call sigma_op(V(1,1), AV(1,1), n)
                H(1,1) = ddot(n, V(1,1), 1, AV(1,1), 1)
                rhs(1) = ddot(n, V(1,1), 1, b, 1)
                cycle
            end if

            ! Expand subspace
            k = k + 1
            call dcopy(n, z, 1, V(1,k), 1)
            call dscal(n, 1.0d0 / znorm, V(1,k), 1)
            call sigma_op(V(1,k), AV(1,k), n)

            ! Update projected matrix (symmetric) and RHS
            do j = 1, k
                H(j,k) = ddot(n, V(1,j), 1, AV(1,k), 1)
                H(k,j) = H(j,k)
            end do
            rhs(k) = ddot(n, V(1,k), 1, b, 1)

        end do

        if (converged) then
            print '(A, E12.4)', '  Converged, ||r|| =', rnorm
        else
            print *, '  WARNING: Davidson linear solve did not converge'
        end if
        print '(A)', ''

        deallocate(V, AV, H, Hcopy, rhs, coeff, proj, r, z, ipiv)
    end subroutine davidson_linear_solve


    !==========================================================================
    ! Davidson eigensolver for symmetric A: find nroots lowest eigenvalues
    !==========================================================================
    subroutine davidson_eigen_solve(sigma_op, precond_op, x, eigenvalues, n, nroots, &
                                     tol, max_iter, max_sub, nconv)
        procedure(matvec_iface)        :: sigma_op
        procedure(precond_eigen_iface) :: precond_op
        integer, intent(in)            :: n, nroots, max_iter, max_sub
        real(c_double), intent(inout)  :: x(n, nroots)
        real(c_double), intent(out)    :: eigenvalues(nroots)
        real(c_double), intent(in)     :: tol
        integer, intent(out)           :: nconv

        real(c_double), allocatable :: V(:,:), AV(:,:)
        real(c_double), allocatable :: H(:,:), Hcopy(:,:)
        real(c_double), allocatable :: theta(:)
        real(c_double), allocatable :: r(:), z(:), proj(:)
        real(c_double), allocatable :: V_tmp(:,:), AV_tmp(:,:)
        real(c_double), allocatable :: Ax_all(:,:), R_all(:,:)
        real(c_double), allocatable :: work(:)
        real(c_double), allocatable :: rnorms(:)
        logical, allocatable        :: conv_flag(:)

        real(c_double) :: rnorm_max, znorm
        integer :: iter, k, i, j, info, lwork, n_uncnv

        if (max_sub < 2 * nroots) then
            print *, '  Error: max_sub must be >= 2*nroots'
            nconv = 0
            return
        end if

        allocate(V(n, max_sub), AV(n, max_sub))
        allocate(H(max_sub, max_sub), Hcopy(max_sub, max_sub))
        allocate(theta(max_sub))
        allocate(r(n), z(n), proj(max_sub))
        allocate(V_tmp(n, nroots), AV_tmp(n, nroots))
        allocate(Ax_all(n, nroots), R_all(n, nroots))
        allocate(rnorms(nroots))
        allocate(conv_flag(nroots))
        lwork = 3 * max_sub
        allocate(work(lwork))

        nconv = 0
        H = 0.0d0
        eigenvalues = 0.0d0

        ! Seed subspace: orthonormalise initial guesses
        k = 0
        do i = 1, nroots
            call dcopy(n, x(1,i), 1, z, 1)
            if (k > 0) then
                call orthogonalise(V, n, k, z, proj, znorm)
            else
                znorm = dnrm2(n, z, 1)
            end if
            if (znorm < 1.0d-15) cycle
            k = k + 1
            call dcopy(n, z, 1, V(1,k), 1)
            call dscal(n, 1.0d0 / znorm, V(1,k), 1)
            call sigma_op(V(1,k), AV(1,k), n)
        end do

        if (k < nroots) then
            print *, '  Error: initial guesses are linearly dependent'
            nconv = 0
            deallocate(V, AV, H, Hcopy, theta, r, z, proj, V_tmp, AV_tmp, &
                       Ax_all, R_all, rnorms, conv_flag, work)
            return
        end if

        ! Build initial projected matrix
        do i = 1, k
            do j = 1, i
                H(j,i) = ddot(n, V(1,j), 1, AV(1,i), 1)
                H(i,j) = H(j,i)
            end do
        end do

        print '(A)', ''
        print '(A)', '  --- Davidson Eigensolver ---'
        print '(A5, A20, A15)', '  Iter', 'Eigenvalue(1)', '||max res||'

        do iter = 1, max_iter

            ! Diagonalise projected matrix
            Hcopy(1:k,1:k) = H(1:k,1:k)
            call dsyev('V', 'U', k, Hcopy, max_sub, theta, work, lwork, info)
            if (info /= 0) then
                print *, '  Davidson: dsyev failed, info =', info
                exit
            end if

            ! Batch Ritz vectors: x = V * S(:,1:nroots)  via dgemm
            call dgemm('N', 'N', n, nroots, k, 1.0d0, V, n, &
                       Hcopy, max_sub, 0.0d0, x, n)

            ! Batch A*Ritz: Ax_all = AV * S(:,1:nroots)  via dgemm
            call dgemm('N', 'N', n, nroots, k, 1.0d0, AV, n, &
                       Hcopy, max_sub, 0.0d0, Ax_all, n)

            ! Compute all residuals: R_all(:,i) = Ax_all(:,i) - theta(i) * x(:,i)
            ! and cache them for the expansion step
            nconv = 0
            conv_flag = .false.
            rnorm_max = 0.0d0

            do i = 1, nroots
                call dcopy(n, Ax_all(1,i), 1, R_all(1,i), 1)
                call daxpy(n, -theta(i), x(1,i), 1, R_all(1,i), 1)
                rnorms(i) = dnrm2(n, R_all(1,i), 1)
                if (rnorms(i) > rnorm_max) rnorm_max = rnorms(i)
                if (rnorms(i) < tol) then
                    conv_flag(i) = .true.
                    nconv = nconv + 1
                end if
            end do

            eigenvalues(1:nroots) = theta(1:nroots)
            print '(I5, F20.10, E15.4)', iter, theta(1), rnorm_max

            if (nconv >= nroots) exit

            ! Count unconverged roots
            n_uncnv = nroots - nconv

            ! Restart if subspace cannot accommodate new vectors
            if (k + n_uncnv > max_sub) then
                ! Collapse to nroots Ritz vectors without recomputing sigma
                ! V_tmp already = x, AV_tmp already = Ax_all (computed above)
                do j = 1, nroots
                    call dcopy(n, x(1,j), 1, V(1,j), 1)
                    call dcopy(n, Ax_all(1,j), 1, AV(1,j), 1)
                end do
                k = nroots
                ! Rebuild projected matrix (diagonal = theta after restart)
                H(1:k,1:k) = 0.0d0
                do i = 1, k
                    do j = 1, i
                        H(j,i) = ddot(n, V(1,j), 1, AV(1,i), 1)
                        H(i,j) = H(j,i)
                    end do
                end do
                print *, '  Davidson: restarting'
                cycle
            end if

            ! Expand subspace with preconditioned residuals of unconverged roots
            do i = 1, nroots
                if (conv_flag(i)) cycle

                ! Precondition cached residual
                call precond_op(R_all(1,i), z, n, theta(i))

                ! Orthogonalise via dgemv (two passes)
                call orthogonalise(V, n, k, z, proj, znorm)

                if (znorm < 1.0d-15) cycle

                k = k + 1
                call dcopy(n, z, 1, V(1,k), 1)
                call dscal(n, 1.0d0 / znorm, V(1,k), 1)
                call sigma_op(V(1,k), AV(1,k), n)

                ! Update projected matrix (symmetric)
                do j = 1, k
                    H(j,k) = ddot(n, V(1,j), 1, AV(1,k), 1)
                    H(k,j) = H(j,k)
                end do
            end do

        end do

        if (nconv >= nroots) then
            print '(A, I3, A)', '  All ', nroots, ' roots converged'
        else
            print '(A, I3, A, I3, A)', '  WARNING: ', nconv, ' of ', nroots, ' roots converged'
        end if
        print '(A)', ''

        deallocate(V, AV, H, Hcopy, theta, r, z, proj, V_tmp, AV_tmp, &
                   Ax_all, R_all, rnorms, conv_flag, work)
    end subroutine davidson_eigen_solve

end module davidson_module
