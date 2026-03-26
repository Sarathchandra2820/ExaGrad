module scf_module
    use iso_c_binding
    use omp_lib, only: omp_get_max_threads
    use one_eints, only: activate_molecule
    use molecule_t, only: molecule
    use math_utils
    use fock_builder_module, only: nuclear_repulsion, normalize_fock_method, fock_method_banner, &
                                   initialize_fock_backend, build_fock, clear_fock_backend, cholesky_naux, true_df_naux
    implicit none

    integer,  parameter :: MAX_DIIS = 8
    real(c_double), parameter :: SCF_TOL_E = 1.0d-9
    real(c_double), parameter :: SCF_TOL_GRD = 3.16d-5
    integer,  parameter :: MAX_ITER = 100
    real(c_double), parameter :: ADIIS_THRESH = 0.1d0

contains

    subroutine compute_diis_error(F, P, S, err, err_norm)
        implicit none
        real(c_double), intent(in)  :: F(:,:), P(:,:), S(:,:)
        real(c_double), intent(out) :: err(:,:), err_norm
        ! Persistent workspace: allocated once, reused across SCF iterations
        real(c_double), allocatable, save :: tmp1(:,:), tmp2(:,:)
        integer, save :: n_saved = 0
        integer :: n

        n = size(F, 1)
        if (n /= n_saved) then
            if (allocated(tmp1)) deallocate(tmp1, tmp2)
            allocate(tmp1(n,n), tmp2(n,n))
            n_saved = n
        end if
        call dgemm('N','N',n,n,n, 1.0d0, F,n, P,n, 0.0d0, tmp1,n)
        call dgemm('N','N',n,n,n, 1.0d0, S,n, P,n, 0.0d0, tmp2,n)
        call dgemm('N','N',n,n,n,  1.0d0, tmp1,n, S,n, 0.0d0, err,n)
        call dgemm('N','N',n,n,n, -1.0d0, tmp2,n, F,n, 1.0d0, err,n)
        err_norm = maxval(abs(err))
    end subroutine compute_diis_error

    subroutine update_density(F, X, nocc, P, e_orb, C_out)
        implicit none
        real(c_double), intent(in)  :: F(:,:), X(:,:)
        integer, intent(in) :: nocc
        real(c_double), intent(out) :: P(:,:), e_orb(:)
        real(c_double), intent(out), optional :: C_out(:,:)
        real(c_double), allocatable :: Fo(:,:), Co(:,:), C(:,:), tmp(:,:), wk(:)
        real(c_double) :: wk1(1)
        integer(c_int) :: lwork, info
        integer :: n

        n = size(F, 1)
        allocate(Fo(n,n), Co(n,n), C(n,n), tmp(n,n))
        call dgemm('N','N',n,n,n, 1.0d0, F,n, X,n, 0.0d0, tmp,n)
        call dgemm('T','N',n,n,n, 1.0d0, X,n, tmp,n, 0.0d0, Fo,n)
        Co = Fo
        lwork = -1
        call dsyev('V','L',n,Co,n, e_orb, wk1, lwork, info)
        lwork = int(wk1(1))
        allocate(wk(lwork))
        call dsyev('V','L',n,Co,n, e_orb, wk, lwork, info)
        if (info /= 0) stop 'dsyev failed in update_density'
        deallocate(wk)
        call dgemm('N','N',n,n,n, 1.0d0, X,n, Co,n, 0.0d0, C,n)
        if (present(C_out)) C_out = C
        ! Single dgemm: P = 2 * C_occ * C_occ^T  (replaces the old per-column loop)
        call dgemm('N','T', n, n, nocc, 2.0d0, C(:,1:nocc), n, C(:,1:nocc), n, 0.0d0, P, n)
        deallocate(Fo, Co, C, tmp)
    end subroutine update_density

    subroutine diis_extrapolate(ns, F_hist, e_hist, F_new)
        implicit none
        integer, intent(in) :: ns
        real(c_double), intent(in)  :: F_hist(:,:,:)
        real(c_double), intent(in)  :: e_hist(:,:,:)
        real(c_double), intent(out) :: F_new(:,:)
        integer :: i, j, ndim
        real(c_double), allocatable :: B(:,:), rhs(:), c(:)
        integer(c_int), allocatable :: ipiv(:)
        integer(c_int) :: info_d

        ndim = ns + 1
        allocate(B(ndim,ndim), rhs(ndim), c(ndim), ipiv(ndim))
        B = 0.0d0
        do i = 1, ns
            do j = 1, ns
                B(i,j) = sum(e_hist(:,:,i) * e_hist(:,:,j))
            end do
        end do
        B(ns+1, 1:ns) = -1.0d0
        B(1:ns, ns+1) = -1.0d0
        B(ns+1, ns+1) = 0.0d0
        rhs = 0.0d0
        rhs(ns+1) = -1.0d0

        call dgesv(int(ndim,c_int), 1_c_int, B, int(ndim,c_int), ipiv, rhs, int(ndim,c_int), info_d)

        F_new = 0.0d0
        do i = 1, ns
            F_new = F_new + rhs(i) * F_hist(:,:,i)
        end do
        deallocate(B, rhs, c, ipiv)
    end subroutine diis_extrapolate

    subroutine adiis_extrapolate(ns, F_hist, P_hist, F_n, P_n, F_new)
        implicit none
        integer, intent(in) :: ns
        real(c_double), intent(in)  :: F_hist(:,:,:)
        real(c_double), intent(in)  :: P_hist(:,:,:)
        real(c_double), intent(in)  :: F_n(:,:), P_n(:,:)
        real(c_double), intent(out) :: F_new(:,:)
        real(c_double), allocatable :: z(:), A(:,:), c(:), grad(:), d(:)
        real(c_double), allocatable :: dP(:,:,:), dF(:,:,:)
        real(c_double) :: gamma, numer, denom
        integer :: i, j, k_fw, iter_fw, n

        n = size(F_n, 1)
        allocate(z(ns), A(ns,ns), c(ns), grad(ns), d(ns))

        ! Precompute dP_i = P_hist_i - P_n and dF_i = F_hist_i - F_n once
        ! avoids re-computing n^2 subtractions in every (i,j) inner-loop step
        allocate(dP(n,n,ns), dF(n,n,ns))
        do i = 1, ns
            dP(:,:,i) = P_hist(:,:,i) - P_n
            dF(:,:,i) = F_hist(:,:,i) - F_n
        end do
        do i = 1, ns
            z(i) = sum(dP(:,:,i) * F_n)
            do j = 1, ns
                A(i,j) = sum(dP(:,:,i) * dF(:,:,j))
            end do
        end do
        deallocate(dP, dF)

        c = 1.0d0 / dble(ns)
        do iter_fw = 1, 200
            grad = z + matmul(A, c)
            ! Exit when gradient is flat (all entries equal) — no improvement possible
            if (maxval(grad) - minval(grad) < 1.0d-12) exit
            k_fw = minloc(grad, 1)
            d = -c
            d(k_fw) = d(k_fw) + 1.0d0
            numer = -dot_product(grad, d)
            denom = dot_product(d, matmul(A, d))
            if (denom > 1.0d-14) then
                gamma = min(1.0d0, max(0.0d0, numer/denom))
            else
                gamma = 1.0d0
            end if
            c = c + gamma * d
            if (gamma < 1.0d-10) exit
        end do

        F_new = 0.0d0
        do i = 1, ns
            F_new = F_new + c(i) * F_hist(:,:,i)
        end do
        deallocate(z, A, c, grad, d)
    end subroutine adiis_extrapolate

    subroutine run_scf(mol, S, Hcore, X, P_in, C_out)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double), intent(in)    :: S(:,:), Hcore(:,:), X(:,:)
        real(c_double), intent(inout) :: P_in(:,:)
        real(c_double), intent(out), optional :: C_out(:,:)

        real(c_double), allocatable :: F(:,:), P(:,:), err(:,:)
        real(c_double), allocatable :: F_hist(:,:,:), P_hist(:,:,:), e_hist(:,:,:)
        real(c_double), allocatable :: F_extrap(:,:), Q(:,:), e_orb(:), C_last(:,:)
        real(c_double) :: E_elec, E_tot, E_prev, err_norm, E_nuc
        integer :: iter, ns, nocc, env_stat, env_len, n, nb
        character(len=32) :: method_env, fock_method
        character(len=40) :: banner

        call activate_molecule(mol)
        n = size(S, 1)
        nb = int(mol%basis%nbas)

        allocate(F(n,n), P(n,n), err(n,n), F_extrap(n,n))
        allocate(Q(nb,nb), e_orb(n), C_last(n,n))
        allocate(F_hist(n,n,MAX_DIIS), P_hist(n,n,MAX_DIIS), e_hist(n,n,MAX_DIIS))

        nocc = mol%total_electrons / 2
        E_nuc = nuclear_repulsion(mol)
        P = P_in
        E_prev = huge(1.0d0)
        ns = 0
        C_last = 0.0d0
        fock_method = 'direct'

        call get_environment_variable('EXAGRAD_FOCK_BUILDER', method_env, length=env_len, status=env_stat)
        if (env_stat == 0) then
            method_env = adjustl(method_env(1:env_len))
            call normalize_fock_method(trim(method_env), fock_method)
        end if
        banner = fock_method_banner(trim(fock_method))

        print *, ''
        print *, '  ------ ' // trim(banner) // ' ------'
        print '(A, I4)', '  OpenMP max threads =', omp_get_max_threads()
        print *, '  E_nuc =', E_nuc, ' Ha'
        print '(A5, A20, A18, A12)', '  Iter', 'E_total (Ha)', '|dE|', '|FPS-SPF|'

        call initialize_fock_backend(trim(fock_method), Q)
        if (trim(fock_method) == 'cholesky' .or. trim(fock_method) == 'block_cholesky') then
            print '(A, I8)', '  Cholesky auxiliary rank =', cholesky_naux
        else if (trim(fock_method) == 'true_df') then
            print '(A, I8)', '  True-DF auxiliary rank =', true_df_naux
        end if

        do iter = 1, MAX_ITER
            call build_fock(trim(fock_method), P, Hcore, Q, F)
            E_elec = 0.5d0 * sum(P * (Hcore + F))
            E_tot = E_elec + E_nuc
            call compute_diis_error(F, P, S, err, err_norm)
            print '(I5, F20.10, E18.4, E12.3)', iter, E_tot, abs(E_tot-E_prev), err_norm

            if (abs(E_tot - E_prev) < SCF_TOL_E .and. err_norm < SCF_TOL_GRD) then
                print *, ''
                print *, '  SCF converged!'
                print '(A, F18.10, A)', '  E_total =', E_tot, ' Ha'
                P_in = P
                if (present(C_out)) C_out = C_last
                call clear_fock_backend(trim(fock_method))
                deallocate(F, P, err, F_extrap, Q, e_orb, C_last, F_hist, P_hist, e_hist)
                return
            end if

            if (ns < MAX_DIIS) then
                ns = ns + 1
            else
                F_hist(:,:,1:MAX_DIIS-1) = F_hist(:,:,2:MAX_DIIS)
                P_hist(:,:,1:MAX_DIIS-1) = P_hist(:,:,2:MAX_DIIS)
                e_hist(:,:,1:MAX_DIIS-1) = e_hist(:,:,2:MAX_DIIS)
            end if
            F_hist(:,:,ns) = F
            P_hist(:,:,ns) = P
            e_hist(:,:,ns) = err

            if (ns >= 2 .and. err_norm > ADIIS_THRESH) then
                call adiis_extrapolate(ns, F_hist, P_hist, F, P, F_extrap)
            else if (ns >= 2) then
                call diis_extrapolate(ns, F_hist, e_hist, F_extrap)
            else
                F_extrap = F
            end if

            E_prev = E_tot
            call update_density(F_extrap, X, nocc, P, e_orb, C_last)
        end do

        print *, '  WARNING: SCF did not converge in', MAX_ITER, 'cycles'
        P_in = P
        if (present(C_out)) C_out = C_last
        call clear_fock_backend(trim(fock_method))
        deallocate(F, P, err, F_extrap, Q, e_orb, C_last, F_hist, P_hist, e_hist)
    end subroutine run_scf

end module scf_module
