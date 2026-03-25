module scf_module
    use iso_c_binding
    use one_eints
    use libcint_interface
    use math_utils
    implicit none

    integer,  parameter :: MAX_DIIS    = 8
    real(c_double), parameter :: SCF_TOL_E   = 1.0d-9
    real(c_double), parameter :: SCF_TOL_GRD = 3.16d-5
    integer,  parameter :: MAX_ITER    = 100
    real(c_double), parameter :: SCHWARZ_TOL = 1.0d-12
    real(c_double), parameter :: ADIIS_THRESH = 0.1d0   ! switch to DIIS below this

    logical :: df_ready = .false.
    integer :: df_naux = 0
    real(c_double), allocatable :: df_B(:,:,:)

contains

    ! ================================================================
    ! Nuclear repulsion energy
    ! ================================================================
    function nuclear_repulsion() result(E_nuc)
        implicit none
        real(c_double) :: E_nuc
        integer :: A, B, pA, pB
        real(c_double) :: ZA, ZB, dR(3)
        E_nuc = 0.0d0
        do A = 1, natm
            ZA = dble(atm((A-1)*6 + 1))
            pA = atm((A-1)*6 + 2) + 1
            do B = 1, A-1
                ZB = dble(atm((B-1)*6 + 1))
                pB = atm((B-1)*6 + 2) + 1
                dR = env(pA:pA+2) - env(pB:pB+2)
                E_nuc = E_nuc + ZA * ZB / sqrt(sum(dR**2))
            end do
        end do
    end function nuclear_repulsion

    ! ================================================================
    ! Pre-compute Cauchy-Schwarz bounds Q(I,J) = max|(IJ|IJ)|^{1/2}
    ! ================================================================
    subroutine compute_schwarz(Q)
        implicit none
        real(c_double), intent(out) :: Q(nbas, nbas)
        integer :: shI, shJ, di, dj
        integer(c_int) :: shls(4), cdims(4), info
        real(c_double), allocatable :: buf(:)
        Q = 0.0d0
        do shI = 1, nbas
            di = ao_loc(shI+1) - ao_loc(shI)
            do shJ = 1, shI
                dj = ao_loc(shJ+1) - ao_loc(shJ)
                allocate(buf(di*dj*di*dj))
                buf = 0.0d0
                cdims = [int(di,c_int), int(dj,c_int), int(di,c_int), int(dj,c_int)]
                shls = [shI-1, shJ-1, shI-1, shJ-1]
                info = cint2e_sph(buf, cdims, shls, atm, natm, bas, nbas, &
                                  env, c_null_ptr, c_null_ptr)
                Q(shI,shJ) = sqrt(maxval(abs(buf)))
                Q(shJ,shI) = Q(shI,shJ)
                deallocate(buf)
            end do
        end do
    end subroutine compute_schwarz

    ! ================================================================
    ! Direct JK build with shell-quartet symmetry and OpenMP.
    !
    ! Coulomb and exchange must use different canonical shell pairings:
    !   J  from unique quartets (IJ|KL)
    !   K  from unique quartets (IK|JL)
    ! Re-using (IJ|KL) for exchange is incorrect because (IJ|KL) is not
    ! related to (IK|JL) by the 8-fold permutation symmetry.
    !
    ! F = Hcore + J - K    (K already includes factor 0.5)
    ! ================================================================

    subroutine build_fock_direct(P, Hcore, Q, F)
        implicit none
        real(c_double), intent(in)  :: P(nao,nao), Hcore(nao,nao), Q(nbas,nbas)
        real(c_double), intent(out) :: F(nao,nao)
        integer :: shI, shJ, shK, shL, ai, aj, ak, al
        integer :: mi, mj, mk, ml, di, dj, dk, dl
        logical :: need_j, need_k
        integer(c_int) :: shls(4), cdims(4), info
        real(c_double), allocatable :: J_mat(:,:), K_mat(:,:)
        real(c_double) :: Pmax, g
        integer :: max_ao_per_shell, bufsize
        real(c_double), allocatable :: buf_j(:), buf_k(:)

        ! Find largest shell dimension for pre-allocating workspace
        max_ao_per_shell = 0
        do shI = 1, nbas
            max_ao_per_shell = max(max_ao_per_shell, ao_loc(shI+1) - ao_loc(shI))
        end do
        bufsize = max_ao_per_shell**4

        allocate(J_mat(nao,nao), K_mat(nao,nao))
        J_mat = 0.0d0;  K_mat = 0.0d0
        Pmax = maxval(abs(P))

        !$omp parallel default(none) &
        !$omp   shared(nbas, nao, ao_loc, atm, natm, bas, env, P, Q, Pmax, bufsize) &
        !$omp   private(shI, shJ, shK, shL, di, dj, dk, dl, ai, aj, ak, al, &
        !$omp           mi, mj, mk, ml, need_j, need_k, shls, cdims, info, buf_j, buf_k, g) &
        !$omp   reduction(+: J_mat, K_mat)
        allocate(buf_j(bufsize), buf_k(bufsize))
        !$omp do schedule(dynamic)
        do shI = 1, nbas
            di = ao_loc(shI+1) - ao_loc(shI)
            do shJ = 1, nbas
                dj = ao_loc(shJ+1) - ao_loc(shJ)
                do shK = 1, nbas
                    dk = ao_loc(shK+1) - ao_loc(shK)
                    do shL = 1, nbas
                        dl = ao_loc(shL+1) - ao_loc(shL)
                        need_j = (Q(shI,shJ)*Q(shK,shL)*Pmax >= SCHWARZ_TOL)
                        need_k = (Q(shI,shK)*Q(shJ,shL)*Pmax >= SCHWARZ_TOL)
                        if (.not. need_j .and. .not. need_k) cycle

                        if (need_j) then
                            buf_j(1:di*dj*dk*dl) = 0.0d0
                            cdims = [int(di,c_int), int(dj,c_int), int(dk,c_int), int(dl,c_int)]
                            shls = [shI-1, shJ-1, shK-1, shL-1]
                            info = cint2e_sph(buf_j, cdims, shls, atm, natm, bas, nbas, &
                                              env, c_null_ptr, c_null_ptr)

                            do al = 1, dl
                                ml = ao_loc(shL) + al - 1
                                do ak = 1, dk
                                    mk = ao_loc(shK) + ak - 1
                                    do aj = 1, dj
                                        mj = ao_loc(shJ) + aj - 1
                                        do ai = 1, di
                                            mi = ao_loc(shI) + ai - 1
                                            g = buf_j((((al-1)*dk + ak-1)*dj + aj-1)*di + ai)
                                            J_mat(mi,mj) = J_mat(mi,mj) + P(mk,ml) * g
                                        end do
                                    end do
                                end do
                            end do
                        end if

                        if (need_k) then
                            buf_k(1:di*dk*dj*dl) = 0.0d0
                            cdims = [int(di,c_int), int(dk,c_int), int(dj,c_int), int(dl,c_int)]
                            shls = [shI-1, shK-1, shJ-1, shL-1]
                            info = cint2e_sph(buf_k, cdims, shls, atm, natm, bas, nbas, &
                                              env, c_null_ptr, c_null_ptr)

                            do al = 1, dl
                                ml = ao_loc(shL) + al - 1
                                do aj = 1, dj
                                    mj = ao_loc(shJ) + aj - 1
                                    do ak = 1, dk
                                        mk = ao_loc(shK) + ak - 1
                                        do ai = 1, di
                                            mi = ao_loc(shI) + ai - 1
                                            g = buf_k((((al-1)*dj + aj-1)*dk + ak-1)*di + ai)
                                            K_mat(mi,mj) = K_mat(mi,mj) + 0.5d0 * P(mk,ml) * g
                                        end do
                                    end do
                                end do
                            end do
                        end if
                    end do
                end do
            end do
        end do
        !$omp end do
        deallocate(buf_j, buf_k)
        !$omp end parallel
        F = Hcore + J_mat - K_mat
        deallocate(J_mat, K_mat)
    end subroutine build_fock_direct

    ! ================================================================
    ! Packed AO-pair index for mu >= nu
    ! ================================================================
    pure integer function pair_index(mu, nu) result(p)
        implicit none
        integer, intent(in) :: mu, nu
        p = mu * (mu - 1) / 2 + nu
    end function pair_index

    ! ================================================================
    ! Build DF factors B(mu,nu,Q) from AO ERI metric over unique pairs
    ! (reference implementation; suitable for small/medium systems)
    ! ================================================================
    subroutine initialize_df_factors(chol_tol, max_rank)
        implicit none
        real(c_double), intent(in), optional :: chol_tol
        integer, intent(in), optional :: max_rank

        integer :: shI, shJ, shK, shL
        integer :: di, dj, dk, dl
        integer :: ai, aj, ak, al
        integer :: mi, mj, mk, ml, i
        integer :: p, q, npair, rank_cap, k, piv
        real(c_double) :: tol, g, pivot_val
        integer(c_int) :: shls(4), cdims(4), info
        integer :: max_ao_per_shell, bufsize
        real(c_double), allocatable :: buf(:), metric(:,:), L(:,:), resid(:), work(:)
        integer, allocatable :: pair_mu(:), pair_nu(:)

        if (df_ready) return

        tol = 1.0d-10
        if (present(chol_tol)) tol = chol_tol

        npair = nao * (nao + 1) / 2
        rank_cap = npair
        if (present(max_rank)) rank_cap = min(max_rank, npair)

        allocate(pair_mu(npair), pair_nu(npair))
        p = 0
        do mi = 1, nao
            do mj = 1, mi
                p = p + 1
                pair_mu(p) = mi
                pair_nu(p) = mj
            end do
        end do

        allocate(metric(npair,npair))
        metric = 0.0d0

        max_ao_per_shell = 0
        do shI = 1, nbas
            max_ao_per_shell = max(max_ao_per_shell, ao_loc(shI+1) - ao_loc(shI))
        end do
        bufsize = max_ao_per_shell**4
        allocate(buf(bufsize))

        do shI = 1, nbas
            di = ao_loc(shI+1) - ao_loc(shI)
            do shJ = 1, nbas
                dj = ao_loc(shJ+1) - ao_loc(shJ)
                do shK = 1, nbas
                    dk = ao_loc(shK+1) - ao_loc(shK)
                    do shL = 1, nbas
                        dl = ao_loc(shL+1) - ao_loc(shL)
                        cdims = [int(di,c_int), int(dj,c_int), int(dk,c_int), int(dl,c_int)]
                        shls  = [shI-1, shJ-1, shK-1, shL-1]
                        info  = cint2e_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                        do al = 1, dl
                            ml = ao_loc(shL) + al - 1
                            do ak = 1, dk
                                mk = ao_loc(shK) + ak - 1
                                q = pair_index(max(mk,ml), min(mk,ml))
                                do aj = 1, dj
                                    mj = ao_loc(shJ) + aj - 1
                                    do ai = 1, di
                                        mi = ao_loc(shI) + ai - 1
                                        p = pair_index(max(mi,mj), min(mi,mj))
                                        g = buf((((al-1)*dk + ak-1)*dj + aj-1)*di + ai)
                                        metric(p,q) = g
                                        metric(q,p) = g
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
        deallocate(buf)

        allocate(L(npair, rank_cap), resid(npair), work(npair))
        L = 0.0d0
        do i = 1, npair
            resid(i) = max(0.0d0, metric(i,i))
        end do
        df_naux = 0

        do k = 1, rank_cap
            piv = maxloc(resid, 1)
            pivot_val = resid(piv)
            if (pivot_val < tol) exit

            work = metric(:,piv)
            if (k > 1) work = work - matmul(L(:,1:k-1), L(piv,1:k-1))

            L(:,k) = work / sqrt(pivot_val)
            resid = resid - L(:,k) * L(:,k)
            where (resid < 0.0d0) resid = 0.0d0
            df_naux = k
        end do

        if (allocated(df_B)) deallocate(df_B)
        allocate(df_B(nao,nao,df_naux))
        df_B = 0.0d0

        do k = 1, df_naux
            do p = 1, npair
                mi = pair_mu(p)
                mj = pair_nu(p)
                g = L(p,k)
                df_B(mi,mj,k) = g
                df_B(mj,mi,k) = g
            end do
        end do

        df_ready = .true.

        deallocate(metric, L, resid, work, pair_mu, pair_nu)
    end subroutine initialize_df_factors

    ! ================================================================
    ! Density-fitted Fock build
    ! ERI approx: (mu nu|la si) ~= sum_Q B(mu,nu,Q) B(la,si,Q)
    ! ================================================================
    subroutine build_fock_df(P, Hcore, F)
        implicit none
        real(c_double), intent(in)  :: P(nao,nao), Hcore(nao,nao)
        real(c_double), intent(out) :: F(nao,nao)

        integer :: q
        real(c_double), allocatable :: J_mat(:,:), K_mat(:,:), T(:,:), z(:)

        if (.not. df_ready) call initialize_df_factors()
        if (df_naux <= 0) stop "DF initialization failed: no auxiliary rank"

        allocate(J_mat(nao,nao), K_mat(nao,nao), T(nao,nao), z(df_naux))
        J_mat = 0.0d0
        K_mat = 0.0d0

        do q = 1, df_naux
            z(q) = sum(P * df_B(:,:,q))
            J_mat = J_mat + z(q) * df_B(:,:,q)
        end do

        do q = 1, df_naux
            call dgemm('N','N',nao,nao,nao, 1.0d0, df_B(:,:,q),nao, P,nao, 0.0d0, T,nao)
            call dgemm('N','T',nao,nao,nao, 1.0d0, T,nao, df_B(:,:,q),nao, 1.0d0, K_mat,nao)
        end do

        F = Hcore + J_mat - 0.5d0 * K_mat

        deallocate(J_mat, K_mat, T, z)
    end subroutine build_fock_df

    ! ================================================================
    ! Release cached DF tensors
    ! ================================================================
    subroutine clear_df_cache()
        implicit none
        if (allocated(df_B)) deallocate(df_B)
        df_naux = 0
        df_ready = .false.
    end subroutine clear_df_cache

    ! ================================================================
    ! DIIS error: e = FPS - SPF  (return max-norm)
    ! ================================================================
    subroutine compute_diis_error(F, P, S, err, err_norm)
        implicit none
        real(c_double), intent(in)  :: F(nao,nao), P(nao,nao), S(nao,nao)
        real(c_double), intent(out) :: err(nao,nao), err_norm
        real(c_double), allocatable :: tmp1(:,:), tmp2(:,:)
        allocate(tmp1(nao,nao), tmp2(nao,nao))
        ! tmp1 = F*P, tmp2 = S*P
        call dgemm('N','N',nao,nao,nao, 1.0d0, F,nao, P,nao, 0.0d0, tmp1,nao)
        call dgemm('N','N',nao,nao,nao, 1.0d0, S,nao, P,nao, 0.0d0, tmp2,nao)
        ! err = (FP)*S - (SP)*F
        call dgemm('N','N',nao,nao,nao,  1.0d0, tmp1,nao, S,nao, 0.0d0, err,nao)
        call dgemm('N','N',nao,nao,nao, -1.0d0, tmp2,nao, F,nao, 1.0d0, err,nao)
        err_norm = maxval(abs(err))
        deallocate(tmp1, tmp2)
    end subroutine compute_diis_error

    ! ================================================================
    ! Update P from F:  F_orth -> C_orth -> C -> P = 2 C_occ C_occ^T
    ! ================================================================
    subroutine update_density(F, X, nocc, P, e_orb, C_out)
        implicit none
        real(c_double), intent(in)  :: F(nao,nao), X(nao,nao)
        integer, intent(in) :: nocc
        real(c_double), intent(out) :: P(nao,nao), e_orb(nao)
        real(c_double), intent(out), optional :: C_out(nao,nao)
        real(c_double), allocatable :: Fo(:,:), Co(:,:), C(:,:), tmp(:,:), wk(:)
        real(c_double) :: wk1(1)
        integer(c_int) :: lwork, info, i

        allocate(Fo(nao,nao), Co(nao,nao), C(nao,nao), tmp(nao,nao))
        call dgemm('N','N',nao,nao,nao, 1.0d0, F,nao, X,nao, 0.0d0, tmp,nao)
        call dgemm('T','N',nao,nao,nao, 1.0d0, X,nao, tmp,nao, 0.0d0, Fo,nao)
        Co = Fo
        lwork = -1
        call dsyev('V','L',nao,Co,nao, e_orb, wk1, lwork, info)
        lwork = int(wk1(1))
        allocate(wk(lwork))
        call dsyev('V','L',nao,Co,nao, e_orb, wk, lwork, info)
        if (info /= 0) stop "dsyev failed in update_density"
        deallocate(wk)
        call dgemm('N','N',nao,nao,nao, 1.0d0, X,nao, Co,nao, 0.0d0, C,nao)
        if (present(C_out)) C_out = C
        P = 0.0d0
        do i = 1, nocc
            call dgemm('N','T',nao,nao,1, 2.0d0, C(:,i),nao, C(:,i),nao, 1.0d0, P,nao)
        end do
        deallocate(Fo, Co, C, tmp)
    end subroutine update_density

    ! ================================================================
    ! Pulay DIIS extrapolation
    ! ================================================================
    subroutine diis_extrapolate(ns, F_hist, e_hist, F_new)
        implicit none
        integer, intent(in) :: ns
        real(c_double), intent(in)  :: F_hist(nao,nao,MAX_DIIS)
        real(c_double), intent(in)  :: e_hist(nao,nao,MAX_DIIS)
        real(c_double), intent(out) :: F_new(nao,nao)
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
        B(ns+1, ns+1) =  0.0d0
        rhs = 0.0d0
        rhs(ns+1) = -1.0d0

        call dgesv(int(ndim,c_int), 1_c_int, B, int(ndim,c_int), ipiv, &
                   rhs, int(ndim,c_int), info_d)

        F_new = 0.0d0
        do i = 1, ns
            F_new = F_new + rhs(i) * F_hist(:,:,i)
        end do
        deallocate(B, rhs, c, ipiv)
    end subroutine diis_extrapolate

    ! ================================================================
    ! ADIIS (Hu & Yang JCP 2010) — Frank-Wolfe minimization on simplex
    ! E[c] = 2 sum_i c_i Tr[D_i * F_n] + sum_ij c_i c_j Tr[D_i * F_j]
    ! ================================================================
    subroutine adiis_extrapolate(ns, F_hist, P_hist, F_n, P_n, F_new)
        implicit none
        integer, intent(in) :: ns
        real(c_double), intent(in)  :: F_hist(nao,nao,MAX_DIIS)
        real(c_double), intent(in)  :: P_hist(nao,nao,MAX_DIIS)
        real(c_double), intent(in)  :: F_n(nao,nao), P_n(nao,nao)
        real(c_double), intent(out) :: F_new(nao,nao)
        real(c_double), allocatable :: z(:), A(:,:), c(:), grad(:), d(:)
        real(c_double) :: gamma, numer, denom
        integer :: i, j, k_fw, iter_fw

        allocate(z(ns), A(ns,ns), c(ns), grad(ns), d(ns))

        ! Build linear term z_i = Tr[(D_i - D_n) * F_n]
        do i = 1, ns
            z(i) = sum((P_hist(:,:,i) - P_n) * F_n)
        end do

        ! Build quadratic term A_ij = Tr[(D_i - D_n) * (F_j - F_n)]
        do i = 1, ns
            do j = 1, ns
                A(i,j) = sum((P_hist(:,:,i) - P_n) * (F_hist(:,:,j) - F_n))
            end do
        end do

        ! Initialize c = uniform
        c = 1.0d0 / dble(ns)

        ! Frank-Wolfe (conditional gradient) iterations on simplex
        do iter_fw = 1, 200
            grad = z + matmul(A, c)
            k_fw = minloc(grad, 1)
            d = -c;  d(k_fw) = d(k_fw) + 1.0d0   ! e_k - c
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

    ! ================================================================
    ! Main SCF driver
    ! ================================================================
    subroutine run_scf(S, Hcore, X, P_in, C_out)
        implicit none
        real(c_double), intent(in)    :: S(nao,nao), Hcore(nao,nao), X(nao,nao)
        real(c_double), intent(inout) :: P_in(nao,nao)
        real(c_double), intent(out), optional :: C_out(nao,nao)

        real(c_double), allocatable :: F(:,:), P(:,:), err(:,:)
        real(c_double), allocatable :: F_hist(:,:,:), P_hist(:,:,:), e_hist(:,:,:)
        real(c_double), allocatable :: F_extrap(:,:), Q(:,:), e_orb(:)
        real(c_double), allocatable :: C_last(:,:)
        real(c_double) :: E_elec, E_tot, E_prev, err_norm, E_nuc
        integer :: iter, ns, nocc, env_stat, env_len
        logical :: use_df
        character(len=32) :: fock_method

        allocate(F(nao,nao), P(nao,nao), err(nao,nao), F_extrap(nao,nao))
        allocate(Q(nbas,nbas), e_orb(nao), C_last(nao,nao))
        allocate(F_hist(nao,nao,MAX_DIIS), P_hist(nao,nao,MAX_DIIS), e_hist(nao,nao,MAX_DIIS))

        nocc  = total_electrons / 2
        E_nuc = nuclear_repulsion()
        P     = P_in
        E_prev= huge(1.0d0)
        ns    = 0
        C_last = 0.0d0
        use_df = .false.
        fock_method = 'DIRECT'

        call get_environment_variable('EXAGRAD_FOCK_BUILDER', fock_method, length=env_len, status=env_stat)
        if (env_stat == 0) then
            fock_method = adjustl(fock_method(1:env_len))
            if (trim(fock_method) == 'df' .or. trim(fock_method) == 'DF' .or. &
                trim(fock_method) == 'density_fitting' .or. trim(fock_method) == 'DENSITY_FITTING') then
                use_df = .true.
            end if
        end if

        print *, ""
        if (use_df) then
            print *, "  ------ RHF SCF (DF-JK) ------"
        else
            print *, "  ------ RHF SCF (Direct-JK) ------"
        end if
        print *, "  E_nuc =", E_nuc, " Ha"
        print "(A5, A20, A18, A12)", "  Iter", "E_total (Ha)", "|dE|", "|FPS-SPF|"

        if (.not. use_df) then
            call compute_schwarz(Q)
        else
            call initialize_df_factors()
            print "(A, I8)", "  DF auxiliary rank =", df_naux
        end if

        do iter = 1, MAX_ITER
            ! Build Fock from current density
            if (use_df) then
                call build_fock_df(P, Hcore, F)
            else
                call build_fock_direct(P, Hcore, Q, F)
            end if

            ! Electronic energy
            E_elec = 0.5d0 * sum(P * (Hcore + F))
            E_tot  = E_elec + E_nuc

            ! DIIS error vector
            call compute_diis_error(F, P, S, err, err_norm)

            print "(I5, F20.10, E18.4, E12.3)", iter, E_tot, abs(E_tot-E_prev), err_norm

            ! Convergence check
            if (abs(E_tot - E_prev) < SCF_TOL_E .and. err_norm < SCF_TOL_GRD) then
                print *, ""
                print *, "  SCF converged!"
                print "(A, F18.10, A)", "  E_total =", E_tot, " Ha"
                P_in = P
                if (present(C_out)) C_out = C_last
                if (use_df) call clear_df_cache()
                deallocate(F, P, err, F_extrap, Q, e_orb, C_last, F_hist, P_hist, e_hist)
                return
            end if

            ! Push into history (circular buffer)
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

            ! Extrapolate Fock:
            !  - iter <= 2: ADIIS (energy-based, useful far from convergence)
            !  - iter > 2:  pure Pulay DIIS
            if (ns >= 2 .and. err_norm > ADIIS_THRESH) then
                call adiis_extrapolate(ns, F_hist, P_hist, F, P, F_extrap)
            else if (ns >= 2) then
                call diis_extrapolate(ns, F_hist, e_hist, F_extrap)
            else
                F_extrap = F
            end if

            ! Update density from extrapolated Fock
            E_prev = E_tot
            call update_density(F_extrap, X, nocc, P, e_orb, C_last)

        end do

        print *, "  WARNING: SCF did not converge in", MAX_ITER, "cycles"
        P_in = P
        if (present(C_out)) C_out = C_last
        if (use_df) call clear_df_cache()
        deallocate(F, P, err, F_extrap, Q, e_orb, C_last, F_hist, P_hist, e_hist)
    end subroutine run_scf

end module scf_module
