module fock_builder_module
    use iso_c_binding
    use omp_lib, only: omp_get_max_threads, omp_get_num_threads, omp_get_thread_num
    use one_eints
    use molecule_t, only: molecule
    use libcint_interface
    use math_utils
    implicit none

    real(c_double), parameter :: SCHWARZ_TOL = 1.0d-12

    logical :: cholesky_ready = .false.
    logical :: cholesky_omp_reported = .false.
    logical :: cholesky_block_omp_reported = .false.
    integer :: cholesky_naux = 0
    real(c_double), allocatable :: cholesky_B(:,:,:)

    logical :: true_df_ready = .false.
    logical :: true_df_omp_reported = .false.
    integer :: true_df_naux = 0
    real(c_double), allocatable :: true_df_B(:,:,:)

contains

    function nuclear_repulsion(mol) result(E_nuc)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double) :: E_nuc
        integer :: A, B, pA, pB
        real(c_double) :: ZA, ZB, dR(3)

        E_nuc = 0.0d0
        do A = 1, int(mol%basis%natm)
            ZA = dble(mol%basis%atm(1, A))
            pA = int(mol%basis%atm(2, A)) + 1
            do B = 1, A-1
                ZB = dble(mol%basis%atm(1, B))
                pB = int(mol%basis%atm(2, B)) + 1
                dR = mol%basis%env(pA:pA+2) - mol%basis%env(pB:pB+2)
                E_nuc = E_nuc + ZA * ZB / sqrt(sum(dR**2))
            end do
        end do
    end function nuclear_repulsion

    pure integer function pair_index(mu, nu) result(p)
        implicit none
        integer, intent(in) :: mu, nu
        p = mu * (mu - 1) / 2 + nu
    end function pair_index

    subroutine normalize_fock_method(method_in, method_out)
        implicit none
        character(len=*), intent(in) :: method_in
        character(len=*), intent(out) :: method_out
        character(len=64) :: tmp

        tmp = adjustl(trim(method_in))
        select case (trim(tmp))
        case ('cholesky', 'CHOLESKY', 'df', 'DF')
            method_out = 'cholesky'
        case ('block_cholesky', 'BLOCK_CHOLESKY', 'block_df', 'BLOCK_DF', 'df_block', 'DF_BLOCK', 'blocked_df', 'BLOCKED_DF')
            method_out = 'block_cholesky'
        case ('true_df', 'TRUE_DF', 'density_fitting', 'DENSITY_FITTING')
            method_out = 'true_df'
        case default
            method_out = 'direct'
        end select
    end subroutine normalize_fock_method

    pure function fock_method_banner(method) result(label)
        implicit none
        character(len=*), intent(in) :: method
        character(len=40) :: label
        select case (trim(method))
        case ('true_df')
            label = 'RHF SCF (True-DF-JK)'
        case ('block_cholesky')
            label = 'RHF SCF (Cholesky-Block-JK)'
        case ('cholesky')
            label = 'RHF SCF (Cholesky-JK)'
        case default
            label = 'RHF SCF (Direct-JK)'
        end select
    end function fock_method_banner

    subroutine compute_schwarz(Q)
        implicit none
        real(c_double), intent(out) :: Q(nbas, nbas)
        integer :: shI, shJ, di, dj, max_ao
        integer(c_int) :: shls(4), cdims(4), info
        real(c_double), allocatable :: buf(:)

        Q = 0.0d0
        ! Allocate once at max shell-quartet size, reuse across all pairs
        max_ao = maxval(ao_loc(2:nbas+1) - ao_loc(1:nbas))
        allocate(buf(max_ao**4))
        do shI = 1, nbas
            di = ao_loc(shI+1) - ao_loc(shI)
            do shJ = 1, shI
                dj = ao_loc(shJ+1) - ao_loc(shJ)
                buf(1:di*dj*di*dj) = 0.0d0
                cdims = [int(di,c_int), int(dj,c_int), int(di,c_int), int(dj,c_int)]
                shls = [shI-1, shJ-1, shI-1, shJ-1]
                info = cint2e_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                Q(shI,shJ) = sqrt(maxval(abs(buf(1:di*dj*di*dj))))
                Q(shJ,shI) = Q(shI,shJ)
            end do
        end do
        deallocate(buf)
    end subroutine compute_schwarz

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

        max_ao_per_shell = 0
        do shI = 1, nbas
            max_ao_per_shell = max(max_ao_per_shell, ao_loc(shI+1) - ao_loc(shI))
        end do
        bufsize = max_ao_per_shell**4

        allocate(J_mat(nao,nao), K_mat(nao,nao))
        J_mat = 0.0d0
        K_mat = 0.0d0
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
                            info = cint2e_sph(buf_j, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
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
                            info = cint2e_sph(buf_k, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
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

    subroutine initialize_cholesky_factors(chol_tol, max_rank)
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

        if (cholesky_ready) return
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

        print "(A, I4, A)", "  Building Cholesky metric (", omp_get_max_threads(), " threads) …"

        ! Each thread allocates buf once before the shell loop (not per-shI iteration).
        ! Atomics removed: ERI permutation symmetry guarantees concurrent writes of
        ! identical values to the same metric(p,q) element — last write wins safely.
        !$omp parallel default(none) &
        !$omp   shared(nbas, nao, ao_loc, atm, natm, bas, env, metric, bufsize) &
        !$omp   private(shI, shJ, shK, shL, di, dj, dk, dl, ai, aj, ak, al, &
        !$omp           mi, mj, mk, ml, p, q, g, shls, cdims, info, buf)
        allocate(buf(bufsize))
        !$omp do schedule(dynamic)
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
        !$omp end do
        deallocate(buf)
        !$omp end parallel

        allocate(L(npair, rank_cap), resid(npair), work(npair))
        L = 0.0d0
        do i = 1, npair
            resid(i) = max(0.0d0, metric(i,i))
        end do
        cholesky_naux = 0

        do k = 1, rank_cap
            piv = maxloc(resid, 1)
            pivot_val = resid(piv)
            if (pivot_val < tol) exit

            work = metric(:,piv)
            ! dgemv replaces the intrinsic matmul: work -= L(:,1:k-1) * L(piv,1:k-1)
            ! L(piv,1) with stride npair traverses row piv across columns (column-major)
            if (k > 1) call dgemv('N', npair, k-1, -1.0d0, L, npair, L(piv,1), npair, 1.0d0, work, 1)
            L(:,k) = work / sqrt(pivot_val)
            resid = resid - L(:,k) * L(:,k)
            where (resid < 0.0d0) resid = 0.0d0
            cholesky_naux = k
        end do

        if (allocated(cholesky_B)) deallocate(cholesky_B)
        allocate(cholesky_B(nao,nao,cholesky_naux))
        cholesky_B = 0.0d0

        do k = 1, cholesky_naux
            do p = 1, npair
                mi = pair_mu(p)
                mj = pair_nu(p)
                g = L(p,k)
                cholesky_B(mi,mj,k) = g
                cholesky_B(mj,mi,k) = g
            end do
        end do

        cholesky_ready = .true.
        deallocate(metric, L, resid, work, pair_mu, pair_nu)
    end subroutine initialize_cholesky_factors

    subroutine build_fock_cholesky(P, Hcore, F)
        implicit none
        real(c_double), intent(in)  :: P(nao,nao), Hcore(nao,nao)
        real(c_double), intent(out) :: F(nao,nao)
        integer :: q, team_threads
        real(c_double) :: zq
        real(c_double), allocatable :: J_mat(:,:), K_mat(:,:), T_local(:,:)

        if (.not. cholesky_ready) call initialize_cholesky_factors()
        if (cholesky_naux <= 0) stop 'Cholesky initialization failed: no auxiliary rank'

        allocate(J_mat(nao,nao), K_mat(nao,nao))
        J_mat = 0.0d0
        K_mat = 0.0d0
        team_threads = 1

        ! reduction(+:J_mat,K_mat): each thread accumulates privately, then summed —
        ! eliminates the serialising critical section
        !$omp parallel default(none) &
        !$omp   shared(cholesky_naux, nao, P, cholesky_B, team_threads) &
        !$omp   private(q, zq, T_local) &
        !$omp   reduction(+: J_mat, K_mat)
        allocate(T_local(nao,nao))

        !$omp single
        team_threads = omp_get_num_threads()
        !$omp end single

        !$omp do schedule(dynamic)
        do q = 1, cholesky_naux
            zq = sum(P * cholesky_B(:,:,q))
            J_mat = J_mat + zq * cholesky_B(:,:,q)
            call dgemm('N','N',nao,nao,nao, 1.0d0, cholesky_B(:,:,q),nao, P,nao, 0.0d0, T_local,nao)
            call dgemm('N','T',nao,nao,nao, 1.0d0, T_local,nao, cholesky_B(:,:,q),nao, 1.0d0, K_mat,nao)
        end do
        !$omp end do

        deallocate(T_local)
        !$omp end parallel

        if (.not. cholesky_omp_reported) then
            print "(A, I4)", "  Cholesky OpenMP team size =", team_threads
            cholesky_omp_reported = .true.
        end if

        F = Hcore + J_mat - 0.5d0 * K_mat
        deallocate(J_mat, K_mat)
    end subroutine build_fock_cholesky

    subroutine build_fock_cholesky_block(P, Hcore, F, block_q)
        implicit none
        real(c_double), intent(in)  :: P(nao,nao), Hcore(nao,nao)
        real(c_double), intent(out) :: F(nao,nao)
        integer, intent(in), optional :: block_q
        integer :: q, q0, q1, b, nblocks, bs
        integer :: tid, nthreads, team_threads
        real(c_double) :: zq
        real(c_double), allocatable :: J_mat(:,:), K_mat(:,:), T_local(:,:)

        if (.not. cholesky_ready) call initialize_cholesky_factors()
        if (cholesky_naux <= 0) stop 'Cholesky initialization failed: no auxiliary rank'

        bs = 16
        if (present(block_q)) bs = max(1, block_q)
        nblocks = (cholesky_naux + bs - 1) / bs

        allocate(J_mat(nao,nao), K_mat(nao,nao))
        J_mat = 0.0d0
        K_mat = 0.0d0
        team_threads = 1

        ! reduction(+:J_mat,K_mat): eliminates serialising critical section
        !$omp parallel default(none) &
        !$omp   shared(cholesky_naux, nao, P, cholesky_B, nblocks, bs, team_threads) &
        !$omp   private(q, q0, q1, b, tid, nthreads, zq, T_local) &
        !$omp   reduction(+: J_mat, K_mat)
        allocate(T_local(nao,nao))

        tid = omp_get_thread_num()
        nthreads = omp_get_num_threads()

        !$omp single
        team_threads = nthreads
        !$omp end single

        do b = tid + 1, nblocks, nthreads
            q0 = (b - 1) * bs + 1
            q1 = min(cholesky_naux, b * bs)
            do q = q0, q1
                zq = sum(P * cholesky_B(:,:,q))
                J_mat = J_mat + zq * cholesky_B(:,:,q)
                call dgemm('N','N',nao,nao,nao, 1.0d0, cholesky_B(:,:,q),nao, P,nao, 0.0d0, T_local,nao)
                call dgemm('N','T',nao,nao,nao, 1.0d0, T_local,nao, cholesky_B(:,:,q),nao, 1.0d0, K_mat,nao)
            end do
        end do

        deallocate(T_local)
        !$omp end parallel

        if (.not. cholesky_block_omp_reported) then
            print "(A, I4, A, I4)", "  Cholesky-Block OpenMP team size =", team_threads, ", Q-block =", bs
            cholesky_block_omp_reported = .true.
        end if

        F = Hcore + J_mat - 0.5d0 * K_mat
        deallocate(J_mat, K_mat)
    end subroutine build_fock_cholesky_block

    subroutine initialize_true_df_factors()
        implicit none
        integer :: ios, i, j, q, p, ai, aj, ap
        integer :: shI, shJ, shP, shQ
        integer :: di, dj, dp, dq, npair
        integer :: nbas_mol, nbas_aux, nbas_all
        integer :: natm_mol, natm_aux, natm_all
        integer :: env_all_size, nao_mol, nao_aux
        integer :: max_mol_shell, max_aux_shell
        integer :: aux_ao0, ip, mu, nu, info
        integer(c_int) :: shls2(2), shls3(3), cdims2(2), cdims3(3)
        integer(c_int), allocatable :: atm_all(:), bas_all(:)
        integer, allocatable :: ao_loc_all(:), pair_mu(:), pair_nu(:)
        real(c_double), allocatable :: env_all(:), metric(:,:), rhs(:,:), buf2c(:), buf3c(:)

        if (true_df_ready) return

        open(unit=61, file='ints/df_meta_info.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open ints/df_meta_info.txt'
        read(61, *, iostat=ios) natm_mol
        if (ios /= 0) stop 'Failed reading natm_mol from df_meta_info'
        read(61, *, iostat=ios) nbas_mol
        if (ios /= 0) stop 'Failed reading nbas_mol from df_meta_info'
        read(61, *, iostat=ios) natm_aux
        if (ios /= 0) stop 'Failed reading natm_aux from df_meta_info'
        read(61, *, iostat=ios) nbas_aux
        if (ios /= 0) stop 'Failed reading nbas_aux from df_meta_info'
        read(61, *, iostat=ios) natm_all
        if (ios /= 0) stop 'Failed reading natm_all from df_meta_info'
        read(61, *, iostat=ios) nbas_all
        if (ios /= 0) stop 'Failed reading nbas_all from df_meta_info'
        read(61, *, iostat=ios) env_all_size
        if (ios /= 0) stop 'Failed reading env_all_size from df_meta_info'
        read(61, *, iostat=ios) nao_mol
        if (ios /= 0) stop 'Failed reading nao_mol from df_meta_info'
        read(61, *, iostat=ios) nao_aux
        if (ios /= 0) stop 'Failed reading nao_aux from df_meta_info'
        close(61)

        if (nao_mol /= nao) stop 'True DF metadata nao_mol does not match current molecule'
        if (nbas_mol /= nbas) stop 'True DF metadata nbas_mol does not match current molecule'

        true_df_naux = nao_aux
        if (true_df_naux <= 0) stop 'True DF metadata has non-positive auxiliary AO count'

        allocate(atm_all(natm_all*6), bas_all(nbas_all*8), env_all(env_all_size))

        open(unit=62, file='ints/df_meta_atm.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open ints/df_meta_atm.txt'
        read(62, *, iostat=ios) atm_all
        if (ios /= 0) stop 'Failed reading ints/df_meta_atm.txt'
        close(62)

        open(unit=63, file='ints/df_meta_bas.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open ints/df_meta_bas.txt'
        read(63, *, iostat=ios) bas_all
        if (ios /= 0) stop 'Failed reading ints/df_meta_bas.txt'
        close(63)

        open(unit=64, file='ints/df_meta_env.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open ints/df_meta_env.txt'
        read(64, *, iostat=ios) env_all
        if (ios /= 0) stop 'Failed reading ints/df_meta_env.txt'
        close(64)

        allocate(ao_loc_all(nbas_all+1))
        ao_loc_all(1) = 1
        do i = 1, nbas_all
            ai = bas_all((i-1)*8 + 2)
            aj = bas_all((i-1)*8 + 4)
            ao_loc_all(i+1) = ao_loc_all(i) + (2*ai + 1) * aj
        end do

        aux_ao0 = ao_loc_all(nbas_mol+1) - 1
        npair = nao * (nao + 1) / 2

        allocate(metric(true_df_naux,true_df_naux), rhs(true_df_naux,npair))
        metric = 0.0d0
        rhs = 0.0d0

        max_mol_shell = 0
        do shI = 1, nbas_mol
            max_mol_shell = max(max_mol_shell, ao_loc_all(shI+1) - ao_loc_all(shI))
        end do
        max_aux_shell = 0
        do shP = 1, nbas_aux
            max_aux_shell = max(max_aux_shell, ao_loc_all(nbas_mol+shP+1) - ao_loc_all(nbas_mol+shP))
        end do

        allocate(buf2c(max_aux_shell*max_aux_shell))
        allocate(buf3c(max_mol_shell*max_mol_shell*max_aux_shell))

        do shP = 1, nbas_aux
            shI = nbas_mol + shP
            dp = ao_loc_all(shI+1) - ao_loc_all(shI)
            do shQ = 1, nbas_aux
                shJ = nbas_mol + shQ
                dq = ao_loc_all(shJ+1) - ao_loc_all(shJ)
                cdims2 = [int(dp,c_int), int(dq,c_int)]
                shls2  = [int(shI-1,c_int), int(shJ-1,c_int)]
                buf2c(1:dp*dq) = 0.0d0
                info = cint2c2e_sph(buf2c, cdims2, shls2, atm_all, int(natm_all,c_int), bas_all, int(nbas_all,c_int), &
                                   env_all, c_null_ptr, c_null_ptr)
                do ap = 1, dp
                    i = ao_loc_all(shI) + ap - 1 - aux_ao0
                    do aj = 1, dq
                        j = ao_loc_all(shJ) + aj - 1 - aux_ao0
                        metric(i,j) = buf2c((aj-1)*dp + ap)
                    end do
                end do
            end do
        end do

        do shI = 1, nbas_mol
            di = ao_loc_all(shI+1) - ao_loc_all(shI)
            do shJ = 1, nbas_mol
                dj = ao_loc_all(shJ+1) - ao_loc_all(shJ)
                do shP = 1, nbas_aux
                    shQ = nbas_mol + shP
                    dp = ao_loc_all(shQ+1) - ao_loc_all(shQ)
                    cdims3 = [int(di,c_int), int(dj,c_int), int(dp,c_int)]
                    shls3  = [int(shI-1,c_int), int(shJ-1,c_int), int(shQ-1,c_int)]
                    buf3c(1:di*dj*dp) = 0.0d0
                    info = cint3c2e_sph(buf3c, cdims3, shls3, atm_all, int(natm_all,c_int), bas_all, int(nbas_all,c_int), &
                                       env_all, c_null_ptr, c_null_ptr)
                    do ap = 1, dp
                        ip = ao_loc_all(shQ) + ap - 1 - aux_ao0
                        do aj = 1, dj
                            j = ao_loc_all(shJ) + aj - 1
                            do ai = 1, di
                                i = ao_loc_all(shI) + ai - 1
                                p = pair_index(max(i,j), min(i,j))
                                rhs(ip,p) = buf3c(((ap-1)*dj + aj-1)*di + ai)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call dpotrf('L', int(true_df_naux,c_int), metric, int(true_df_naux,c_int), info)
        if (info /= 0) stop 'dpotrf failed while factoring true DF metric'
        call dtrsm('L', 'L', 'N', 'N', int(true_df_naux,c_int), int(npair,c_int), 1.0d0, &
                   metric, int(true_df_naux,c_int), rhs, int(true_df_naux,c_int))

        if (allocated(true_df_B)) deallocate(true_df_B)
        allocate(true_df_B(nao, nao, true_df_naux))
        true_df_B = 0.0d0

        allocate(pair_mu(npair), pair_nu(npair))
        p = 0
        do mu = 1, nao
            do nu = 1, mu
                p = p + 1
                pair_mu(p) = mu
                pair_nu(p) = nu
            end do
        end do

        do q = 1, true_df_naux
            do p = 1, npair
                i = pair_mu(p)
                j = pair_nu(p)
                true_df_B(i,j,q) = rhs(q,p)
                true_df_B(j,i,q) = rhs(q,p)
            end do
        end do

        deallocate(pair_mu, pair_nu)
        deallocate(buf2c, buf3c)
        deallocate(metric, rhs)
        deallocate(ao_loc_all, atm_all, bas_all, env_all)

        true_df_ready = .true.
        print "(A, I8)", "  True-DF auxiliary rank =", true_df_naux
        print "(A)", "  True-DF source: metadata + Fortran-built CDERI"
    end subroutine initialize_true_df_factors

    subroutine build_fock_true_df(P, Hcore, F)
        implicit none
        real(c_double), intent(in)  :: P(nao,nao), Hcore(nao,nao)
        real(c_double), intent(out) :: F(nao,nao)
        integer :: q, team_threads
        real(c_double) :: zq
        real(c_double), allocatable :: J_mat(:,:), K_mat(:,:), T_local(:,:)

        if (.not. true_df_ready) call initialize_true_df_factors()

        allocate(J_mat(nao,nao), K_mat(nao,nao))
        J_mat = 0.0d0
        K_mat = 0.0d0
        team_threads = 1

        ! reduction(+:J_mat,K_mat): eliminates serialising critical section
        !$omp parallel default(none) &
        !$omp   shared(true_df_naux, nao, P, true_df_B, team_threads) &
        !$omp   private(q, zq, T_local) &
        !$omp   reduction(+: J_mat, K_mat)
        allocate(T_local(nao,nao))

        !$omp single
        team_threads = omp_get_num_threads()
        !$omp end single

        !$omp do schedule(dynamic)
        do q = 1, true_df_naux
            zq = sum(P * true_df_B(:,:,q))
            J_mat = J_mat + zq * true_df_B(:,:,q)
            call dgemm('N','N',nao,nao,nao, 1.0d0, true_df_B(:,:,q),nao, P,nao, 0.0d0, T_local,nao)
            call dgemm('N','T',nao,nao,nao, 1.0d0, T_local,nao, true_df_B(:,:,q),nao, 1.0d0, K_mat,nao)
        end do
        !$omp end do

        deallocate(T_local)
        !$omp end parallel

        if (.not. true_df_omp_reported) then
            print "(A, I4)", "  True-DF OpenMP team size =", team_threads
            true_df_omp_reported = .true.
        end if

        F = Hcore + J_mat - 0.5d0 * K_mat
        deallocate(J_mat, K_mat)
    end subroutine build_fock_true_df

    subroutine initialize_fock_backend(method, Q)
        implicit none
        character(len=*), intent(in) :: method
        real(c_double), intent(out) :: Q(nbas,nbas)

        select case (trim(method))
        case ('true_df')
            Q = 0.0d0
            call initialize_true_df_factors()
        case ('block_cholesky', 'cholesky')
            Q = 0.0d0
            call initialize_cholesky_factors()
        case default
            call compute_schwarz(Q)
        end select
    end subroutine initialize_fock_backend

    subroutine build_fock(method, P, Hcore, Q, F)
        implicit none
        character(len=*), intent(in) :: method
        real(c_double), intent(in) :: P(nao,nao), Hcore(nao,nao), Q(nbas,nbas)
        real(c_double), intent(out) :: F(nao,nao)

        select case (trim(method))
        case ('true_df')
            call build_fock_true_df(P, Hcore, F)
        case ('block_cholesky')
            call build_fock_cholesky_block(P, Hcore, F)
        case ('cholesky')
            call build_fock_cholesky(P, Hcore, F)
        case default
            call build_fock_direct(P, Hcore, Q, F)
        end select
    end subroutine build_fock

    subroutine clear_fock_backend(method)
        implicit none
        character(len=*), intent(in) :: method

        select case (trim(method))
        case ('true_df')
            if (allocated(true_df_B)) deallocate(true_df_B)
            true_df_naux = 0
            true_df_ready = .false.
            true_df_omp_reported = .false.
        case ('block_cholesky', 'cholesky')
            if (allocated(cholesky_B)) deallocate(cholesky_B)
            cholesky_naux = 0
            cholesky_ready = .false.
            cholesky_omp_reported = .false.
            cholesky_block_omp_reported = .false.
        case default
        end select
    end subroutine clear_fock_backend

end module fock_builder_module
