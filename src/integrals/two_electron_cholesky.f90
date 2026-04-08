module two_electron_cholesky_module
    use iso_c_binding
    use omp_lib, only: omp_get_max_threads
    use molecule_loader, only: ensure_molecule_loaded, natm, nbas, nao, atm, bas, env, ao_loc
    use libcint_interface
    use math_utils
    implicit none

    logical :: cholesky_ready = .false.
    integer :: cholesky_naux  = 0
    real(c_double), allocatable :: cholesky_B(:,:,:)   ! (nao, nao, cholesky_naux)

contains

    pure integer function pair_index(mu, nu) result(p)
        implicit none
        integer, intent(in) :: mu, nu
        p = mu * (mu - 1) / 2 + nu
    end function pair_index

    subroutine initialize_cholesky_factors(chol_tol, max_rank)
        implicit none
        real(c_double), intent(in), optional :: chol_tol
        integer,        intent(in), optional :: max_rank

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
        call ensure_molecule_loaded()

        tol      = 1.0d-10
        if (present(chol_tol)) tol = chol_tol

        npair    = nao * (nao + 1) / 2
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

        print "(A, I4, A)", "  Building Cholesky metric (", omp_get_max_threads(), " threads) ..."

        ! ERI permutation symmetry guarantees concurrent writes of identical values —
        ! last write wins safely, no atomics needed.
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

        ! Pivoted Cholesky decomposition
        allocate(L(npair, rank_cap), resid(npair), work(npair))
        L = 0.0d0
        do i = 1, npair
            resid(i) = max(0.0d0, metric(i,i))
        end do
        cholesky_naux = 0

        do k = 1, rank_cap
            piv       = maxloc(resid, 1)
            pivot_val = resid(piv)
            if (pivot_val < tol) exit

            work = metric(:,piv)
            if (k > 1) call dgemv('N', npair, k-1, -1.0d0, L, npair, L(piv,1), npair, 1.0d0, work, 1)
            L(:,k) = work / sqrt(pivot_val)
            resid  = resid - L(:,k) * L(:,k)
            where (resid < 0.0d0) resid = 0.0d0
            cholesky_naux = k
        end do

        if (allocated(cholesky_B)) deallocate(cholesky_B)
        allocate(cholesky_B(nao, nao, cholesky_naux))
        cholesky_B = 0.0d0

        do k = 1, cholesky_naux
            do p = 1, npair
                mi = pair_mu(p)
                mj = pair_nu(p)
                g  = L(p,k)
                cholesky_B(mi,mj,k) = g
                cholesky_B(mj,mi,k) = g
            end do
        end do

        cholesky_ready = .true.
        print "(A, I8)", "  Cholesky auxiliary rank =", cholesky_naux
        deallocate(metric, L, resid, work, pair_mu, pair_nu)
    end subroutine initialize_cholesky_factors

    subroutine clear_cholesky_factors()
        implicit none
        if (allocated(cholesky_B)) deallocate(cholesky_B)
        cholesky_naux  = 0
        cholesky_ready = .false.
    end subroutine clear_cholesky_factors

end module two_electron_cholesky_module
