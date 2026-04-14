module jk_contraction_module
    use iso_c_binding
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    use molecule_loader, only: ensure_molecule_loaded, natm, nbas, nao, atm, bas, env, ao_loc
    use libcint_interface
    use math_utils
    implicit none

    real(c_double), parameter :: SCHWARZ_TOL = 1.0d-12

contains

    ! Build J and K from DF/Cholesky auxiliary tensors B(nao,nao,naux) and density P.
    ! J and K are zeroed on entry.
    ! Caller assembles: F = Hcore + J - 0.5*K
    !
    ! Optional block_size: if provided, auxiliary index is chunked into blocks and
    ! distributed manually across threads (better for cache locality at large naux).
    subroutine contract_jk(B, naux, P, J, K, block_size)
        implicit none
        integer,        intent(in)           :: naux
        real(c_double), intent(in)           :: B(nao,nao,naux)
        real(c_double), intent(in)           :: P(nao,nao)
        real(c_double), intent(out)          :: J(nao,nao), K(nao,nao)
        integer,        intent(in), optional :: block_size

        integer :: q, q0, q1, iblk, nblocks, bs
        integer :: tid, nthreads
        real(c_double) :: zq
        real(c_double), allocatable :: T_local(:,:)

        J = 0.0d0
        K = 0.0d0

        if (present(block_size) .and. block_size > 0) then
            bs      = block_size
            nblocks = (naux + bs - 1) / bs

   
            allocate(T_local(nao,nao))
            tid      = omp_get_thread_num()
            nthreads = omp_get_num_threads()
            do iblk = tid + 1, nblocks, nthreads
                q0 = (iblk-1)*bs + 1
                q1 = min(naux, iblk*bs)
                do q = q0, q1
                    zq = sum(P * B(:,:,q))
                    J  = J + zq * B(:,:,q)
                    call dgemm('N','N', nao,nao,nao, 1.0d0, B(:,:,q),nao, P,nao, 0.0d0, T_local,nao)
                    call dgemm('N','T', nao,nao,nao, 1.0d0, T_local,nao, B(:,:,q),nao, 1.0d0, K,nao)
                end do
            end do
            deallocate(T_local)
    

        else

            allocate(T_local(nao,nao))
          
            do q = 1, naux
                ! zq = sum(P * B(:,:,q))
                ! J  = J + zq * B(:,:,q)
                zq = ddot(nao*nao, P, 1, B(1,1,q), 1)
                call daxpy(nao*nao, zq, B(1,1,q), 1, J, 1)
                call dgemm('N','N', nao,nao,nao, 1.0d0, B(:,:,q),nao, P,nao, 0.0d0, T_local,nao)
                call dgemm('N','T', nao,nao,nao, 1.0d0, T_local,nao, B(:,:,q),nao, 1.0d0, K,nao)
            end do

            deallocate(T_local)

        end if
    end subroutine contract_jk

    ! Build J and K via direct 4-index integral evaluation with Schwarz screening.
    ! Q(nbas,nbas) are the Schwarz bounds from compute_schwarz.
    ! J and K are zeroed on entry.
    ! Caller assembles: F = Hcore + J - 0.5*K
    subroutine contract_jk_direct(P, Q, J, K)
        implicit none
        real(c_double), intent(in)  :: P(nao,nao), Q(nbas,nbas)
        real(c_double), intent(out) :: J(nao,nao), K(nao,nao)

        integer :: shI, shJ, shK, shL, ai, aj, ak, al
        integer :: mi, mj, mk, ml, di, dj, dk, dl
        logical :: need_j, need_k
        integer(c_int) :: shls(4), cdims(4), info
        real(c_double) :: Pmax, g
        integer :: max_ao_per_shell, bufsize
        real(c_double), allocatable :: buf_j(:), buf_k(:)

        J = 0.0d0
        K = 0.0d0
        call ensure_molecule_loaded()

        max_ao_per_shell = 0
        do shI = 1, nbas
            max_ao_per_shell = max(max_ao_per_shell, ao_loc(shI+1) - ao_loc(shI))
        end do
        bufsize = max_ao_per_shell**4
        Pmax    = maxval(abs(P))

        !$omp parallel default(none) &
        !$omp   shared(nbas, nao, ao_loc, atm, natm, bas, env, P, Q, Pmax, bufsize) &
        !$omp   private(shI, shJ, shK, shL, di, dj, dk, dl, ai, aj, ak, al, &
        !$omp           mi, mj, mk, ml, need_j, need_k, shls, cdims, info, buf_j, buf_k, g) &
        !$omp   reduction(+: J, K)
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
                            shls  = [shI-1, shJ-1, shK-1, shL-1]
                            info  = cint2e_sph(buf_j, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                            do al = 1, dl
                                ml = ao_loc(shL) + al - 1
                                do ak = 1, dk
                                    mk = ao_loc(shK) + ak - 1
                                    do aj = 1, dj
                                        mj = ao_loc(shJ) + aj - 1
                                        do ai = 1, di
                                            mi = ao_loc(shI) + ai - 1
                                            g = buf_j((((al-1)*dk + ak-1)*dj + aj-1)*di + ai)
                                            J(mi,mj) = J(mi,mj) + P(mk,ml) * g
                                        end do
                                    end do
                                end do
                            end do
                        end if

                        if (need_k) then
                            buf_k(1:di*dk*dj*dl) = 0.0d0
                            cdims = [int(di,c_int), int(dk,c_int), int(dj,c_int), int(dl,c_int)]
                            shls  = [shI-1, shK-1, shJ-1, shL-1]
                            info  = cint2e_sph(buf_k, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                            do al = 1, dl
                                ml = ao_loc(shL) + al - 1
                                do aj = 1, dj
                                    mj = ao_loc(shJ) + aj - 1
                                    do ak = 1, dk
                                        mk = ao_loc(shK) + ak - 1
                                        do ai = 1, di
                                            mi = ao_loc(shI) + ai - 1
                                            g = buf_k((((al-1)*dj + aj-1)*dk + ak-1)*di + ai)
                                            K(mi,mj) = K(mi,mj) + 0.5d0 * P(mk,ml) * g
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
    end subroutine contract_jk_direct

    ! Compute Schwarz screening bounds Q(shI,shJ) = sqrt( max|(shI shJ|shI shJ)| )
    subroutine compute_schwarz(Q)
        implicit none
        real(c_double), intent(out) :: Q(nbas,nbas)

        integer :: shI, shJ, di, dj, max_ao
        integer(c_int) :: shls(4), cdims(4), info
        real(c_double), allocatable :: buf(:)

        Q = 0.0d0
        call ensure_molecule_loaded()
        max_ao = maxval(ao_loc(2:nbas+1) - ao_loc(1:nbas))
        allocate(buf(max_ao**4))
        do shI = 1, nbas
            di = ao_loc(shI+1) - ao_loc(shI)
            do shJ = 1, shI
                dj = ao_loc(shJ+1) - ao_loc(shJ)
                buf(1:di*dj*di*dj) = 0.0d0
                cdims = [int(di,c_int), int(dj,c_int), int(di,c_int), int(dj,c_int)]
                shls  = [shI-1, shJ-1, shI-1, shJ-1]
                info  = cint2e_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                Q(shI,shJ) = sqrt(maxval(abs(buf(1:di*dj*di*dj))))
                Q(shJ,shI) = Q(shI,shJ)
            end do
        end do
        deallocate(buf)
    end subroutine compute_schwarz

end module jk_contraction_module
