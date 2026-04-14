module two_electron_df_module
    use iso_c_binding
    use molecule_loader, only: ensure_molecule_loaded, nbas, nao, active_mol_dir
    use libcint_interface
    use math_utils
    use jk_contraction_module, only: contract_jk
    implicit none

    logical :: true_df_ready = .false.
    logical :: true_df_packed_ready = .false.
    integer :: true_df_naux  = 0
    real(c_double), allocatable :: true_df_B(:,:,:)   ! (nao, nao, true_df_naux)
    real(c_double), allocatable :: true_df_Bp(:,:)     ! (true_df_naux, npair)

    ! Persistent workspace for contract_jk_true_df_blocked; allocated once,
    ! freed by clear_true_df_factors.
    integer :: jk_bs = 0
    real(c_double), allocatable :: jk_Bblk(:,:,:)  ! (nao, nao, jk_bs)
    real(c_double), allocatable :: jk_Jblk(:,:)    ! (nao, nao)
    real(c_double), allocatable :: jk_Kblk(:,:)    ! (nao, nao)

contains

    pure integer function pair_index(mu, nu) result(p)
        implicit none
        integer, intent(in) :: mu, nu
        p = mu * (mu - 1) / 2 + nu
    end function pair_index

    ! Read DF metadata, evaluate 2c and 3c integrals, and Cholesky-solve for CDERIs.
    ! Returns rhs(naux, npair) = L^{-1} * (mu nu | P) and the auxiliary basis count naux.
    subroutine load_df_cderi(rhs, naux)
        implicit none
        real(c_double), allocatable, intent(out) :: rhs(:,:)
        integer, intent(out) :: naux

        integer :: ios, i, j, ai, aj, ap
        integer :: shI, shJ, shP, shQ
        integer :: di, dj, dp, dq, npair
        integer :: nbas_mol, nbas_aux, nbas_all
        integer :: natm_mol, natm_aux, natm_all
        integer :: env_all_size, nao_mol, nao_aux
        integer :: max_mol_shell, max_aux_shell
        integer :: aux_ao0, ip, info, unit_id
        integer(c_int) :: shls2(2), shls3(3), cdims2(2), cdims3(3)
        integer(c_int), allocatable :: atm_all(:), bas_all(:)
        integer, allocatable :: ao_loc_all(:)
        real(c_double), allocatable :: env_all(:), metric(:,:), buf2c(:), buf3c(:)
        character(len=2048) :: ints_dir

        call ensure_molecule_loaded()

        ints_dir = trim(active_mol_dir)//'/ints'

        open(newunit=unit_id, file=trim(ints_dir)//'/df_meta_info.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open df_meta_info.txt in molecule directory'
        read(unit_id, *, iostat=ios) natm_mol;    if (ios /= 0) stop 'Failed reading natm_mol'
        read(unit_id, *, iostat=ios) nbas_mol;    if (ios /= 0) stop 'Failed reading nbas_mol'
        read(unit_id, *, iostat=ios) natm_aux;    if (ios /= 0) stop 'Failed reading natm_aux'
        read(unit_id, *, iostat=ios) nbas_aux;    if (ios /= 0) stop 'Failed reading nbas_aux'
        read(unit_id, *, iostat=ios) natm_all;    if (ios /= 0) stop 'Failed reading natm_all'
        read(unit_id, *, iostat=ios) nbas_all;    if (ios /= 0) stop 'Failed reading nbas_all'
        read(unit_id, *, iostat=ios) env_all_size; if (ios /= 0) stop 'Failed reading env_all_size'
        read(unit_id, *, iostat=ios) nao_mol;     if (ios /= 0) stop 'Failed reading nao_mol'
        read(unit_id, *, iostat=ios) nao_aux;     if (ios /= 0) stop 'Failed reading nao_aux'
        close(unit_id)

        if (nao_mol /= nao)   stop 'True DF metadata nao_mol does not match current molecule'
        if (nbas_mol /= nbas) stop 'True DF metadata nbas_mol does not match current molecule'

        naux = nao_aux
        if (naux <= 0) stop 'True DF metadata has non-positive auxiliary AO count'

        allocate(atm_all(natm_all*6), bas_all(nbas_all*8), env_all(env_all_size))

        open(newunit=unit_id, file=trim(ints_dir)//'/df_meta_atm.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open df_meta_atm.txt in molecule directory'
        read(unit_id, *, iostat=ios) atm_all; if (ios /= 0) stop 'Failed reading df_meta_atm.txt'
        close(unit_id)

        open(newunit=unit_id, file=trim(ints_dir)//'/df_meta_bas.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open df_meta_bas.txt in molecule directory'
        read(unit_id, *, iostat=ios) bas_all; if (ios /= 0) stop 'Failed reading df_meta_bas.txt'
        close(unit_id)

        open(newunit=unit_id, file=trim(ints_dir)//'/df_meta_env.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open df_meta_env.txt in molecule directory'
        read(unit_id, *, iostat=ios) env_all; if (ios /= 0) stop 'Failed reading df_meta_env.txt'
        close(unit_id)

        allocate(ao_loc_all(nbas_all+1))
        ao_loc_all(1) = 1
        do i = 1, nbas_all
            ai = bas_all((i-1)*8 + 2)
            aj = bas_all((i-1)*8 + 4)
            ao_loc_all(i+1) = ao_loc_all(i) + (2*ai + 1) * aj
        end do

        aux_ao0 = ao_loc_all(nbas_mol+1) - 1
        npair   = nao * (nao + 1) / 2

        allocate(metric(naux, naux), rhs(naux, npair))
        metric = 0.0d0
        rhs    = 0.0d0

        max_mol_shell = 0
        do shI = 1, nbas_mol
            max_mol_shell = max(max_mol_shell, ao_loc_all(shI+1) - ao_loc_all(shI))
        end do
        max_aux_shell = 0
        do shP = 1, nbas_aux
            max_aux_shell = max(max_aux_shell, ao_loc_all(nbas_mol+shP+1) - ao_loc_all(nbas_mol+shP))
        end do

        ! 2-centre auxiliary integrals (P|Q)
        ! Parallelise over shP; each thread owns its buf2c to avoid races on the
        ! libcint output buffer. Writes to metric(i,j) are race-free because
        ! distinct (shP,shQ) pairs map to distinct (i,j) index pairs.
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP          PRIVATE(buf2c, shI, dp, shQ, shJ, dq, cdims2, shls2, ap, i, aj, j, info)
        allocate(buf2c(max_aux_shell * max_aux_shell))
        !$OMP DO SCHEDULE(STATIC)
        do shP = 1, nbas_aux
            shI = nbas_mol + shP
            dp  = ao_loc_all(shI+1) - ao_loc_all(shI)
            do shQ = 1, nbas_aux
                shJ = nbas_mol + shQ
                dq  = ao_loc_all(shJ+1) - ao_loc_all(shJ)
                cdims2 = [int(dp,c_int), int(dq,c_int)]
                shls2  = [int(shI-1,c_int), int(shJ-1,c_int)]
                buf2c(1:dp*dq) = 0.0d0
                info = cint2c2e_sph(buf2c, cdims2, shls2, atm_all, int(natm_all,c_int), &
                                    bas_all, int(nbas_all,c_int), env_all, c_null_ptr, c_null_ptr)
                do ap = 1, dp
                    i = ao_loc_all(shI) + ap - 1 - aux_ao0
                    do aj = 1, dq
                        j = ao_loc_all(shJ) + aj - 1 - aux_ao0
                        metric(i,j) = buf2c((aj-1)*dp + ap)
                    end do
                end do
            end do
        end do
        !$OMP END DO
        deallocate(buf2c)
        !$OMP END PARALLEL

        ! 3-centre integrals (mu nu | P)
        ! Parallelise over shI with dynamic scheduling: the triangular inner loop
        ! (shJ = 1, shI) means shI=1 has 1 iteration while shI=nbas_mol has nbas_mol,
        ! so work is heavily skewed and dynamic load-balancing is essential.
        ! Writes to rhs(ip, p) are race-free: distinct (shI,shJ) pairs produce
        ! disjoint pair indices p, and distinct shP values produce disjoint ip rows.
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP          PRIVATE(buf3c, di, shJ, dj, shP, shQ, dp, cdims3, shls3, ap, ip, aj, j, ai, i, info)
        allocate(buf3c(max_mol_shell * max_mol_shell * max_aux_shell))
        !$OMP DO SCHEDULE(DYNAMIC)
        do shI = 1, nbas_mol
            di = ao_loc_all(shI+1) - ao_loc_all(shI)
            do shJ = 1, shI
                dj = ao_loc_all(shJ+1) - ao_loc_all(shJ)
                do shP = 1, nbas_aux
                    shQ = nbas_mol + shP
                    dp  = ao_loc_all(shQ+1) - ao_loc_all(shQ)
                    cdims3 = [int(di,c_int), int(dj,c_int), int(dp,c_int)]
                    shls3  = [int(shI-1,c_int), int(shJ-1,c_int), int(shQ-1,c_int)]
                    buf3c(1:di*dj*dp) = 0.0d0
                    info = cint3c2e_sph(buf3c, cdims3, shls3, atm_all, int(natm_all,c_int), &
                                        bas_all, int(nbas_all,c_int), env_all, c_null_ptr, c_null_ptr)
                    do ap = 1, dp
                        ip = ao_loc_all(shQ) + ap - 1 - aux_ao0
                        do aj = 1, dj
                            j = ao_loc_all(shJ) + aj - 1
                            do ai = 1, di
                                i = ao_loc_all(shI) + ai - 1
                                rhs(ip, pair_index(max(i,j), min(i,j))) = &
                                    buf3c(((ap-1)*dj + aj-1)*di + ai)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$OMP END DO
        deallocate(buf3c)
        !$OMP END PARALLEL

  

        ! Apply inverse metric: B = L^{-1} * (mu nu | P) via Cholesky factorisation of (P|Q)
        call dpotrf('L', int(naux,c_int), metric, int(naux,c_int), info)
        if (info /= 0) stop 'dpotrf failed while factoring true DF metric'
        call dtrsm('L', 'L', 'N', 'N', int(naux,c_int), int(npair,c_int), 1.0d0, &
                   metric, int(naux,c_int), rhs, int(naux,c_int))

        deallocate(metric, ao_loc_all, atm_all, bas_all, env_all)
    end subroutine load_df_cderi

    subroutine initialize_true_df_factors()
        implicit none

        integer :: q, p, mu, nu, npair
        real(c_double), allocatable :: rhs(:,:)

        if (true_df_ready) return

        if (true_df_packed_ready) then
            npair = nao * (nao + 1) / 2
            if (allocated(true_df_B)) deallocate(true_df_B)
            allocate(true_df_B(nao, nao, true_df_naux))
            true_df_B = 0.0d0

            do q = 1, true_df_naux
                p = 0
                do mu = 1, nao
                    do nu = 1, mu
                        p = p + 1
                        true_df_B(mu,nu,q) = true_df_Bp(q,p)
                        true_df_B(nu,mu,q) = true_df_Bp(q,p)
                    end do
                end do
            end do

            true_df_ready = .true.
            print "(A, I8)", "  True-DF auxiliary rank =", true_df_naux
            print "(A)",     "  True-DF source: packed factors -> dense tensor"
            return
        end if

        call load_df_cderi(rhs, true_df_naux)

        npair = nao * (nao + 1) / 2
        if (allocated(true_df_B)) deallocate(true_df_B)
        allocate(true_df_B(nao, nao, true_df_naux))
        true_df_B = 0.0d0

        do q = 1, true_df_naux
            p = 0
            do mu = 1, nao
                do nu = 1, mu
                    p = p + 1
                    true_df_B(mu,nu,q) = rhs(q,p)
                    true_df_B(nu,mu,q) = rhs(q,p)
                end do
            end do
        end do

        true_df_ready = .true.
        if (allocated(true_df_Bp)) deallocate(true_df_Bp)
        true_df_packed_ready = .false.
        print "(A, I8)", "  True-DF auxiliary rank =", true_df_naux
        print "(A)",     "  True-DF source: metadata + Fortran-built CDERI"

        deallocate(rhs)
    end subroutine initialize_true_df_factors

    subroutine initialize_true_df_packed_factors()
        implicit none

        integer :: q, p, mu, nu, npair
        real(c_double), allocatable :: rhs(:,:)

        if (true_df_packed_ready) return

        if (true_df_ready) then
            npair = nao * (nao + 1) / 2
            if (allocated(true_df_Bp)) deallocate(true_df_Bp)
            allocate(true_df_Bp(true_df_naux, npair))
            do q = 1, true_df_naux
                p = 0
                do mu = 1, nao
                    do nu = 1, mu
                        p = p + 1
                        true_df_Bp(q,p) = true_df_B(mu,nu,q)
                    end do
                end do
            end do
            true_df_packed_ready = .true.
            print "(A, I8)", "  True-DF auxiliary rank =", true_df_naux
            print "(A)",     "  True-DF source: dense tensor -> packed factors"
            return
        end if

        call load_df_cderi(rhs, true_df_naux)

        if (allocated(true_df_Bp)) deallocate(true_df_Bp)
        allocate(true_df_Bp(true_df_naux, nao * (nao + 1) / 2))
        true_df_Bp = rhs

        if (allocated(true_df_B)) deallocate(true_df_B)
        true_df_ready = .false.
        true_df_packed_ready = .true.

        print "(A, I8)", "  True-DF auxiliary rank =", true_df_naux
        print "(A)",     "  True-DF source: metadata + packed factors"

        deallocate(rhs)
    end subroutine initialize_true_df_packed_factors

    subroutine contract_jk_true_df_blocked(P, J, K, block_size)
        implicit none
        real(c_double), intent(in)  :: P(nao,nao)
        real(c_double), intent(out) :: J(nao,nao), K(nao,nao)
        integer, intent(in), optional :: block_size

        integer :: bs, npair, q0, q1, qb, q, qq, mu, nu, ipair

        call initialize_true_df_packed_factors()

       

        if (.not. true_df_packed_ready) stop 'Blocked True-DF initialisation failed'

        

        bs = 16
        if (present(block_size)) bs = max(1, block_size)
        npair = nao * (nao + 1) / 2

        ! Lazy allocation: only (re)allocate when block size changes.
        if (.not. allocated(jk_Bblk) .or. jk_bs /= bs) then
            if (allocated(jk_Bblk)) deallocate(jk_Bblk)
            allocate(jk_Bblk(nao, nao, bs))
            jk_bs = bs
        end if
        if (.not. allocated(jk_Jblk)) allocate(jk_Jblk(nao, nao))
        if (.not. allocated(jk_Kblk)) allocate(jk_Kblk(nao, nao))

        J = 0.0d0
        K = 0.0d0

        do q0 = 1, true_df_naux, bs
            q1 = min(true_df_naux, q0 + bs - 1)
            qb = q1 - q0 + 1

            jk_Bblk(:,:,1:qb) = 0.0d0
            do q = 1, qb
                qq = q0 + q - 1
                ipair = 0
                do mu = 1, nao
                    do nu = 1, mu
                        ipair = ipair + 1
                        jk_Bblk(mu,nu,q) = true_df_Bp(qq,ipair)
                        jk_Bblk(nu,mu,q) = true_df_Bp(qq,ipair)
                    end do
                end do
            end do
            
            

            call contract_jk(jk_Bblk(:,:,1:qb), qb, P, jk_Jblk, jk_Kblk)
            
            J = J + jk_Jblk
            K = K + jk_Kblk
        end do
    end subroutine contract_jk_true_df_blocked

    subroutine clear_true_df_factors()
        implicit none
        if (allocated(true_df_B))    deallocate(true_df_B)
        if (allocated(true_df_Bp))   deallocate(true_df_Bp)
        if (allocated(jk_Bblk))      deallocate(jk_Bblk)
        if (allocated(jk_Jblk))      deallocate(jk_Jblk)
        if (allocated(jk_Kblk))      deallocate(jk_Kblk)
        true_df_naux  = 0
        jk_bs         = 0
        true_df_ready = .false.
        true_df_packed_ready = .false.
    end subroutine clear_true_df_factors

end module two_electron_df_module
