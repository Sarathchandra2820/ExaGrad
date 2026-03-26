module two_electron_df_module
    use iso_c_binding
    use one_eints
    use libcint_interface
    use math_utils
    implicit none

    logical :: true_df_ready = .false.
    integer :: true_df_naux  = 0
    real(c_double), allocatable :: true_df_B(:,:,:)   ! (nao, nao, true_df_naux)

contains

    pure integer function pair_index(mu, nu) result(p)
        implicit none
        integer, intent(in) :: mu, nu
        p = mu * (mu - 1) / 2 + nu
    end function pair_index

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
        read(61, *, iostat=ios) natm_mol;    if (ios /= 0) stop 'Failed reading natm_mol'
        read(61, *, iostat=ios) nbas_mol;    if (ios /= 0) stop 'Failed reading nbas_mol'
        read(61, *, iostat=ios) natm_aux;    if (ios /= 0) stop 'Failed reading natm_aux'
        read(61, *, iostat=ios) nbas_aux;    if (ios /= 0) stop 'Failed reading nbas_aux'
        read(61, *, iostat=ios) natm_all;    if (ios /= 0) stop 'Failed reading natm_all'
        read(61, *, iostat=ios) nbas_all;    if (ios /= 0) stop 'Failed reading nbas_all'
        read(61, *, iostat=ios) env_all_size; if (ios /= 0) stop 'Failed reading env_all_size'
        read(61, *, iostat=ios) nao_mol;     if (ios /= 0) stop 'Failed reading nao_mol'
        read(61, *, iostat=ios) nao_aux;     if (ios /= 0) stop 'Failed reading nao_aux'
        close(61)

        if (nao_mol /= nao)   stop 'True DF metadata nao_mol does not match current molecule'
        if (nbas_mol /= nbas) stop 'True DF metadata nbas_mol does not match current molecule'

        true_df_naux = nao_aux
        if (true_df_naux <= 0) stop 'True DF metadata has non-positive auxiliary AO count'

        allocate(atm_all(natm_all*6), bas_all(nbas_all*8), env_all(env_all_size))

        open(unit=62, file='ints/df_meta_atm.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open ints/df_meta_atm.txt'
        read(62, *, iostat=ios) atm_all; if (ios /= 0) stop 'Failed reading ints/df_meta_atm.txt'
        close(62)

        open(unit=63, file='ints/df_meta_bas.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open ints/df_meta_bas.txt'
        read(63, *, iostat=ios) bas_all; if (ios /= 0) stop 'Failed reading ints/df_meta_bas.txt'
        close(63)

        open(unit=64, file='ints/df_meta_env.txt', status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Could not open ints/df_meta_env.txt'
        read(64, *, iostat=ios) env_all; if (ios /= 0) stop 'Failed reading ints/df_meta_env.txt'
        close(64)

        allocate(ao_loc_all(nbas_all+1))
        ao_loc_all(1) = 1
        do i = 1, nbas_all
            ai = bas_all((i-1)*8 + 2)
            aj = bas_all((i-1)*8 + 4)
            ao_loc_all(i+1) = ao_loc_all(i) + (2*ai + 1) * aj
        end do

        aux_ao0 = ao_loc_all(nbas_mol+1) - 1
        npair   = nao * (nao + 1) / 2

        allocate(metric(true_df_naux,true_df_naux), rhs(true_df_naux,npair))
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

        allocate(buf2c(max_aux_shell*max_aux_shell))
        allocate(buf3c(max_mol_shell*max_mol_shell*max_aux_shell))

        ! 2-centre auxiliary integrals (P|Q)
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

        ! 3-centre integrals (mu nu | P)
        do shI = 1, nbas_mol
            di = ao_loc_all(shI+1) - ao_loc_all(shI)
            do shJ = 1, nbas_mol
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
                                p = pair_index(max(i,j), min(i,j))
                                rhs(ip,p) = buf3c(((ap-1)*dj + aj-1)*di + ai)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        ! Apply inverse metric: B = L^{-1} * (mu nu | P)  via Cholesky factorisation of (P|Q)
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

        true_df_ready = .true.
        print "(A, I8)", "  True-DF auxiliary rank =", true_df_naux
        print "(A)",     "  True-DF source: metadata + Fortran-built CDERI"

        deallocate(pair_mu, pair_nu, buf2c, buf3c, metric, rhs)
        deallocate(ao_loc_all, atm_all, bas_all, env_all)
    end subroutine initialize_true_df_factors

    subroutine clear_true_df_factors()
        implicit none
        if (allocated(true_df_B)) deallocate(true_df_B)
        true_df_naux  = 0
        true_df_ready = .false.
    end subroutine clear_true_df_factors

end module two_electron_df_module
