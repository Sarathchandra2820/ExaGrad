module cpks
    use iso_c_binding, only: c_double, c_int, c_null_ptr
    use libcint_interface, only: cint1e_r_sph
    use one_eints, only: nao, ao_loc, atm, natm, bas, nbas, env
    use math_utils
    implicit none

contains

    subroutine transform_dipole_integrals(C_mo, dip_ao, dip_mo)
        implicit none
        real(c_double), intent(in)  :: C_mo(nao,nao)
        real(c_double), intent(out) :: dip_ao(nao,nao,3)
        real(c_double), intent(out) :: dip_mo(3,nao,nao)

        integer :: shI, shJ, di, dj
        integer :: i, j, k, nbuf, max_ao_per_shell
        integer(c_int) :: stat, shls(2), cdims(2)
        real(c_double), allocatable :: buf(:), tmp(:,:), mo_blk(:,:)

        dip_ao = 0.0d0
        dip_mo = 0.0d0

        max_ao_per_shell = 0
        do shI = 1, nbas
            max_ao_per_shell = max(max_ao_per_shell, ao_loc(shI+1) - ao_loc(shI))
        end do

        allocate(buf(3*max_ao_per_shell*max_ao_per_shell))
        allocate(tmp(nao,nao), mo_blk(nao,nao))

        do shI = 1, nbas
            di = ao_loc(shI+1) - ao_loc(shI)
            do shJ = 1, nbas
                dj = ao_loc(shJ+1) - ao_loc(shJ)

                buf(1:3*di*dj) = 0.0d0
                shls  = [int(shI-1,c_int), int(shJ-1,c_int)]
                cdims = [int(di,c_int), int(dj,c_int)]

                stat = cint1e_r_sph(buf, cdims, shls, atm, int(natm,c_int), bas, int(nbas,c_int), &
                                   env, c_null_ptr, c_null_ptr)

                if (stat /= 0_c_int) then
                    do k = 1, 3
                        do j = 1, dj
                            do i = 1, di
                                nbuf = ((k-1)*dj + (j-1))*di + i
                                dip_ao(ao_loc(shI)+i-1, ao_loc(shJ)+j-1, k) = buf(nbuf)
                            end do
                        end do
                    end do
                end if
            end do
        end do

        do k = 1, 3
            call dgemm('N','N', int(nao,c_int), int(nao,c_int), int(nao,c_int), 1.0d0, &
                       dip_ao(:,:,k), int(nao,c_int), C_mo, int(nao,c_int), 0.0d0, tmp, int(nao,c_int))
            call dgemm('T','N', int(nao,c_int), int(nao,c_int), int(nao,c_int), 1.0d0, &
                       C_mo, int(nao,c_int), tmp, int(nao,c_int), 0.0d0, mo_blk, int(nao,c_int))
            dip_mo(k,:,:) = mo_blk(:,:)
        end do

        deallocate(buf, tmp, mo_blk)
    end subroutine transform_dipole_integrals

    subroutine prepare_initial_trial_density(D, U, C_mo, nocc, nvirt)
        implicit none
        integer, intent(in) :: nocc, nvirt
        real(c_double), intent(in)  :: U(nvirt,nocc,3), C_mo(nao,nao)
        real(c_double), intent(out) :: D(nao,nao,3)

        integer :: k
        real(c_double), allocatable :: C_occ(:,:), C_vir(:,:), X(:,:)

        allocate(C_occ(nao,nocc), C_vir(nao,nvirt), X(nao,nocc))
        C_occ = C_mo(:,1:nocc)
        C_vir = C_mo(:,nocc+1:nocc+nvirt)

        D = 0.0d0
        do k = 1, 3
            call dgemm('N','N', int(nao,c_int), int(nocc,c_int), int(nvirt,c_int), 1.0d0, &
                       C_vir, int(nao,c_int), U(:,:,k), int(nvirt,c_int), 0.0d0, X, int(nao,c_int))

            call dgemm('N','T', int(nao,c_int), int(nao,c_int), int(nocc,c_int), 1.0d0, &
                       X, int(nao,c_int), C_occ, int(nao,c_int), 0.0d0, D(:,:,k), int(nao,c_int))
            call dgemm('N','T', int(nao,c_int), int(nao,c_int), int(nocc,c_int), 1.0d0, &
                       C_occ, int(nao,c_int), X, int(nao,c_int), 1.0d0, D(:,:,k), int(nao,c_int))
        end do

        deallocate(C_occ, C_vir, X)
    end subroutine prepare_initial_trial_density

    subroutine build_trial_vector(C_mo, C_mo_occ, trial_vector)
        implicit none
        real(c_double), intent(in) :: C_mo(:,:), C_mo_occ(:)
        real(c_double), intent(out) :: trial_vector(:)
        print *, "CPKS computation placeholder"
    end subroutine build_trial_vector

end module cpks