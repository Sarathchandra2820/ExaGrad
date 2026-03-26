module cpks
    use iso_c_binding, only: c_double, c_int, c_null_ptr
    use libcint_interface, only: cint1e_r_sph
    use one_eints, only: activate_molecule, atm, bas, env
    use molecule_t, only: molecule
    use math_utils
    implicit none

contains

    subroutine transform_dipole_integrals(mol, C_mo, dip_ao, dip_mo)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double), intent(in)  :: C_mo(:,:)
        real(c_double), intent(out) :: dip_ao(:,:,:)
        real(c_double), intent(out) :: dip_mo(:,:,:)

        integer :: shI, shJ, di, dj
        integer :: i, j, k, nbuf, max_ao_per_shell, n, nbas_local
        integer(c_int) :: stat, shls(2), cdims(2)
        real(c_double), allocatable :: buf(:), tmp(:,:), mo_blk(:,:)

        call activate_molecule(mol)
        n = size(C_mo, 1)
        nbas_local = int(mol%basis%nbas)

        dip_ao = 0.0d0
        dip_mo = 0.0d0

        max_ao_per_shell = 0
        do shI = 1, nbas_local
            max_ao_per_shell = max(max_ao_per_shell, int(mol%basis%ao_loc(shI+1) - mol%basis%ao_loc(shI)))
        end do

        allocate(buf(3*max_ao_per_shell*max_ao_per_shell))
        allocate(tmp(n,n), mo_blk(n,n))

        do shI = 1, nbas_local
            di = int(mol%basis%ao_loc(shI+1) - mol%basis%ao_loc(shI))
            do shJ = 1, nbas_local
                dj = int(mol%basis%ao_loc(shJ+1) - mol%basis%ao_loc(shJ))

                buf(1:3*di*dj) = 0.0d0
                shls  = [int(shI-1,c_int), int(shJ-1,c_int)]
                cdims = [int(di,c_int), int(dj,c_int)]

                stat = cint1e_r_sph(buf, cdims, shls, atm, mol%basis%natm, bas, mol%basis%nbas, &
                                   env, c_null_ptr, c_null_ptr)

                if (stat /= 0_c_int) then
                    do k = 1, 3
                        do j = 1, dj
                            do i = 1, di
                                nbuf = ((k-1)*dj + (j-1))*di + i
                                dip_ao(int(mol%basis%ao_loc(shI))+i-1, int(mol%basis%ao_loc(shJ))+j-1, k) = buf(nbuf)
                            end do
                        end do
                    end do
                end if
            end do
        end do

        do k = 1, 3
            call dgemm('N','N', int(n,c_int), int(n,c_int), int(n,c_int), 1.0d0, &
                       dip_ao(:,:,k), int(n,c_int), C_mo, int(n,c_int), 0.0d0, tmp, int(n,c_int))
            call dgemm('T','N', int(n,c_int), int(n,c_int), int(n,c_int), 1.0d0, &
                       C_mo, int(n,c_int), tmp, int(n,c_int), 0.0d0, mo_blk, int(n,c_int))
            dip_mo(k,:,:) = mo_blk(:,:)
        end do

        deallocate(buf, tmp, mo_blk)
    end subroutine transform_dipole_integrals

    subroutine prepare_initial_trial_density(D, U, C_mo, nocc, nvirt)
        implicit none
        integer, intent(in) :: nocc, nvirt
        real(c_double), intent(in)  :: U(:,:,:), C_mo(:,:)
        real(c_double), intent(out) :: D(:,:,:)

        integer :: k, n
        real(c_double), allocatable :: C_occ(:,:), C_vir(:,:), X(:,:)

        n = size(C_mo, 1)
        allocate(C_occ(n,nocc), C_vir(n,nvirt), X(n,nocc))
        C_occ = C_mo(:,1:nocc)
        C_vir = C_mo(:,nocc+1:nocc+nvirt)

        D = 0.0d0
        do k = 1, 3
            call dgemm('N','N', int(n,c_int), int(nocc,c_int), int(nvirt,c_int), 1.0d0, &
                       C_vir, int(n,c_int), U(:,:,k), int(nvirt,c_int), 0.0d0, X, int(n,c_int))

            call dgemm('N','T', int(n,c_int), int(n,c_int), int(nocc,c_int), 1.0d0, &
                       X, int(n,c_int), C_occ, int(n,c_int), 0.0d0, D(:,:,k), int(n,c_int))
            call dgemm('N','T', int(n,c_int), int(n,c_int), int(nocc,c_int), 1.0d0, &
                       C_occ, int(n,c_int), X, int(n,c_int), 1.0d0, D(:,:,k), int(n,c_int))
        end do

        deallocate(C_occ, C_vir, X)
    end subroutine prepare_initial_trial_density

    subroutine build_trial_vector(C_mo, C_mo_occ, trial_vector)
        implicit none
        real(c_double), intent(in) :: C_mo(:,:), C_mo_occ(:)
        real(c_double), intent(out) :: trial_vector(:)
        trial_vector = 0.0d0
        if (size(trial_vector) > 0 .and. size(C_mo, 1) > 0 .and. size(C_mo_occ) > 0) then
            trial_vector(1) = trial_vector(1) + 0.0d0 * (C_mo(1,1) + C_mo_occ(1))
        end if
        print *, "CPKS computation placeholder"
    end subroutine build_trial_vector

end module cpks
