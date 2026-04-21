module polarisability_init
    use math_utils
    use iso_c_binding
    use molecule_t, only: molecule
    use one_eints, only: activate_molecule
    use molecule_loader, only: atm, bas, env
    use sigma_build, only: transform_sigma_mo_ao, transform_sigma_ao_mo
    use jk_contraction_module
    use two_electron_df_module, only: true_df_B, true_df_naux
    implicit none



contains

    subroutine compute_dipole_ao(mol, dip_ao)
        implicit none
        type(molecule), intent(in)  :: mol
        real(c_double), intent(out) :: dip_ao(:,:,:)

        integer :: shI, shJ, di, dj
        integer :: i, j, k, nbuf, max_ao_per_shell, nbas_local
        integer(c_int) :: stat, shls(2), cdims(2)
        real(c_double), allocatable :: buf(:)

        call activate_molecule(mol)
        nbas_local = int(mol%basis%nbas)

        dip_ao = 0.0d0

        max_ao_per_shell = 0
        do shI = 1, nbas_local
            max_ao_per_shell = max(max_ao_per_shell, int(mol%basis%ao_loc(shI+1) - mol%basis%ao_loc(shI)))
        end do

        allocate(buf(3*max_ao_per_shell*max_ao_per_shell))

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

        deallocate(buf)

    end subroutine compute_dipole_ao

    ! ================================================================
    ! Transform dipole integrals from AO to MO basis using C_mo:
    subroutine transform_dipole_integrals(mol, C_mo, nocc, dip_ao, dip_mo)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double), intent(in) :: C_mo(:,:)
        integer, intent(in) :: nocc
        real(c_double), intent(in) :: dip_ao(:,:,:)
        real(c_double), allocatable, intent(out) :: dip_mo(:,:,:)

        integer :: nao, nvir, k
        real(c_double), allocatable :: temp(:,:)

        nao = size(C_mo, 1)
        nvir = size(C_mo, 2) - nocc

        allocate(dip_mo(nvir, nocc, 3))
        allocate(temp(nao, max(nocc, nvir)))

        do k = 1, 3
            ! temp(:,1:nocc) = dip_ao(:,:,k) @ C_mo(:,1:nocc)
            call dgemm('N', 'N', nao, nocc, nao, 1.0d0, dip_ao(:,:,k), nao, &
                       C_mo(:, 1:nocc), nao, 0.0d0, temp(:, 1:nocc), nao)

            ! dip_mo(:,:,k) = C_mo(:,nocc+1:nocc+nvir)^T @ temp(:,1:nocc)
            call dgemm('T', 'N', nvir, nocc, nao, 1.0d0, C_mo(:, nocc+1:nocc+nvir), nao, &
                       temp(:, 1:nocc), nao, 0.0d0, dip_mo(:,:,k), nvir)
        end do

        deallocate(temp)

    end subroutine transform_dipole_integrals

    ! ================================================================
    ! transform_dipole_ao_to_lmo
    !
    ! Transforms AO dipole integrals directly into the LMO basis.
    ! Forms C_lmo = C_mo_block @ U, then contracts:
    !   dip_lmo(:,:,k) = C_lmo^T @ dip_ao(:,:,k) @ C_lmo
    !
    ! Returns three blocks:
    !   dip_oo(nocc, nocc, 3)  -- occupied-occupied
    !   dip_vv(nvir, nvir, 3)  -- virtual-virtual
    !   dip_ov(nocc, nvir, 3)  -- occupied-virtual
    ! ================================================================
    subroutine transform_dipole_ao_to_lmo(C_lmo_occ, C_lmo_vir, nocc, nvir, dip_ao, &
                                          dip_oo, dip_vv, dip_ov)
        implicit none
        real(c_double), intent(in)               :: C_lmo_occ(:,:), C_lmo_vir(:,:)
        integer,        intent(in)               :: nocc, nvir
        real(c_double), intent(in)               :: dip_ao(:,:,:)
        real(c_double), allocatable, intent(out) :: dip_oo(:,:,:)
        real(c_double), allocatable, intent(out) :: dip_vv(:,:,:)
        real(c_double), allocatable, intent(out) :: dip_ov(:,:,:)

        real(c_double), allocatable :: temp(:,:)
        integer :: nao, k

        nao = size(C_lmo_occ, 1)

        allocate(dip_oo(nocc, nocc, 3))
        allocate(dip_vv(nvir, nvir, 3))
        allocate(dip_ov(nocc, nvir, 3))
        allocate(temp(nao, max(nocc, nvir)))

        do k = 1, 3
            ! occ-occ: C_lmo_occ^T @ dip_ao @ C_lmo_occ
            call dgemm('N', 'N', nao, nocc, nao, 1.0d0, dip_ao(:,:,k), nao, &
                       C_lmo_occ, nao, 0.0d0, temp(:, 1:nocc), nao)
            call dgemm('T', 'N', nocc, nocc, nao, 1.0d0, C_lmo_occ, nao, &
                       temp(:, 1:nocc), nao, 0.0d0, dip_oo(:,:,k), nocc)

            ! vir-vir: C_lmo_vir^T @ dip_ao @ C_lmo_vir
            call dgemm('N', 'N', nao, nvir, nao, 1.0d0, dip_ao(:,:,k), nao, &
                       C_lmo_vir, nao, 0.0d0, temp(:, 1:nvir), nao)
            call dgemm('T', 'N', nvir, nvir, nao, 1.0d0, C_lmo_vir, nao, &
                       temp(:, 1:nvir), nao, 0.0d0, dip_vv(:,:,k), nvir)

            ! occ-vir: C_lmo_occ^T @ dip_ao @ C_lmo_vir
            ! temp(:,1:nvir) already holds dip_ao @ C_lmo_vir from above
            call dgemm('T', 'N', nocc, nvir, nao, 1.0d0, C_lmo_occ, nao, &
                       temp(:, 1:nvir), nao, 0.0d0, dip_ov(:,:,k), nocc)
        end do

        deallocate(temp)

    end subroutine transform_dipole_ao_to_lmo

    subroutine prepare_initial_trial_density(D, U, C_mo)
        implicit none
        real(c_double), intent(in)  :: U(:,:,:), C_mo(:,:)
        real(c_double), intent(out) :: D(:,:,:)

        integer :: k
        real(c_double), allocatable :: D_k(:,:)

        D = 0.0d0
        do k = 1, 3
            call transform_sigma_mo_ao(U(:,:,k), C_mo, D_k)
            D(:,:,k) = D_k
        end do

        if (allocated(D_k)) deallocate(D_k)
    end subroutine prepare_initial_trial_density




end module polarisability_init
