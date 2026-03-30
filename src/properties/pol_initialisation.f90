module polarisability_init
    use math_utils
    use iso_c_binding
    use davidson_module
    use jk_contraction_module
    use two_electron_df_module, only: true_df_B, true_df_naux
    implicit none



contains

    subroutine transform_dipole_integrals(mol, C_mo, nocc, dip_ao, dip_mo)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double), intent(in)  :: C_mo(:,:)
        integer, intent(in) :: nocc
        real(c_double), intent(out) :: dip_ao(:,:,:)
        real(c_double), intent(out) :: dip_mo(:,:,:)

        integer :: shI, shJ, di, dj
        integer :: i, j, k, nbuf, max_ao_per_shell, n, nbas_local
        integer(c_int) :: stat, shls(2), cdims(2)
        real(c_double), allocatable :: buf(:), dip_mo_k(:,:)

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
            call transform_sigma_ao_mo(dip_ao(:,:,k), C_mo, nocc, dip_mo_k)
            dip_mo(:,:,k) = dip_mo_k
        end do

        deallocate(buf)
        if (allocated(dip_mo_k)) deallocate(dip_mo_k)
    end subroutine transform_dipole_integrals

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