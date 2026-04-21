module sigma_build

    use math_utils
    use iso_c_binding
    use jk_contraction_module
    use libcint_interface
    use two_electron_df_module, only: true_df_B, true_df_naux

    implicit none

    real(c_double), allocatable, save :: sigma_ao_scratch(:,:), J_scratch(:,:), K_scratch(:,:)
    real(c_double), allocatable, save :: mo2ao_tmp1(:,:), mo2ao_tmp2(:,:), ao2mo_tmp(:,:)
!$omp threadprivate(sigma_ao_scratch, J_scratch, K_scratch, mo2ao_tmp1, mo2ao_tmp2, ao2mo_tmp)


    contains


    subroutine ensure_sigma_scratch(nao, nocc)

        integer, intent(in) :: nao, nocc

        if (allocated(sigma_ao_scratch)) then
            if (size(sigma_ao_scratch,1) /= nao .or. size(sigma_ao_scratch,2) /= nao) then
                deallocate(sigma_ao_scratch, J_scratch, K_scratch, mo2ao_tmp1, mo2ao_tmp2, ao2mo_tmp)
            end if
        end if

        if (.not. allocated(sigma_ao_scratch)) then
            allocate(sigma_ao_scratch(nao, nao), J_scratch(nao, nao), K_scratch(nao, nao))
            allocate(mo2ao_tmp1(nao, nocc), mo2ao_tmp2(nao, nao), ao2mo_tmp(nao, nocc))
        else if (size(mo2ao_tmp1,2) /= nocc .or. size(ao2mo_tmp,2) /= nocc) then
            deallocate(mo2ao_tmp1, ao2mo_tmp)
            allocate(mo2ao_tmp1(nao, nocc), ao2mo_tmp(nao, nocc))
        end if

    end subroutine ensure_sigma_scratch


    ! Transform U_{ai} (nvir x nocc, MO basis) -> sigma_ao (nao x nao, AO basis)
    ! sigma_ao = C_vir * U * C_occ^T + (C_vir * U * C_occ^T)^T  (symmetrised)
    subroutine transform_sigma_mo_ao(sigma_mo, C_mo, sigma_ao)

        integer :: nao, nocc, nvir
        real(c_double), intent(in)  :: sigma_mo(:,:), C_mo(:,:)
        real(c_double), allocatable, intent(out) :: sigma_ao(:,:)
        real(c_double), allocatable :: temp1(:,:), temp2(:,:)
        integer :: i, j

        nao  = size(C_mo, 1)
        nocc = size(sigma_mo, 2)
        nvir = size(sigma_mo, 1)

        allocate(temp1(nao, nocc), temp2(nao, nao), sigma_ao(nao, nao))

        ! temp1(nao,nocc) = C_vir(nao,nvir) * U(nvir,nocc)
        call dgemm('N','N', nao, nocc, nvir, 1.0d0, C_mo(:,nocc+1:), nao, &
                   sigma_mo, nvir, 0.0d0, temp1, nao)
        ! temp2(nao,nao) = temp1(nao,nocc) * C_occ(nao,nocc)^T
        call dgemm('N','T', nao, nao, nocc, 1.0d0, temp1, nao, &
                   C_mo(:,1:nocc), nao, 0.0d0, temp2, nao)

        do j = 1, nao
            do i = 1, nao
                sigma_ao(i,j) = temp2(i,j) + temp2(j,i)
            end do
        end do

        deallocate(temp1, temp2)

    end subroutine transform_sigma_mo_ao


    ! Transform sigma_ao (nao x nao, AO basis) -> sigma_mo (nvir x nocc, MO basis)
    ! sigma_mo = C_vir^T * sigma_ao * C_occ
    subroutine transform_sigma_ao_mo(sigma_ao, C_mo, nocc, sigma_mo)

        integer :: nao, nvir
        integer, intent(in) :: nocc
        real(c_double), intent(in)  :: sigma_ao(:,:), C_mo(:,:)
        real(c_double), allocatable, intent(out) :: sigma_mo(:,:)
        real(c_double), allocatable :: temp1(:,:)

        nao  = size(C_mo, 1)
        nvir = size(C_mo, 2) - nocc

        allocate(temp1(nao, nocc), sigma_mo(nvir, nocc))

        ! Step 1: temp1(nao,nocc) = sigma_ao(nao,nao) * C_occ(nao,nocc)
        call dgemm('N','N', nao, nocc, nao, 1.0d0, sigma_ao, nao, &
                   C_mo(:,1:nocc), nao, 0.0d0, temp1, nao)
        ! Step 2: sigma_mo(nvir,nocc) = C_vir(nao,nvir)^T * temp1(nao,nocc)
        call dgemm('T','N', nvir, nocc, nao, 1.0d0, C_mo(:,nocc+1:), nao, &
                   temp1, nao, 0.0d0, sigma_mo, nvir)

        deallocate(temp1)

    end subroutine transform_sigma_ao_mo



    ! Apply the two-electron kernel to sigma_init (nvir x nocc, MO basis).
    ! Returns sigma_out (nvir x nocc, MO basis) = C_vir^T [J - 0.5*K][D_sym] C_occ
    ! where D_sym = C_vir*sigma_init*C_occ^T + transpose.
    ! For CPHF use: multiply result by 2 to get the full G_ai response.
    subroutine contract_sigma_over_kernel(method, sigma_init, C_mo, nocc, B, naux, sigma_out, Q)

        character(len=*), intent(in)  :: method
        integer, intent(in) :: nocc, naux
        real(c_double), intent(in) :: sigma_init(:,:)
        real(c_double), intent(in) :: C_mo(:,:)
        real(c_double), intent(in) :: B(:,:,:)
        real(c_double), intent(out) :: sigma_out(:,:)
        real(c_double), intent(in), optional :: Q(:,:)

        integer :: nao_local, nvir, i, j

        nao_local = size(C_mo, 1)
        nvir = size(sigma_init, 1)

        call ensure_sigma_scratch(nao_local, nocc)

        ! MO -> AO (reused scratch)
        call dgemm('N','N', nao_local, nocc, nvir, 1.0d0, C_mo(:,nocc+1:), nao_local, &
                   sigma_init, nvir, 0.0d0, mo2ao_tmp1, nao_local)
        call dgemm('N','T', nao_local, nao_local, nocc, 1.0d0, mo2ao_tmp1, nao_local, &
                   C_mo(:,1:nocc), nao_local, 0.0d0, mo2ao_tmp2, nao_local)

        do j = 1, nao_local
            do i = 1, nao_local
                sigma_ao_scratch(i,j) = mo2ao_tmp2(i,j) + mo2ao_tmp2(j,i)
            end do
        end do

        select case (trim(method))
        case ('true_df', 'block_cholesky', 'cholesky')
            if (naux <= 0) stop 'Error: no auxiliary rank for DF/Cholesky contraction'
            call contract_jk(B, naux, sigma_ao_scratch, J_scratch, K_scratch)
            ! J <- J - 0.5*K
            call daxpy(nao_local*nao_local, -0.5d0, K_scratch, 1, J_scratch, 1)
        case default
            if (.not. present(Q)) stop 'Error: Schwarz bounds Q required for direct contraction'
            call contract_jk_direct(sigma_ao_scratch, Q, J_scratch, K_scratch)
            ! For direct path K already carries the 0.5 factor internally
            call daxpy(nao_local*nao_local, -1.0d0, K_scratch, 1, J_scratch, 1)
        end select

        ! AO -> MO (reused scratch)
        call dgemm('N','N', nao_local, nocc, nao_local, 1.0d0, J_scratch, nao_local, &
                   C_mo(:,1:nocc), nao_local, 0.0d0, ao2mo_tmp, nao_local)
        call dgemm('T','N', nvir, nocc, nao_local, 1.0d0, C_mo(:,nocc+1:), nao_local, &
                   ao2mo_tmp, nao_local, 0.0d0, sigma_out, nvir)

    end subroutine contract_sigma_over_kernel


    ! Build the orbital energy denominator matrix: denom(a,i) = 1/(e_a - e_i)
    ! where a = virtual (nocc+1..nmo), i = occupied (1..nocc).
    subroutine build_energy_denominator(mo_energies, nocc, nvir, energy_denom)

        real(c_double), intent(in) :: mo_energies(:)
        integer, intent(in) :: nocc, nvir
        real(c_double), allocatable, intent(out) :: energy_denom(:,:)
        integer :: a, i

        allocate(energy_denom(nvir, nocc))
        do i = 1, nocc
            do a = 1, nvir
                energy_denom(a, i) = 1.0d0 / (mo_energies(nocc + a) - mo_energies(i))
            end do
        end do

    end subroutine build_energy_denominator


    ! Element-wise multiply: sigma_out(a,i) = sigma_in(a,i) * energy_denom(a,i)
    subroutine contract_wt_denominator(sigma_in, energy_denom, sigma_out)

        real(c_double), intent(in) :: sigma_in(:,:), energy_denom(:,:)
        real(c_double), allocatable, intent(inout) :: sigma_out(:,:)
        integer :: a, i, nocc, nvir

        nvir = size(sigma_in, 1)
        nocc = size(sigma_in, 2)

        if (.not. allocated(sigma_out)) allocate(sigma_out(nvir, nocc))

        do i = 1, nocc
            do a = 1, nvir
                sigma_out(a, i) = sigma_in(a, i) * energy_denom(a, i)
            end do
        end do

    end subroutine contract_wt_denominator


end module sigma_build
