module davidon

    use math_utils
    use iso_c_binding
    use jk_contraction_module
    use libcint_interface
    use two_electron_df_module, only: true_df_B, true_df_naux


    implicit none

    contains

    ! The routine for the cpks is very simple this is merely a sketch of what i want.
    ! We start with the U_{ai} (MO) -> (AO)


    subroutine transform_sigma_mo_ao(sigma_mo, C_mo, sigma_ao)

        integer :: nao,  nocc, nvir
        real(c_double), intent(in) :: sigma_mo(:,:), C_mo(:,:)
        real(c_double), allocatable, intent(out) :: sigma_ao(:,:)
        real(c_double), allocatable :: temp1(:,:), temp2(:,:)


        ! Now we need to transform U_ai from MO to AO basis
        ! U_ai is nvir x nocc, we want to get sigma_ij in AO basis
        ! sigma_ij = sum_{a,i} C_ia * U_ai
        nao = size(C_mo, 1)
        nocc = size(sigma_mo, 2)
        nvir = size(sigma_mo, 1)

        allocate(temp1(nao,nocc))
        allocate(temp2(nao,nao))
        ! First transform U_ai to AO basis

        ! temp = C_virt * U_ai * C_occ^T
        ! sigma = temp + temp^T (since sigma is symmetric)

        !calculate temp = C_vir * U_ai
        call dgemm('N','N', nao,nocc,nvir, 1.0d0, C_mo(1, nocc+1), nao, sigma_mo, nvir, 0.0d0, temp1, nao)
        call dgemm('N','T', nao,nao,nocc, 1.0d0, temp1, nao, C_mo(1,1), nao, 0.0d0, temp2, nao)
        sigma_ao = temp2 + transpose(temp2)

        deallocate(temp1, temp2)


    end subroutine transform_sigma_mo_ao

    subroutine transform_sigma_ao_mo(sigma_ao, C_mo, nocc, sigma_mo)
        integer :: nao, nvir
        integer, intent(in) :: nocc
        real(c_double), intent(in) :: sigma_ao(:,:), C_mo(:,:)
        real(c_double), allocatable, intent(out) :: sigma_mo(:,:)

        real(c_double), allocatable :: temp1(:,:)

        nao = size(C_mo, 1)
        nvir = size(C_mo, 2) - nocc

        allocate(temp1(nao,nocc))
        allocate(sigma_mo(nvir,nocc))

        ! Transform sigma from AO to MO basis
        ! sigma_ai = C_vir^T * sigma_ao * C_occ
        ! Step 1: temp(nao,nocc) = sigma_ao * C_occ
        call dgemm('N','N', nao, nocc, nao, 1.0d0, sigma_ao, nao, C_mo(1,1), nao, 0.0d0, temp1, nao)
        ! Step 2: sigma_mo(nvir,nocc) = C_vir^T * temp
        call dgemm('T','N', nvir, nocc, nao, 1.0d0, C_mo(1, nocc+1), nao, temp1, nao, 0.0d0, sigma_mo, nvir)

        deallocate(temp1)
    end subroutine transform_sigma_ao_mo











end module davidon
