module sigma_build 

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

        integer :: nao, nocc, nvir
        real(c_double), intent(in) :: sigma_mo(:,:), C_mo(:,:)
        real(c_double), allocatable, intent(out) :: sigma_ao(:,:)
        real(c_double), allocatable :: temp1(:,:), temp2(:,:)


        integer :: i, j


        ! Now we need to transform U_ai from MO to AO basis
        ! U_ai is nvir x nocc, we want to get sigma_ij in AO basis
        ! sigma_ij = sum_{a,i} C_ia * U_ai
        nao = size(C_mo, 1)
        nocc = size(sigma_mo, 2)
        nvir = size(sigma_mo, 1)

       
        ! First transform U_ai to AO basis

        ! temp = C_virt * U_ai * C_occ^T
        ! sigma = temp + temp^T (since sigma is symmetric)

        !calculate temp = C_vir * U_ai
        call dgemm('N','N', nao,nocc,nvir, 1.0d0, C_mo(:, nocc+1:), nao, sigma_mo, nvir, 0.0d0, temp1, nao)
        call dgemm('N','T', nao,nao,nocc, 1.0d0, temp1, nao, C_mo(:,1:nocc), nao, 0.0d0, temp2, nao)

        !Writing the code to be memory efficient, we are allocating 3*NAO^2 for temp1 and temp2, but we can do the transformation in place to save memory.
        do i = 1, nao
            do j = 1, nao
                sigma_ao(i,j) = temp2(i,j) + temp2(j,i)
            end do
        end do

      


    end subroutine transform_sigma_mo_ao

    subroutine transform_sigma_ao_mo(sigma_ao, C_mo, nocc, sigma_mo)
        integer :: nao, nvir
        integer, intent(in) :: nocc
        real(c_double), intent(in) :: sigma_ao(:,:), C_mo(:,:)
        real(c_double), allocatable, intent(out) :: sigma_mo(:,:)

        real(c_double), allocatable :: temp1(:,:)

        nao = size(C_mo, 1)
        nvir = size(C_mo, 2) - nocc


        ! Transform sigma from AO to MO basis
        ! sigma_ai = C_vir^T * sigma_ao * C_occ
        ! Step 1: temp(nao,nocc) = sigma_ao * C_occ
        call dgemm('N','N', nao, nocc, nao, 1.0d0, sigma_ao, nao, C_mo(:,1:nocc), nao, 0.0d0, temp1, nao)
        ! Step 2: sigma_mo(nvir,nocc) = C_vir^T * temp
        call dgemm('T','N', nvir, nocc, nao, 1.0d0, C_mo(:, nocc+1:), nao, temp1, nao, 0.0d0, sigma_mo, nvir)


    end subroutine transform_sigma_ao_mo

    subroutine contract_sigma_over_kernel(method, sigma_init, C_mo, nocc, B, naux, sigma_out, Q)
        character(len=*), intent(in)  :: method
        integer, intent(in) :: nocc, naux
        integer :: nao
        real(c_double), intent(in) :: sigma_init(:,:)
        real(c_double), intent(in) :: C_mo(:,:)
        real(c_double), intent(in) :: B(:,:,:)
        real(c_double), allocatable, intent(out) :: sigma_out(:,:)
        real(c_double), intent(in), optional :: Q(:,:)

        real(c_double), allocatable :: sigma_ao(:,:), J(:,:), K(:,:)

        ! MO -> AO
        call transform_sigma_mo_ao(sigma_init, C_mo, sigma_ao)

        ! Contract over kernel in AO basis
        nao = size(sigma_ao, 1)

        select case (trim(method))
        case ('true_df', 'block_cholesky', 'cholesky')
            if (naux <= 0) stop 'Error: no auxiliary rank for DF/Cholesky contraction'
            call contract_jk(B, naux, sigma_ao, J, K)
            call daxpy(nao*nao, -0.5d0, K, 1, J, 1)
            sigma_ao = J
        case default
            if (.not. present(Q)) stop 'Error: Schwarz bounds Q required for direct contraction'
            call contract_jk_direct(sigma_ao, Q, J, K)
            call daxpy(nao*nao, -1.0d0, K, 1, J, 1)
            sigma_ao = J
        end select


        ! AO -> MO
        call transform_sigma_ao_mo(sigma_ao, C_mo, nocc, sigma_out)


    end subroutine contract_sigma_over_kernel

    subroutine build_energy_denominator(mo_energies, nocc, nvir, energy_denom)
        real(c_double), intent(in) :: mo_energies(:)
        integer, intent(in) :: nocc, nvir
        real(c_double), allocatable, intent(out) :: energy_denom(:,:)

        integer :: i, j


        do j = 1, nvir
            do i = 1, nocc
                energy_denom(i,j) = 1.0d0/ (mo_energies(nocc+j) - mo_energies(i))
            end do
        end do

    end subroutine build_energy_denominator

    subroutine contract_wt_denominator(sigma_in, energy_denom, sigma_out)
        implicit none 
        real(c_double), intent(in) :: sigma_in(:,:), energy_denom(:,:)
        real(c_double), allocatable, intent(inout) :: sigma_out(:,:)

        integer :: i, j
        integer :: nocc, nvir


        nocc = size(sigma_in, 2)
        nvir = size(sigma_in, 1)


        do j = 1, nocc
            do i = 1, nvir
                sigma_out(i,j) = sigma_in(i,j) * energy_denom(i,j)
            end do
        end do

    end subroutine contract_wt_denominator

    

end module sigma_build

