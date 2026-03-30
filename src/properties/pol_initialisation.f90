module polarisability_init
    use math_utils_module
    use iso_c_binding
    use davidson_module
    use jk_contraction_module
    use two_electron_df_module, only: true_df_B, true_df_naux
    implicit none



    subroutine contract_over_kernel(sigma_ao, B, naux, method, Q)
        character(len=*), intent(in)  :: method
        integer, intent(in) :: naux
        integer :: nao
        real(c_double), intent(in) :: B(:,:,:)
        real(c_double), intent(inout) :: sigma_ao(:,:)
        real(c_double), intent(in), optional :: Q(:,:)
        real(c_double), allocatable :: J(:,:), K(:,:)

        nao = size(sigma_ao, 1)

        allocate(J(nao,nao), K(nao,nao))

        select case (trim(method))
        case ('true_df', 'block_cholesky', 'cholesky')
            if (naux <= 0) stop 'Error: no auxiliary rank for DF/Cholesky contraction'
            call contract_jk(B, naux, sigma_ao, J, K)
            sigma_ao = J - 0.5d0 * K
        case default
            if (.not. present(Q)) stop 'Error: Schwarz bounds Q required for direct contraction'
            call contract_jk_direct(sigma_ao, Q, J, K)
            sigma_ao = J - K
        end select

        deallocate(J, K)

    end subroutine contract_over_kernel

    subroutine build_sigma(method, sigma_init, C_mo, nocc, B, naux, sigma_out, Q)
        character(len=*), intent(in)  :: method
        integer, intent(in) :: nocc, naux
        real(c_double), intent(in) :: sigma_init(:,:)
        real(c_double), intent(in) :: C_mo(:,:)
        real(c_double), intent(in) :: B(:,:,:)
        real(c_double), allocatable, intent(out) :: sigma_out(:,:)
        real(c_double), intent(in), optional :: Q(:,:)

        real(c_double), allocatable :: sigma_ao(:,:)

        !Here we assume sigma_init is in MO basis, we need to transform it to AO basis, contract over the kernel and then transform back to MO basis

        call transform_sigma_mo_ao(sigma_init, C_mo, sigma_ao)
        call contract_over_kernel(sigma_ao, B, naux, method, Q)
        call transform_sigma_ao_mo(sigma_ao, C_mo, nocc, sigma_out)

        deallocate(sigma_ao)

    end subroutine build_sigma