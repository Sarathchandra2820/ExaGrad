module fock_builder_module
    use iso_c_binding
    use molecule_loader, only: nbas, nao
    use molecule_t, only: molecule
    use libcint_interface
    use math_utils
    use jk_contraction_module
    use two_electron_cholesky_module
    use two_electron_df_module
    implicit none

    integer, parameter :: DEFAULT_BLOCK_SIZE = 16

contains

    function nuclear_repulsion(mol) result(E_nuc)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double) :: E_nuc
        integer :: A, B, pA, pB
        real(c_double) :: ZA, ZB, dR(3)

        E_nuc = 0.0d0
        do A = 1, int(mol%basis%natm)
            ZA = dble(mol%basis%atm(1, A))
            pA = int(mol%basis%atm(2, A)) + 1
            do B = 1, A-1
                ZB = dble(mol%basis%atm(1, B))
                pB = int(mol%basis%atm(2, B)) + 1
                dR = mol%basis%env(pA:pA+2) - mol%basis%env(pB:pB+2)
                E_nuc = E_nuc + ZA * ZB / sqrt(sum(dR**2))
            end do
        end do
    end function nuclear_repulsion

    subroutine normalize_fock_method(method_in, method_out)
        implicit none
        character(len=*), intent(in)  :: method_in
        character(len=*), intent(out) :: method_out
        character(len=64) :: tmp

        tmp = adjustl(trim(method_in))
        select case (trim(tmp))
        case ('cholesky', 'CHOLESKY', 'df', 'DF')
            method_out = 'cholesky'
        case ('block_cholesky', 'BLOCK_CHOLESKY', 'block_df', 'BLOCK_DF', &
              'df_block', 'DF_BLOCK', 'blocked_df', 'BLOCKED_DF')
            method_out = 'block_cholesky'
        case ('true_df', 'TRUE_DF', 'density_fitting', 'DENSITY_FITTING')
            method_out = 'true_df'
        case default
            method_out = 'direct'
        end select
    end subroutine normalize_fock_method

    pure function fock_method_banner(method) result(label)
        implicit none
        character(len=*), intent(in) :: method
        character(len=40) :: label
        select case (trim(method))
        case ('true_df')
            label = 'RHF SCF (True-DF-JK)'
        case ('block_cholesky')
            label = 'RHF SCF (Cholesky-Block-JK)'
        case ('cholesky')
            label = 'RHF SCF (Cholesky-JK)'
        case default
            label = 'RHF SCF (Direct-JK)'
        end select
    end function fock_method_banner

    subroutine initialize_fock_backend(method, Q)
        implicit none
        character(len=*), intent(in)  :: method
        real(c_double),   intent(out) :: Q(nbas,nbas)

        select case (trim(method))
        case ('true_df')
            Q = 0.0d0
            call initialize_true_df_factors()
        case ('block_cholesky', 'cholesky')
            Q = 0.0d0
            call initialize_cholesky_factors()
        case default
            call compute_schwarz(Q)
        end select
    end subroutine initialize_fock_backend

    subroutine build_fock(method, P, Hcore, Q, F)
        implicit none
        character(len=*), intent(in)  :: method
        real(c_double),   intent(in)  :: P(nao,nao), Hcore(nao,nao), Q(nbas,nbas)
        real(c_double),   intent(out) :: F(nao,nao)

        real(c_double), allocatable :: J(:,:), K(:,:)

        allocate(J(nao,nao), K(nao,nao))

        select case (trim(method))

        case ('true_df')
            if (.not. true_df_ready) call initialize_true_df_factors()
            if (true_df_naux <= 0)   stop 'True-DF initialisation failed: no auxiliary rank'
            call contract_jk(true_df_B, true_df_naux, P, J, K)

        case ('block_cholesky')
            if (.not. cholesky_ready) call initialize_cholesky_factors()
            if (cholesky_naux <= 0)   stop 'Cholesky initialisation failed: no auxiliary rank'
            call contract_jk(cholesky_B, cholesky_naux, P, J, K, block_size=DEFAULT_BLOCK_SIZE)

        case ('cholesky')
            if (.not. cholesky_ready) call initialize_cholesky_factors()
            if (cholesky_naux <= 0)   stop 'Cholesky initialisation failed: no auxiliary rank'
            call contract_jk(cholesky_B, cholesky_naux, P, J, K)

        case default
            call contract_jk_direct(P, Q, J, K)

        end select

        ! For the DF/Cholesky paths K is the raw exchange sum; apply the 1/2 prefactor here.
        ! For the direct path the 1/2 is already inside contract_jk_direct (accumulated per ERI).
        select case (trim(method))
        case ('true_df', 'block_cholesky', 'cholesky')
            F = Hcore + J - 0.5d0 * K
        case default
            F = Hcore + J - K
        end select

        deallocate(J, K)
    end subroutine build_fock

    subroutine clear_fock_backend(method)
        implicit none
        character(len=*), intent(in) :: method

        select case (trim(method))
        case ('true_df')
            call clear_true_df_factors()
        case ('block_cholesky', 'cholesky')
            call clear_cholesky_factors()
        case default
        end select
    end subroutine clear_fock_backend

end module fock_builder_module
