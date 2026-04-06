module cpks

    use iso_c_binding,                only: c_double
    use davidson_module,              only: davidson_linear_solve
    use sigma_build,                  only: contract_sigma_over_kernel
    use two_electron_df_module,       only: true_df_B, true_df_naux, initialize_true_df_factors
    use two_electron_cholesky_module, only: cholesky_B, cholesky_naux
    use jk_contraction_module,        only: compute_schwarz
    use one_eints,                    only: nbas
    use math_utils
    implicit none


    ! -----------------------------------------------------------------------
    ! Module-level context variables — set by solve_cpks before Davidson runs.
    ! These act as a "thunk" so that the fixed-interface callbacks (cphf_sigma,
    ! cphf_precond) can access all necessary data without extra arguments.
    ! -----------------------------------------------------------------------
    integer :: cpks_nocc, cpks_nvir, cpks_nao, cpks_naux
    real(c_double), allocatable :: cpks_C_mo(:,:)       ! (nao, nmo)
    real(c_double), allocatable :: cpks_eps_diff(:,:)   ! (nvir, nocc): e_a - e_i
    real(c_double), allocatable :: cpks_B(:,:,:)        ! DF/Cholesky tensor
    real(c_double), allocatable :: cpks_Q(:,:)          ! Schwarz bounds (direct only)
    character(len=32) :: cpks_method
    integer :: cpks_method_id = -1    ! 0=direct, 1=cholesky, 2=true_df
    real(c_double), allocatable, save :: U_scratch(:,:), G_scratch(:,:)
    real(c_double), allocatable, save :: R_scratch(:,:), Z_scratch(:,:)
!$omp threadprivate(U_scratch, G_scratch, R_scratch, Z_scratch)


contains


    ! -----------------------------------------------------------------------
    ! CPHF matrix-vector product:  Ax = [diag(e_a - e_i) + 2*G] * x
    !
    ! G is the two-electron response from contract_sigma_over_kernel, which
    ! returns (J[D_sym] - 0.5*K[D_sym]) in MO basis.  Multiplying by 2 gives
    ! (2*J - K) = the full 4-index CPHF response (see derivation in plan).
    ! -----------------------------------------------------------------------
    subroutine cphf_sigma(x, Ax, n)
        integer,        intent(in)  :: n
        real(c_double), intent(in)  :: x(n)
        real(c_double), intent(out) :: Ax(n)

        if (.not. allocated(U_scratch)) then
            allocate(U_scratch(cpks_nvir, cpks_nocc), G_scratch(cpks_nvir, cpks_nocc))
        else if (size(U_scratch,1) /= cpks_nvir .or. size(U_scratch,2) /= cpks_nocc) then
            deallocate(U_scratch, G_scratch)
            allocate(U_scratch(cpks_nvir, cpks_nocc), G_scratch(cpks_nvir, cpks_nocc))
        end if

        call dcopy(n, x, 1, U_scratch, 1)

        ! Two-electron response; factor of 2 converts (J-0.5K) -> (2J-K)
        if (cpks_method_id == 0) then
            call contract_sigma_over_kernel(cpks_method, U_scratch, cpks_C_mo, cpks_nocc, &
                                             cpks_B, cpks_naux, G_scratch, cpks_Q)
        else
            call contract_sigma_over_kernel(cpks_method, U_scratch, cpks_C_mo, cpks_nocc, &
                                             cpks_B, cpks_naux, G_scratch)
        end if

        ! Full orbital Hessian action: (e_a - e_i)*U + 2*G
        G_scratch = cpks_eps_diff * U_scratch + 2.0d0 * G_scratch
        call dcopy(n, G_scratch, 1, Ax, 1)

    end subroutine cphf_sigma


    ! -----------------------------------------------------------------------
    ! Diagonal preconditioner:  z_ai = r_ai / (e_a - e_i)
    ! -----------------------------------------------------------------------
    subroutine cphf_precond(r, z, n)
        integer,        intent(in)  :: n
        real(c_double), intent(in)  :: r(n)
        real(c_double), intent(out) :: z(n)

        if (.not. allocated(R_scratch)) then
            allocate(R_scratch(cpks_nvir, cpks_nocc), Z_scratch(cpks_nvir, cpks_nocc))
        else if (size(R_scratch,1) /= cpks_nvir .or. size(R_scratch,2) /= cpks_nocc) then
            deallocate(R_scratch, Z_scratch)
            allocate(R_scratch(cpks_nvir, cpks_nocc), Z_scratch(cpks_nvir, cpks_nocc))
        end if

        call dcopy(n, r, 1, R_scratch, 1)
        Z_scratch = R_scratch / cpks_eps_diff
        call dcopy(n, Z_scratch, 1, z, 1)

    end subroutine cphf_precond


    ! -----------------------------------------------------------------------
    ! Solve the CPHF linear equations for all 3 field directions and compute
    ! the static dipole polarizability tensor.
    !
    ! Equations solved: A * U_k = -d_k    for k = x, y, z
    ! where A_{ai,bj} = delta(e_a - e_i) + [4*(ai|jb) - (aj|ib) - (ab|ij)]
    ! and d_k_{ai} = <a|r_k|i> are the dipole MO integrals.
    !
    ! Polarizability: alpha(k,l) = -4 * Tr(d_k * U_l)
    ! Closed-shell first-order density carries spin factor 2, and D^(1)
    ! also contains symmetric occ-vir + vir-occ contributions.
    ! -----------------------------------------------------------------------
    subroutine solve_cpks(C_mo, mo_energies, nocc, dip_mo, method, alpha)

        real(c_double), intent(in)    :: C_mo(:,:)
        real(c_double), intent(in)    :: mo_energies(:)
        integer,        intent(in)    :: nocc
        real(c_double), intent(in)    :: dip_mo(:,:,:)   ! (nvir, nocc, 3)
        character(len=*), intent(in)  :: method
        real(c_double), intent(out)   :: alpha(3,3)

        real(c_double), allocatable :: b(:), U_flat(:), U_x(:,:,:)
        integer :: n, k, l, a, i
        logical :: converged
        logical :: converged_dir(3)
        integer :: max_sub

        ! --- Populate module-level context ----------------------------------
        cpks_nao  = size(C_mo, 1)
        cpks_nocc = nocc
        cpks_nvir = size(C_mo, 2) - nocc
        n         = cpks_nvir * cpks_nocc
        cpks_method = trim(method)

        allocate(cpks_C_mo(cpks_nao, size(C_mo, 2)))
        cpks_C_mo = C_mo

        ! Build (e_a - e_i) matrix, shape (nvir, nocc)
        allocate(cpks_eps_diff(cpks_nvir, cpks_nocc))
        do i = 1, cpks_nocc
            do a = 1, cpks_nvir
                cpks_eps_diff(a, i) = mo_energies(nocc + a) - mo_energies(i)
            end do
        end do

        ! Backend tensor / Schwarz bounds
        select case (trim(method))
        case ('true_df')
            cpks_method_id = 2
            call initialize_true_df_factors()
            cpks_naux = true_df_naux
            allocate(cpks_B(cpks_nao, cpks_nao, cpks_naux))
            cpks_B = true_df_B

        case ('block_cholesky', 'cholesky')
            cpks_method_id = 1
            cpks_naux = cholesky_naux
            allocate(cpks_B(cpks_nao, cpks_nao, cpks_naux))
            cpks_B = cholesky_B

        case default   ! direct
            cpks_method_id = 0
            cpks_naux = 0
            allocate(cpks_B(1, 1, 1))
            cpks_B = 0.0d0
            allocate(cpks_Q(nbas, nbas))
            call compute_schwarz(cpks_Q)

        end select

        ! --- Solve CPHF for each Cartesian direction ------------------------
        allocate(U_x(cpks_nvir, cpks_nocc, 3))
        converged_dir = .false.

        ! Subspace size: cap at problem dimension, use at most 30
        max_sub = min(n, 30)

        print *, ''
        print *, '  ===== CPHF/CPKS Polarizability Solver ====='
        print '(A,I0,A,I0,A,I0)', '  Dimensions: nocc=', cpks_nocc, &
              '  nvir=', cpks_nvir, '  n=', n

!$omp parallel do default(shared) private(k, b, U_flat, converged) schedule(static,1)
        do k = 1, 3
            allocate(b(n), U_flat(n))
            b      = -reshape(dip_mo(:,:,k), [n])
            U_flat = 0.0d0
            call davidson_linear_solve(cphf_sigma, cphf_precond, b, U_flat, n, &
                                        1.0d-8, 150, max_sub, converged)
            converged_dir(k) = converged
            U_x(:,:,k) = reshape(U_flat, [cpks_nvir, cpks_nocc])
            deallocate(b, U_flat)
        end do
!$omp end parallel do

        do k = 1, 3
            if (.not. converged_dir(k)) then
                print '(A,I1)', '  WARNING: CPHF did not converge for direction ', k
            end if
        end do

        ! --- Polarizability tensor ------------------------------------------
        ! alpha(k,l) = -4 * sum_{a,i} d_k(a,i) * U_l(a,i)
        do k = 1, 3
            do l = 1, 3
                alpha(k,l) = -4.0d0 * sum(dip_mo(:,:,k) * U_x(:,:,l))
            end do
        end do

        ! --- Cleanup --------------------------------------------------------
        deallocate(cpks_C_mo, cpks_eps_diff, cpks_B, U_x)
        if (allocated(cpks_Q)) deallocate(cpks_Q)

    end subroutine solve_cpks


end module cpks
