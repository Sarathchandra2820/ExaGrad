program rhf_main
    use iso_c_binding
    use one_eints
    use scf_module
    use cpks
    implicit none

    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)
    real(c_double), allocatable :: C_mo(:,:), dip_ao(:,:,:), dip_mo(:,:,:)
    real(c_double) :: mu_elec(3)
    integer :: k

    print *, "=========================================="
    print *, " Fortran Direct RHF SCF powered by libcint"
    print *, "=========================================="

    call init_molecule()
    call build_hcore_overlap(S, Hcore)

    allocate(S_inv_sqrt(nao, nao))
    call build_s_inv(nao, S, S_inv_sqrt)

    allocate(P(nao, nao))
    call build_initial_guess(nao, total_electrons, Hcore, S_inv_sqrt, P)

    allocate(C_mo(nao,nao), dip_ao(nao,nao,3), dip_mo(3,nao,nao))

    call run_scf(S, Hcore, S_inv_sqrt, P, C_mo)
    call transform_dipole_integrals(C_mo, dip_ao, dip_mo)

    do k = 1, 3
        mu_elec(k) = -sum(P * dip_ao(:,:,k))
    end do
    print '(A,3F20.10)', '  MU_ELEC_AU =', mu_elec(1), mu_elec(2), mu_elec(3)

end program rhf_main
