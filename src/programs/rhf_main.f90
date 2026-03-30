program rhf_main
    use iso_c_binding
    use one_eints
    use molecule_t, only: molecule
    use scf_module
    use polarisability_init
    implicit none

    type(molecule) :: mol
    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)
    real(c_double), allocatable :: C_mo(:,:), dip_ao(:,:,:), dip_mo(:,:,:)
    real(c_double) :: mu_elec(3)
    integer :: k, nocc, nvir

    print *, "=========================================="
    print *, " Fortran Direct RHF SCF powered by libcint"
    print *, "=========================================="

    call init_molecule(mol)
    call build_hcore_overlap(mol, S, Hcore)

    call build_s_inv(S, S_inv_sqrt)

    call build_initial_guess(mol, Hcore, S_inv_sqrt, P)

    nocc = mol%total_electrons / 2
    nvir = int(mol%basis%nao) - nocc

    allocate(C_mo(int(mol%basis%nao),int(mol%basis%nao)))
    allocate(dip_ao(int(mol%basis%nao),int(mol%basis%nao),3))
    allocate(dip_mo(nvir,nocc,3))

    call run_scf(mol, S, Hcore, S_inv_sqrt, P, C_mo)
    call transform_dipole_integrals(mol, C_mo, nocc, dip_ao, dip_mo)

    do k = 1, 3
        mu_elec(k) = -sum(P * dip_ao(:,:,k))
    end do
    print '(A,3F20.10)', '  MU_ELEC_AU =', mu_elec(1), mu_elec(2), mu_elec(3)

end program rhf_main
