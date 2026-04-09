program rhf_main
    use iso_c_binding
    use one_eints
    use molecule_t,        only: molecule
    use molecule_loader
    use scf_module
    use fock_builder_module, only: normalize_fock_method
    use polarisability_init
    use cpks
    implicit none

    type(molecule) :: mol
    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)
    real(c_double), allocatable :: C_mo(:,:), dip_ao(:,:,:), dip_mo(:,:,:), mo_energies(:)
    real(c_double) :: mu_elec(3), alpha(3,3)
    integer :: k, nocc, nvir
    character(len=64) :: method_env
    character(len=32) :: fock_method
    integer :: env_len, env_stat

    print *, "=========================================="
    print *, " Fortran Direct RHF SCF powered by libcint"
    print *, "=========================================="

    ! Determine integral backend (same env var as SCF driver)
    fock_method = 'direct'
    call get_environment_variable('EXAGRAD_FOCK_BUILDER', method_env,length=env_len, status=env_stat)
    if (env_stat == 0) then
        method_env = adjustl(method_env(1:env_len))
        call normalize_fock_method(trim(method_env), fock_method)
    end if

    call init_molecule(mol)
    call build_hcore_overlap(mol, S, Hcore)
    call build_s_inv(S, S_inv_sqrt)
    call build_initial_guess(mol, Hcore, S_inv_sqrt, P)

    nocc = mol%total_electrons / 2
    nvir = int(mol%basis%nao) - nocc

    allocate(C_mo(int(mol%basis%nao), int(mol%basis%nao)))
    allocate(mo_energies(nocc + nvir))
    allocate(dip_ao(int(mol%basis%nao), int(mol%basis%nao), 3))
    allocate(dip_mo(nvir, nocc, 3))

    call run_scf(mol, S, Hcore, S_inv_sqrt, P, C_mo, mo_energies)



    ! Dipole integrals in AO and MO basis
    call transform_dipole_integrals(mol, C_mo, nocc, dip_ao, dip_mo)

    print *, ''
    print *, 'The occupied MO energies (eV) are:'
    do k = 1, nocc
        print '(F20.10)', mo_energies(k) * 27.211386245988d0
    end do
    print *, 'The virtual MO energies (eV) are:'
    do k = nocc+1, nocc+nvir
        print '(F20.10)', mo_energies(k) * 27.211386245988d0
    end do

    do k = 1, 3
        mu_elec(k) = -sum(P * dip_ao(:,:,k))
    end do
    print '(A,3F20.10)', '  MU_ELEC_AU =', mu_elec(1), mu_elec(2), mu_elec(3)

    ! Solve CPHF and compute static polarizability
    call solve_cpks(C_mo, mo_energies, nocc, dip_mo, fock_method, alpha)

    print *, ''
    print *, '  ===== Static Dipole Polarizability Tensor (a.u.) ====='
    print '(5X,A12,A12,A12)', 'x', 'y', 'z'
    do k = 1, 3
        print '(A5,3F12.6)', char(119+k)//' ', alpha(k,1), alpha(k,2), alpha(k,3)
    end do
    print '(A,F12.6)', '  Isotropic alpha (a.u.) =', &
          (alpha(1,1) + alpha(2,2) + alpha(3,3)) / 3.0d0

end program rhf_main
