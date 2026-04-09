program rhf_rose_main

    use iso_c_binding
    use molecule_t,          only: molecule, fragment_scf_info
    use molecule_loader,     only: init_molecule, activate_molecule
    use one_eints,           only: build_hcore_overlap, build_s_inv, build_initial_guess
    use scf_module,          only: run_scf
    use fock_builder_module,  only: normalize_fock_method
    use rose_interface_module, only: perform_frag_rhf, run_rose_localization

    implicit none

    ! --- supersystem ---
    type(molecule) :: mol
    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)
    real(c_double), allocatable :: C_mo(:,:), mo_energies(:)
    real(c_double), allocatable :: lmo_occ(:,:), lmo_vir(:,:)
    integer :: nocc, nvir, k, nfrags
    ! --- fragments ---
    type(molecule), allocatable :: mol_frag(:)
    type(fragment_scf_info), allocatable :: frag_scf_info(:)

    
    ! --- I/O and config ---
    character(len=1024) :: mol_dir
    character(len=64) :: method_env
    character(len=32) :: fock_method
    integer :: env_len, env_stat

    print *, "=========================================="
    print *, " ROSE-RHF: Supersystem + Fragment SCF"
    print *, "=========================================="

    ! ------------------------------------------------------------------
    ! Determine Fock builder method from environment
    ! ------------------------------------------------------------------
    fock_method = 'direct'
    call get_environment_variable('EXAGRAD_FOCK_BUILDER', method_env, &
                                  length=env_len, status=env_stat)
    if (env_stat == 0) then
        method_env = adjustl(method_env(1:env_len))
        call normalize_fock_method(trim(method_env), fock_method)
    end if

    ! ------------------------------------------------------------------
    ! Determine molecule directory from environment
    ! ------------------------------------------------------------------
    mol_dir = '.'
    call get_environment_variable('EXAGRAD_MOL_DIR', mol_dir, &
                                  length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) mol_dir = adjustl(mol_dir(1:env_len))
    

    call perform_frag_rhf(mol_dir, nfrags, mol_frag, frag_scf_info)

    nocc = frag_scf_info(0)%nocc
    nvir = frag_scf_info(0)%nvir
    C_mo = frag_scf_info(0)%C_mo
    mo_energies = frag_scf_info(0)%mo_energies
    
    print *, ''
    print *, 'Supersystem occupied MO energies (eV):'
    do k = 1, nocc
        print '(F20.10)', frag_scf_info(0)%mo_energies(k) * 27.211386245988d0
    end do
    
    print *, 'Supersystem virtual MO energies (eV):'
    do k = nocc+1, nocc+nvir
        print '(F20.10)', frag_scf_info(0)%mo_energies(k) * 27.211386245988d0
    end do



    call run_rose_localization(frag_scf_info, nfrags, lmo_occ, lmo_vir)

    print *, ''
    print *, '==> ROSE-RHF Done.'

end program rhf_rose_main
