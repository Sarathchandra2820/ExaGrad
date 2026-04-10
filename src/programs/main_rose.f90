program rhf_rose_main

    use iso_c_binding
    use molecule_t,          only: molecule, fragment_scf_info
    use molecule_loader,     only: init_molecule, activate_molecule
    use one_eints,           only: build_hcore_overlap, build_s_inv, build_initial_guess
    use scf_module,          only: run_scf
    use fock_builder_module,  only: normalize_fock_method
    use rose_interface_module, only: perform_frag_rhf, run_rose_localization, normalize_rose_orbital_ordering

    implicit none

    ! --- supersystem ---
    type(molecule) :: mol
    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)
    real(c_double), allocatable :: C_mo(:,:), mo_energies(:)
    real(c_double), allocatable :: lmo_occ(:,:), lmo_vir(:,:), lmo_vir_hard(:,:)
    real(c_double) :: t0, t1, t_scf, t_loc, t_total
    integer :: nocc, nvir, k, nfrags
    ! --- fragments ---
    type(molecule), allocatable :: mol_frag(:)
    type(fragment_scf_info), allocatable :: frag_scf_info(:)

    
    ! --- I/O and config ---
    character(len=1024) :: mol_dir
    character(len=64) :: method_env, rose_order_env
    character(len=32) :: fock_method, rose_ordering
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

    rose_ordering = 'fragment'
    call get_environment_variable('EXAGRAD_ROSE_ORBITAL_ORDER', rose_order_env, &
                                  length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
        rose_order_env = adjustl(rose_order_env(1:env_len))
        call normalize_rose_orbital_ordering(trim(rose_order_env), rose_ordering)
    end if

    ! ------------------------------------------------------------------
    ! Determine molecule directory from environment
    ! ------------------------------------------------------------------
    mol_dir = '.'
    call get_environment_variable('EXAGRAD_MOL_DIR', mol_dir, &
                                  length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) mol_dir = adjustl(mol_dir(1:env_len))
    

    call cpu_time(t0)
    call perform_frag_rhf(mol_dir, nfrags, mol_frag, frag_scf_info)
    call cpu_time(t1)
    t_scf = t1 - t0

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

    call cpu_time(t0)
    call run_rose_localization(frag_scf_info, nfrags, lmo_occ, lmo_vir, rose_ordering, lmo_vir_hard)
    call cpu_time(t1)
    t_loc = t1 - t0
    t_total = t_scf + t_loc

    print *, ''
    print *, 'Timing summary (CPU seconds):'
    print '(A,F10.4)', '  SCF/setup phase        : ', t_scf
    print '(A,F10.4)', '  Localization phase     : ', t_loc
    print '(A,F10.4)', '  SCF + localization     : ', t_total

    print *, ''
    print *, '==> ROSE-RHF Done.'

end program rhf_rose_main
