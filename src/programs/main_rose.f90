program rhf_rose_main

    use iso_c_binding
    use molecule_t,          only: molecule
    use molecule_loader,     only: init_molecule, activate_molecule
    use one_eints,           only: build_hcore_overlap, build_s_inv, build_initial_guess
    use scf_module,          only: run_scf
    use fock_builder_module,  only: normalize_fock_method
    use rose_interface_module, only: localize_mo_spaces

    implicit none

    ! --- supersystem ---
    type(molecule) :: mol
    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)
    real(c_double), allocatable :: C_mo(:,:), mo_energies(:)
    integer :: nocc, nvir

    ! --- fragments ---
    integer :: nfrags, k
    type(molecule), allocatable :: mol_frag(:)
    real(c_double), allocatable :: S_frag(:,:), Hcore_frag(:,:)
    real(c_double), allocatable :: S_inv_sqrt_frag(:,:), P_frag(:,:)
    real(c_double), allocatable :: C_mo_frag(:,:,:), mo_energies_frag(:,:)
    integer :: nocc_frag, nvir_frag, nao_frag

    ! --- I/O and config ---
    character(len=1024) :: mol_dir, frag_dir_name
    character(len=2048) :: frag_path
    character(len=64) :: method_env
    character(len=32) :: fock_method
    integer :: env_len, env_stat, unit_id, ios

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
    if (env_stat == 0) mol_dir = adjustl(mol_dir(1:env_len))

    ! ------------------------------------------------------------------
    ! Read rose_info.txt manifest to discover fragments
    ! ------------------------------------------------------------------
    open(newunit=unit_id, file=trim(mol_dir)//'/rose_info.txt', &
         status='old', iostat=ios)
    if (ios /= 0) then
        print *, "ERROR: Cannot open ", trim(mol_dir)//'/rose_info.txt'
        print *, "  Run the Python integral exporter first."
        stop 1
    end if
    read(unit_id, *) nfrags
    allocate(mol_frag(nfrags))

    print '(A,I0,A)', ' Found ', nfrags, ' fragment(s)'
    print *, ''

    ! ------------------------------------------------------------------
    ! Fragment SCF loop (sequential — each activates its own globals)
    ! ------------------------------------------------------------------
    do k = 1, nfrags
        read(unit_id, '(A)') frag_dir_name
        frag_path = trim(mol_dir)//'/'//trim(frag_dir_name)

        print '(A,I0,A,A)', ' --- Fragment ', k, ': ', trim(frag_dir_name)
        print *, '------------------------------------------------------------'

        call init_molecule(mol_frag(k), trim(frag_path))
        call activate_molecule(mol_frag(k), trim(frag_path))

        call build_hcore_overlap(mol_frag(k), S_frag, Hcore_frag)
        call build_s_inv(S_frag, S_inv_sqrt_frag)
        call build_initial_guess(mol_frag(k), Hcore_frag, S_inv_sqrt_frag, P_frag)

        nao_frag  = int(mol_frag(k)%basis%nao)
        nocc_frag = mol_frag(k)%total_electrons / 2
        nvir_frag = nao_frag - nocc_frag

        allocate(C_mo_frag(nao_frag, nao_frag,nfrags))
        allocate(mo_energies_frag(nao_frag,nfrags))

        call run_scf(mol_frag(k), S_frag, Hcore_frag, S_inv_sqrt_frag, &
                     P_frag, C_mo_frag(:,:,k), mo_energies_frag(:,k))

        print '(A,I0,A,F20.10)', ' Fragment ', k, ' converged.'
        print *, ''

        deallocate(C_mo_frag, mo_energies_frag)
        deallocate(S_frag, Hcore_frag, S_inv_sqrt_frag, P_frag)
    end do
    close(unit_id)

    ! ------------------------------------------------------------------
    ! Supersystem SCF
    ! ------------------------------------------------------------------
    print *, '============================================================'
    print *, ' --- Supersystem SCF ---'
    print *, '============================================================'

    call init_molecule(mol, trim(mol_dir))
    call activate_molecule(mol, trim(mol_dir))

    call build_hcore_overlap(mol, S, Hcore)
    call build_s_inv(S, S_inv_sqrt)
    call build_initial_guess(mol, Hcore, S_inv_sqrt, P)

    nocc = mol%total_electrons / 2
    nvir = int(mol%basis%nao) - nocc

    allocate(C_mo(int(mol%basis%nao), int(mol%basis%nao)))
    allocate(mo_energies(nocc + nvir))

    call run_scf(mol, S, Hcore, S_inv_sqrt, P, C_mo, mo_energies)

    print *, ''
    print *, 'Supersystem occupied MO energies (eV):'
    do k = 1, nocc
        print '(F20.10)', mo_energies(k) * 27.211386245988d0
    end do
    

    print *, ''
    print *, '==> ROSE-RHF Done.'

end program rhf_rose_main
