program rhf_rose_main

    use iso_c_binding
    use utils, only: write_matrix_to_file
    use molecule_t,          only: molecule, fragment_scf_info
    use molecule_loader,     only: init_molecule, activate_molecule
    use one_eints,           only: build_hcore_overlap, build_s_inv, build_initial_guess
    use scf_module,          only: run_scf
    use fock_builder_module,  only: normalize_fock_method
    use rose_interface_module, only: perform_frag_rhf, run_rose_localization, normalize_rose_orbital_ordering
    use polarisability_init, only: compute_dipole_ao, transform_dipole_ao_to_lmo, &
                                   transform_dipole_integrals
    use cpks,                only: solve_cpks
    use rose_davidson,       only: partition_rose_clmo


    implicit none

    ! --- supersystem ---
    type(molecule) :: mol
    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)
    real(c_double), allocatable :: C_mo(:,:), mo_energies(:)
    integer(c_int), allocatable :: nlo_frag_occ(:), nlo_frag_vir(:)
    real(c_double), allocatable :: c_lmo_occ(:,:), c_lmo_vir(:,:), lmo_vir_hard(:,:), U_occ(:,:), U_vir(:,:)
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

    mol = mol_frag(0)  ! Assuming the first fragment is the supersystem; adjust if needed
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
    call run_rose_localization(frag_scf_info, nfrags, nlo_frag_occ, nlo_frag_vir, c_lmo_occ, U_occ, c_lmo_vir, U_vir, rose_ordering, lmo_vir_hard)
    call cpu_time(t1)

    print*, 'The number of LMO per fragment (occupied and virtual) are:'
    do k = 1, nfrags
        print '(I4, 2I8)', k, nlo_frag_occ(k), nlo_frag_vir(k)
    end do
    t_loc = t1 - t0
    t_total = t_scf + t_loc

    print *, ''
    print *, 'Timing summary (CPU seconds):'
    print '(A,F10.4)', '  SCF/setup phase        : ', t_scf
    print '(A,F10.4)', '  Localization phase     : ', t_loc
    print '(A,F10.4)', '  SCF + localization     : ', t_total

    print *, ''
    print *, '==> ROSE-RHF Done.'

    block
        integer :: unit_lmo, ios
        character(len=1024) :: lmo_path

        lmo_path = trim(mol_dir)//'/lmo_coefficients.txt'
        call write_matrix_to_file(lmo_path, c_lmo_occ, status_msg='replace')
        call write_matrix_to_file(lmo_path, c_lmo_vir, status_msg='old')
        print '(A,A)', ' LMO coefficients written to ', trim(lmo_path)
    end block


    block
        real(c_double), allocatable :: dip_ao(:,:,:), dip_oo(:,:,:), dip_vv(:,:,:), dip_ov(:,:,:)
        real(c_double), allocatable :: dip_mo(:,:,:), dip_mo_k(:,:,:)
        real(c_double) :: alpha_total(3,3), alpha_k(3,3), alpha_sum(3,3)
        character(len=16) :: cx
        character(len=1024) :: component_path
        integer :: fk, nocc_k, nvir_k, xi

        ! Re-activate supersystem with its directory so that active_mol_dir
        ! is correct for DF integral reads (fragment SCF leaves it pointing
        ! at the last fragment's directory).
        call activate_molecule(mol, trim(mol_dir))

        allocate(dip_ao(size(C_mo, 1), size(C_mo, 1), 3))
        call compute_dipole_ao(mol, dip_ao)

        ! --- LMO dipole integrals (existing) ---
        call transform_dipole_ao_to_lmo(c_lmo_occ, c_lmo_vir, nocc, nvir, dip_ao, &
                                        dip_oo, dip_vv, dip_ov)
        do k = 1, 3
            write(cx, '(I0)') k
            component_path = trim(mol_dir)//'/dip_component_'//trim(cx)//'.txt'
            call write_matrix_to_file(component_path, dip_oo(:,:,k), status_msg='old')
            call write_matrix_to_file(component_path, dip_vv(:,:,k), status_msg='old')
            call write_matrix_to_file(component_path, dip_ov(:,:,k), status_msg='old')
        end do

        ! --- Total CPHF polarizability in canonical MO basis ---
        print *, ''
        print *, '=========================================='
        print *, ' CPHF Polarizability: Canonical MO basis'
        print *, '=========================================='
        call transform_dipole_integrals(mol, C_mo, nocc, dip_ao, dip_mo)
        call solve_cpks(C_mo, mo_energies, nocc, dip_mo, fock_method, alpha_total)
        deallocate(dip_mo)

        print *, ''
        print *, 'Total polarizability tensor (a.u.):'
        do xi = 1, 3
            print '(3F14.6)', alpha_total(xi,1), alpha_total(xi,2), alpha_total(xi,3)
        end do
        print '(A,F14.6)', '  Isotropic = ', &
            (alpha_total(1,1) + alpha_total(2,2) + alpha_total(3,3)) / 3.0d0

        ! --- Partition c_lmo into per-fragment C_mo blocks ---
        call partition_rose_clmo(frag_scf_info, nfrags, nlo_frag_occ, nlo_frag_vir, &
                                 c_lmo_occ, c_lmo_vir)

        ! --- Per-fragment CPHF polarizability ---
        print *, ''
        print *, '=========================================='
        print *, ' CPHF Polarizability: Per-fragment LMO basis'
        print *, '=========================================='
        alpha_sum = 0.0d0
        do fk = 1, nfrags
            nocc_k = frag_scf_info(fk)%nocc
            nvir_k = frag_scf_info(fk)%nvir
            call transform_dipole_integrals(mol, frag_scf_info(fk)%C_mo, nocc_k, dip_ao, dip_mo_k)
            call solve_cpks(frag_scf_info(fk)%C_mo, frag_scf_info(fk)%mo_energies, &
                            nocc_k, dip_mo_k, fock_method, alpha_k)
            deallocate(dip_mo_k)
            alpha_sum = alpha_sum + alpha_k
            print '(A,I0,A,F14.6)', '  Fragment ', fk, ' isotropic alpha = ', &
                (alpha_k(1,1) + alpha_k(2,2) + alpha_k(3,3)) / 3.0d0
        end do

        ! --- Additivity test ---
        print *, ''
        print *, '=========================================='
        print *, ' Additivity Test: sum(frags) vs total'
        print *, '=========================================='
        print *, 'Sum of fragment polarizability tensors (a.u.):'
        do xi = 1, 3
            print '(3F14.6)', alpha_sum(xi,1), alpha_sum(xi,2), alpha_sum(xi,3)
        end do
        print '(A,F14.6)', '  Isotropic (fragment sum) = ', &
            (alpha_sum(1,1) + alpha_sum(2,2) + alpha_sum(3,3)) / 3.0d0
        print '(A,F14.6)', '  Isotropic (total CPHF)   = ', &
            (alpha_total(1,1) + alpha_total(2,2) + alpha_total(3,3)) / 3.0d0
        print '(A,F14.6)', '  Abs difference           = ', &
            abs((alpha_sum(1,1)+alpha_sum(2,2)+alpha_sum(3,3)) &
               -(alpha_total(1,1)+alpha_total(2,2)+alpha_total(3,3))) / 3.0d0

    end block

end program rhf_rose_main

