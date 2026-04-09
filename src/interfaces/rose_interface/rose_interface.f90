module rose_interface_module

    use iso_c_binding
    use molecule_t,          only: molecule, fragment_scf_info
    use molecule_loader,     only: init_molecule, activate_molecule
    use one_eints,           only: build_hcore_overlap, build_s_inv, build_initial_guess
    use scf_module,          only: run_scf
    use fock_builder_module,  only: normalize_fock_method
    use make_ibo,            only: Localization, sort_on_fragments

    implicit none

contains

    subroutine perform_frag_rhf(mol_dir, nfrags, mol_frag, frag_scf_info)

        implicit none

        character(len=*), intent(in) :: mol_dir
        integer, intent(out) :: nfrags
        type(molecule), allocatable, intent(out) :: mol_frag(:)
        type(fragment_scf_info), allocatable, intent(out) :: frag_scf_info(:)
        real(c_double), allocatable :: C_mo_frag(:,:), mo_energies_frag(:)
        real(c_double), allocatable :: S_inv_sqrt_frag(:,:), P_frag(:,:), Hcore_frag(:,:),S_frag(:,:)
        integer :: nao_frag, nocc_frag, nvir_frag, nao_sm, nocc_sm, nvir_sm, k
        character(len=64) :: method_env
        character(len=32) :: fock_method
        integer :: env_len, env_stat, unit_id, ios
        character(len=1024) :: frag_dir_name
        character(len=2048) :: frag_path
        character(len=2048) :: resolved_mol_dir



        ! ------------------------------------------------------------------
        ! Determine molecule directory from environment
        ! ------------------------------------------------------------------
        resolved_mol_dir = trim(mol_dir)
        call get_environment_variable('EXAGRAD_MOL_DIR', resolved_mol_dir, &
                                    length=env_len, status=env_stat)
        if (env_stat == 0) resolved_mol_dir = adjustl(resolved_mol_dir(1:env_len))




        ! ------------------------------------------------------------------
        ! Read rose_info.txt manifest to discover fragments
        ! ------------------------------------------------------------------
        open(newunit=unit_id, file=trim(resolved_mol_dir)//'/rose_info.txt', &
            status='old', iostat=ios)
        if (ios /= 0) then
            print *, "ERROR: Cannot open ", trim(resolved_mol_dir)//'/rose_info.txt'
            print *, "  Run the Python integral exporter first."
            stop 1
        end if

        read(unit_id, *) nfrags
        allocate(mol_frag(0:nfrags))
        allocate(frag_scf_info(0:nfrags))


        call init_molecule(mol_frag(0), trim(resolved_mol_dir))
        call activate_molecule(mol_frag(0), trim(resolved_mol_dir))
        frag_scf_info(0)%mol = mol_frag(0)

        call build_hcore_overlap(mol_frag(0), S_frag, Hcore_frag)
        frag_scf_info(0)%S = S_frag
        call build_s_inv(S_frag, S_inv_sqrt_frag)
        call build_initial_guess(mol_frag(0), Hcore_frag, S_inv_sqrt_frag, P_frag)
        nao_sm  = int(mol_frag(0)%basis%nao)
        nocc_sm = mol_frag(0)%total_electrons / 2
        nvir_sm = nao_sm - nocc_sm
        allocate(C_mo_frag(nao_sm, nao_sm))
        allocate(mo_energies_frag(nao_sm))
        call run_scf(mol_frag(0), S_frag, Hcore_frag, S_inv_sqrt_frag, &
                    P_frag, C_mo_frag, mo_energies_frag)
        frag_scf_info(0)%C_mo = C_mo_frag(:,:)
        frag_scf_info(0)%mo_energies = mo_energies_frag(:)
        frag_scf_info(0)%nocc = nocc_sm
        frag_scf_info(0)%nvir = nvir_sm
        deallocate(C_mo_frag, mo_energies_frag)
        deallocate(S_frag, Hcore_frag, S_inv_sqrt_frag, P_frag)
        


        print '(A,I0,A)', ' Found ', nfrags, ' fragment(s)'
        print *, ''
        ! ------------------------------------------------------------------
        ! Fragment SCF loop (sequential — each activates its own globals)
        ! ------------------------------------------------------------------

        

        do k = 1, nfrags
            read(unit_id, '(A)') frag_dir_name
            frag_path = trim(resolved_mol_dir)//'/'//trim(frag_dir_name)

            print '(A,I0,A,A)', ' --- Fragment ', k, ': ', trim(frag_dir_name)
            print *, '------------------------------------------------------------'

            call init_molecule(mol_frag(k), trim(frag_path))
            call activate_molecule(mol_frag(k), trim(frag_path))

            frag_scf_info(k)%mol = mol_frag(k)
            

            call build_hcore_overlap(mol_frag(k), S_frag, Hcore_frag)
            frag_scf_info(k)%S = S_frag

            call build_s_inv(S_frag, S_inv_sqrt_frag)
            call build_initial_guess(mol_frag(k), Hcore_frag, S_inv_sqrt_frag, P_frag)

            nao_frag  = int(mol_frag(k)%basis%nao)
            nocc_frag = mol_frag(k)%total_electrons / 2
            nvir_frag = nao_frag - nocc_frag

            allocate(C_mo_frag(nao_frag, nao_frag))
            allocate(mo_energies_frag(nao_frag))

            

            call run_scf(mol_frag(k), S_frag, Hcore_frag, S_inv_sqrt_frag, &
                        P_frag, C_mo_frag, mo_energies_frag)

            print '(A,I0,A,F20.10)', ' Fragment ', k, ' converged.'
            print *, ''
            frag_scf_info(k)%C_mo = C_mo_frag(:,:)
            frag_scf_info(k)%mo_energies = mo_energies_frag(:)
            frag_scf_info(k)%nocc = nocc_frag
            frag_scf_info(k)%nvir = nvir_frag

            deallocate(C_mo_frag, mo_energies_frag)
            deallocate(S_frag, Hcore_frag, S_inv_sqrt_frag, P_frag)
        end do
        close(unit_id)

    end subroutine perform_frag_rhf

    

        ! ================================================================
    ! build_ao_map: map fragment AO indices -> supersystem AO indices
    ! ================================================================
    subroutine build_ao_map(mol_sm, mol_frag, ao_map)
        implicit none
        type(molecule), intent(in) :: mol_sm, mol_frag
        integer, allocatable, intent(out) :: ao_map(:)

        integer :: ia_frag, ia_sm, ish_frag, ish_sm
        integer :: ptr_frag, ptr_sm
        integer :: nao_frag, natm_frag, natm_sm, nbas_frag, nbas_sm
        integer :: nao_shl, i, cnt, target_count
        integer :: ia_frag_0based, ia_sm_0based
        integer, allocatable :: atom_map(:), shell_count_sm(:)
        real(c_double) :: dx, dy, dz, dist
        real(c_double), parameter :: TOL = 1.0d-8
        logical :: found

        nao_frag  = int(mol_frag%basis%nao)
        natm_frag = int(mol_frag%basis%natm)
        natm_sm   = int(mol_sm%basis%natm)
        nbas_frag = int(mol_frag%basis%nbas)
        nbas_sm   = int(mol_sm%basis%nbas)

        allocate(ao_map(nao_frag))
        allocate(atom_map(natm_frag))
        allocate(shell_count_sm(natm_sm))
        ao_map = 0
        atom_map = -1
        shell_count_sm = 0

        ! Step 1: Match fragment atoms to supersystem atoms by coordinate
        do ia_frag = 1, natm_frag
            ptr_frag = mol_frag%basis%atm(2, ia_frag) + 1   ! 0-based ptr -> 1-based
            found = .false.
            do ia_sm = 1, natm_sm
                ptr_sm = mol_sm%basis%atm(2, ia_sm) + 1
                dx = mol_frag%basis%env(ptr_frag)   - mol_sm%basis%env(ptr_sm)
                dy = mol_frag%basis%env(ptr_frag+1) - mol_sm%basis%env(ptr_sm+1)
                dz = mol_frag%basis%env(ptr_frag+2) - mol_sm%basis%env(ptr_sm+2)
                dist = sqrt(dx*dx + dy*dy + dz*dz)
                if (dist < TOL) then
                    atom_map(ia_frag) = ia_sm - 1            ! store 0-based
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, 'ERROR: Fragment atom ', ia_frag, ' not found in supersystem'
                stop 1
            end if
        end do

        ! Step 2: For each fragment shell, find the matching supersystem shell.
        !         Shells on the same atom appear in the same order (same basis).
        do ish_frag = 1, nbas_frag
            ia_frag_0based = mol_frag%basis%bas(1, ish_frag)     ! 0-based atom
            ia_sm_0based   = atom_map(ia_frag_0based + 1)        ! mapped sm atom

            target_count = shell_count_sm(ia_sm_0based + 1) + 1
            cnt = 0
            do ish_sm = 1, nbas_sm
                if (mol_sm%basis%bas(1, ish_sm) == ia_sm_0based) then
                    cnt = cnt + 1
                    if (cnt == target_count) then
                        nao_shl = mol_frag%basis%ao_loc(ish_frag+1) &
                                - mol_frag%basis%ao_loc(ish_frag)
                        do i = 0, nao_shl - 1
                            ao_map(mol_frag%basis%ao_loc(ish_frag) + i) = &
                                mol_sm%basis%ao_loc(ish_sm) + i
                        end do
                        shell_count_sm(ia_sm_0based + 1) = target_count
                        exit
                    end if
                end if
            end do
        end do

        deallocate(atom_map, shell_count_sm)

        
    end subroutine build_ao_map

    ! ================================================================
    ! run_rose_localization: fragment-based Pipek-Mezey localization
    !
    ! Builds the coupling matrix  A_k = C_frag_occ_k^T @ S12_k @ C_sm_occ
    ! for each fragment k, stacks them, then calls the existing ROSE
    ! Localization routine from make_ibo.
    ! ================================================================
    subroutine run_rose_localization(frag_scf_info, nfrags, lmo_occ, lmo_vir)
        implicit none
        type(fragment_scf_info), intent(in) :: frag_scf_info(0:)
        integer, intent(in) :: nfrags
        real(c_double), intent(out) :: lmo_occ(:,:), lmo_vir(:,:)

        
        real(c_double), allocatable :: CMO_frag(:,:,:), LMO_frag(:,:,:)
        real(c_double), allocatable :: S12(:,:), temp(:,:), coupling(:,:)
        integer, allocatable :: nmo_frag(:), ao_map(:), nlo_frag(:)
        real(c_double), allocatable :: frag_bias(:)
        integer :: k, nao_sm, nocc_sm, nvir_sm, nao_frag, nocc_k, nvir_k
        integer :: total_frag_occ, total_frag_vir, offset, i

        nao_sm  = size(frag_scf_info(0)%C_mo, 1)
        nocc_sm = int(frag_scf_info(0)%nocc)
        nvir_sm = int(frag_scf_info(0)%nvir)

        ! Count total fragment occupied and virtual MOs
        total_frag_occ = 0
        total_frag_vir = 0
        do k = 1, nfrags
            total_frag_occ = total_frag_occ + int(frag_scf_info(k)%nocc)
            total_frag_vir = total_frag_vir + int(frag_scf_info(k)%nvir)
        end do

        print *, ''
        print *, '=========================================='
        print *, ' ROSE Fragment Localization'
        print *, '=========================================='

        ! ==================================================================
        !  OCCUPIED LOCALIZATION
        ! ==================================================================
        print *, ''
        print *, '--- Occupied space ---'
        print '(A,I0)', '  Supersystem occupied MOs : ', nocc_sm
        print '(A,I0)', '  Total fragment occ MOs   : ', total_frag_occ

        allocate(CMO_frag(total_frag_occ, nocc_sm, 1))
        allocate(LMO_frag(total_frag_occ, nocc_sm, 1))
        allocate(lmo_occ(nocc_sm, nocc_sm))
        allocate(nmo_frag(nfrags))
        CMO_frag = 0.0d0

        offset = 0
        do k = 1, nfrags
            nao_frag = size(frag_scf_info(k)%C_mo, 1)
            nocc_k   = int(frag_scf_info(k)%nocc)
            nmo_frag(k) = nocc_k

            call build_ao_map(frag_scf_info(0)%mol, frag_scf_info(k)%mol, ao_map)

            ! S12 = S_sm(ao_map(:), :)   (nao_frag x nao_sm)
            allocate(S12(nao_frag, nao_sm))
            do i = 1, nao_frag
                S12(i, :) = frag_scf_info(0)%S(ao_map(i), :)
            end do

            ! temp = S12 @ C_sm_occ   (nao_frag x nocc_sm)
            allocate(temp(nao_frag, nocc_sm))
            call dgemm('N', 'N', nao_frag, nocc_sm, nao_sm, 1.0d0, &
                       S12, nao_frag, &
                       frag_scf_info(0)%C_mo, nao_sm, &
                       0.0d0, temp, nao_frag)

            ! coupling = C_frag_occ^T @ temp   (nocc_k x nocc_sm)
            allocate(coupling(nocc_k, nocc_sm))
            call dgemm('T', 'N', nocc_k, nocc_sm, nao_frag, 1.0d0, &
                       frag_scf_info(k)%C_mo, nao_frag, &
                       temp, nao_frag, &
                       0.0d0, coupling, nocc_k)

            CMO_frag(offset+1:offset+nocc_k, :, 1) = coupling(:, :)

            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nocc_k, ' occ MOs'

            offset = offset + nocc_k
            deallocate(ao_map, S12, temp, coupling)
        end do

        call Localization(CMO_frag, LMO_frag, nmo_frag, 2)

        allocate(nlo_frag(nfrags))
        allocate(frag_bias(nfrags))
        frag_bias = 0.0d0
        call sort_on_fragments(LMO_frag, nmo_frag, nlo_frag, frag_bias)

        print *, ''
        print *, 'Occupied LMOs assigned to fragments:'
        do k = 1, nfrags
            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nlo_frag(k), ' LMOs'
        end do
        call print_orbital_populations(CMO_frag, LMO_frag, nmo_frag, nfrags, 'OCC')

        lmo_occ = LMO_frag(:,:,1)

        deallocate(CMO_frag, LMO_frag, nmo_frag, nlo_frag, frag_bias)

        ! ==================================================================
        !  VIRTUAL LOCALIZATION
        ! ==================================================================
        print *, ''
        print *, '--- Virtual space ---'
        print '(A,I0)', '  Supersystem virtual MOs  : ', nvir_sm
        print '(A,I0)', '  Total fragment vir MOs   : ', total_frag_vir

        allocate(CMO_frag(total_frag_vir, nvir_sm, 1))
        allocate(LMO_frag(total_frag_vir, nvir_sm, 1))
        allocate(lmo_vir(nvir_sm, nvir_sm))
        allocate(nmo_frag(nfrags))
        CMO_frag = 0.0d0

        offset = 0
        do k = 1, nfrags
            nao_frag = size(frag_scf_info(k)%C_mo, 1)
            nocc_k   = int(frag_scf_info(k)%nocc)
            nvir_k   = int(frag_scf_info(k)%nvir)
            nmo_frag(k) = nvir_k

            call build_ao_map(frag_scf_info(0)%mol, frag_scf_info(k)%mol, ao_map)

            ! S12 = S_sm(ao_map(:), :)   (nao_frag x nao_sm)
            allocate(S12(nao_frag, nao_sm))
            do i = 1, nao_frag
                S12(i, :) = frag_scf_info(0)%S(ao_map(i), :)
            end do

            ! temp = S12 @ C_sm_vir   (nao_frag x nvir_sm)
            ! C_sm_vir starts at column nocc_sm+1
            allocate(temp(nao_frag, nvir_sm))
            call dgemm('N', 'N', nao_frag, nvir_sm, nao_sm, 1.0d0, &
                       S12, nao_frag, &
                       frag_scf_info(0)%C_mo(1, nocc_sm+1), nao_sm, &
                       0.0d0, temp, nao_frag)

            ! coupling = C_frag_vir^T @ temp   (nvir_k x nvir_sm)
            ! C_frag_vir starts at column nocc_k+1
            allocate(coupling(nvir_k, nvir_sm))
            call dgemm('T', 'N', nvir_k, nvir_sm, nao_frag, 1.0d0, &
                       frag_scf_info(k)%C_mo(1, nocc_k+1), nao_frag, &
                       temp, nao_frag, &
                       0.0d0, coupling, nvir_k)

            CMO_frag(offset+1:offset+nvir_k, :, 1) = coupling(:, :)

            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nvir_k, ' vir MOs'

            offset = offset + nvir_k

            lmo_vir = LMO_frag(:,:,1)
            deallocate(ao_map, S12, temp, coupling)
        end do

        call Localization(CMO_frag, LMO_frag, nmo_frag, 2)

        allocate(nlo_frag(nfrags))
        allocate(frag_bias(nfrags))
        frag_bias = 0.0d0
        call sort_on_fragments(LMO_frag, nmo_frag, nlo_frag, frag_bias)

        print *, ''
        print *, 'Virtual LMOs assigned to fragments:'
        do k = 1, nfrags
            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nlo_frag(k), ' LMOs'
        end do
        call print_orbital_populations(CMO_frag, LMO_frag, nmo_frag, nfrags, 'VIR')

        deallocate(CMO_frag, LMO_frag, nmo_frag, nlo_frag, frag_bias)
    end subroutine run_rose_localization

    ! ================================================================
    ! print_orbital_populations: per-orbital fragment population table
    !
    ! For each LMO j, prints the population on each fragment:
    !   pop(k, j) = sum_i  LMO(offset_k+1 : offset_k+nmo_frag(k), j, 1)^2
    ! A well-localized orbital has ~1.0 on one fragment and ~0.0 elsewhere.
    ! ================================================================
    subroutine print_orbital_populations(CMO, LMO, nmo_frag, nfrags, label)
        implicit none
        real(c_double), intent(in) :: CMO(:,:,:)
        real(c_double), intent(in) :: LMO(:,:,:)
        integer, intent(in) :: nmo_frag(:)
        integer, intent(in) :: nfrags
        character(len=*), intent(in) :: label

        integer :: j, k, nact, imo, nmoi
        real(c_double) :: pop, total_pop
        character(len=256) :: header_fmt, row_fmt

        nact = size(LMO, 2)

        ! Build format strings dynamically based on nfrags
        write(header_fmt, '(A,I0,A)') '(A6,', nfrags, '(A12),A12)'
        write(row_fmt, '(A,I0,A)') '(I6,', nfrags, '(F12.6),F12.6)'

        print *, ''
        write(*, '(A,A,A)') '  Orbital populations (', trim(label), '):'
        ! Print header
        write(*, header_fmt, advance='no') '   LMO'
        do k = 1, nfrags
            write(*, '(A8,I0,A3)', advance='no') '  Frag ', k, '  '
        end do
        write(*, '(A12)') '     Total'

       ! Print each orbital
        do j = 1, nact
            write(*, '(I6)', advance='no') j
            total_pop = 0.0d0
            imo = 0
            do k = 1, nfrags
                nmoi = nmo_frag(k)
                pop = sum(CMO(imo+1:imo+nmoi, j, 1)**2)
                total_pop = total_pop + pop
                write(*, '(F12.6)', advance='no') pop
                imo = imo + nmoi
            end do
            write(*, '(F12.6)') total_pop
        end do

        ! Print each orbital
        do j = 1, nact
            write(*, '(I6)', advance='no') j
            total_pop = 0.0d0
            imo = 0
            do k = 1, nfrags
                nmoi = nmo_frag(k)
                pop = sum(LMO(imo+1:imo+nmoi, j, 1)**2)
                total_pop = total_pop + pop
                write(*, '(F12.6)', advance='no') pop
                imo = imo + nmoi
            end do
            write(*, '(F12.6)') total_pop
        end do
    end subroutine print_orbital_populations


end module rose_interface_module
