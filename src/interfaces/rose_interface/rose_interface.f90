module rose_interface_module

    use iso_c_binding
    use molecule_t,          only: molecule, fragment_scf_info
    use molecule_loader,     only: init_molecule, activate_molecule
    use one_eints,           only: build_hcore_overlap, build_s_inv, build_initial_guess
    use scf_module,          only: run_scf
    use make_iao,            only: MakeIAO_Simple_2013
    use make_ibo,            only: Localization, sort_on_fragments, recanonize, sort_LMO_by_energy
    use linear_algebra,      only: SVD

    implicit none

contains

    subroutine normalize_rose_orbital_ordering(order_in, order_out)
        implicit none

        character(len=*), intent(in) :: order_in
        character(len=*), intent(out) :: order_out

        character(len=64) :: lowered

        lowered = lowercase_copy(trim(adjustl(order_in)))

        select case (trim(lowered))
        case ('', 'default', 'fragment', 'fragments', 'current')
            order_out = 'fragment'
        case ('fragment_recanon', 'recanon', 'recanonicalize', 'fragment_semicanonical', 'semicanonical')
            order_out = 'fragment_recanon'
        case ('energy', 'global_energy', 'energy_sorted', 'sorted_energy')
            order_out = 'energy'
        case default
            print *, 'WARNING: Unknown EXAGRAD_ROSE_ORBITAL_ORDER = ', trim(order_in)
            print *, '         Falling back to fragment ordering.'
            order_out = 'fragment'
        end select
    end subroutine normalize_rose_orbital_ordering

    pure function lowercase_copy(text) result(lowered)
        implicit none

        character(len=*), intent(in) :: text
        character(len=len(text)) :: lowered
        integer :: idx, code

        lowered = text
        do idx = 1, len(text)
            code = iachar(text(idx:idx))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lowered(idx:idx) = achar(code + 32)
            end if
        end do
    end function lowercase_copy

    subroutine rotate_localized_iao(cmo_iao, lmo_rot, lmo_iao)
        implicit none

        real(c_double), intent(in) :: cmo_iao(:,:,:), lmo_rot(:,:,:)
        real(c_double), intent(inout) :: lmo_iao(:,:,:)

        integer :: nref, nact

        nref = size(cmo_iao, 1)
        nact = size(cmo_iao, 2)

        call dgemm('N', 'N', nref, nact, nact, 1.0d0, cmo_iao(:,:,1), nref, lmo_rot(:,:,1), nact, 0.0d0, lmo_iao(:,:,1), nref)
    end subroutine rotate_localized_iao

    subroutine apply_localization_ordering(cmo_iao, lmo_iao, lmo_rot, cmo_energies, nlo_frag, ordering_mode, space_label)
        implicit none

        real(c_double), intent(in) :: cmo_iao(:,:,:), cmo_energies(:)
        real(c_double), intent(inout) :: lmo_iao(:,:,:)
        real(c_double), intent(inout) :: lmo_rot(:,:,:)
        integer, intent(in) :: nlo_frag(:)
        character(len=*), intent(in) :: ordering_mode, space_label

        real(c_double), allocatable :: lmo_energies(:)
        character(len=32) :: normalized_mode
        integer :: nact

        nact = size(lmo_iao, 2)
        if (nact == 0) return

        call normalize_rose_orbital_ordering(ordering_mode, normalized_mode)
        if (trim(normalized_mode) == 'fragment') return

        allocate(lmo_energies(nact))

        select case (trim(normalized_mode))
        case ('fragment_recanon')
            call recanonize(lmo_rot, cmo_energies, lmo_energies, nlo_frag)
            call rotate_localized_iao(cmo_iao, lmo_rot, lmo_iao)
            print '(A,A,A)', '  ', trim(space_label), ' LMOs recanonicalized within fragment blocks.'
        case ('energy')
            call recanonize(lmo_rot, cmo_energies, lmo_energies, nlo_frag)
            call sort_LMO_by_energy(lmo_rot, cmo_energies, lmo_energies, .true.)
            call rotate_localized_iao(cmo_iao, lmo_rot, lmo_iao)
            print '(A,A,A)', '  ', trim(space_label), ' LMOs globally sorted by semicanonical energy.'
        end select

        deallocate(lmo_energies)
    end subroutine apply_localization_ordering

    subroutine compute_lmo_energies(lmo_rot, cmo_energies, lmo_energies)
        implicit none

        real(c_double), intent(in) :: lmo_rot(:,:,:)
        real(c_double), intent(in) :: cmo_energies(:)
        real(c_double), allocatable, intent(out) :: lmo_energies(:)

        real(c_double), allocatable :: hf_fock(:,:), temp(:,:), f_lmo(:,:)
        integer :: j, nact

        nact = size(lmo_rot, 2)
        allocate(lmo_energies(nact))
        allocate(hf_fock(nact, nact), temp(nact, nact), f_lmo(nact, nact))

        hf_fock = 0.0d0
        do j = 1, nact
            hf_fock(j, j) = cmo_energies(j)
        end do

        call dgemm('N', 'N', nact, nact, nact, 1.0d0, hf_fock, nact, lmo_rot(:,:,1), nact, 0.0d0, temp, nact)
        call dgemm('T', 'N', nact, nact, nact, 1.0d0, lmo_rot(:,:,1), nact, temp, nact, 0.0d0, f_lmo, nact)

        do j = 1, nact
            lmo_energies(j) = f_lmo(j, j)
        end do

        deallocate(hf_fock, temp, f_lmo)
    end subroutine compute_lmo_energies

    subroutine perform_frag_rhf(mol_dir, nfrags, mol_frag, frag_scf_info)

        implicit none

        character(len=*), intent(in) :: mol_dir
        integer, intent(out) :: nfrags
        type(molecule), allocatable, intent(out) :: mol_frag(:)
        type(fragment_scf_info), allocatable, intent(out) :: frag_scf_info(:)
        real(c_double), allocatable :: C_mo_frag(:,:), mo_energies_frag(:)
        real(c_double), allocatable :: S_inv_sqrt_frag(:,:), P_frag(:,:), Hcore_frag(:,:),S_frag(:,:)
        integer :: nao_frag, nocc_frag, nvir_frag, nao_sm, nocc_sm, nvir_sm, k
        integer :: env_len, env_stat, unit_id, ios, n_val_tmp
        character(len=1024) :: frag_dir_name
        character(len=2048) :: frag_path, line_buf
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
        call activate_molecule(mol_frag(0))
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
            ! Read fragment dir; optionally followed by n_val_vir on the same line:
            !   frag_1          <- use all fragment virtuals as valence (default)
            !   frag_1  3       <- use only first 3 fragment virtuals as valence
            read(unit_id, '(A)') line_buf
            n_val_tmp = -1
            read(line_buf, *, iostat=ios) frag_dir_name, n_val_tmp
           
            if (ios /= 0) read(line_buf, *, iostat=ios) frag_dir_name
            frag_path = trim(resolved_mol_dir)//'/'//trim(frag_dir_name)

        

            print '(A,I0,A,A)', ' --- Fragment ', k, ': ', trim(frag_dir_name), ' ---  n_val_vir = ', n_val_tmp
            print *, '------------------------------------------------------------'

            call init_molecule(mol_frag(k), trim(frag_path))
            call activate_molecule(mol_frag(k))

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

            print '(A,I0,A)', ' Fragment ', k, ' converged.'
            print *, ''
            frag_scf_info(k)%C_mo = C_mo_frag(:,:)
            frag_scf_info(k)%mo_energies = mo_energies_frag(:)
            frag_scf_info(k)%nocc = nocc_frag
            frag_scf_info(k)%nvir = nvir_frag
            frag_scf_info(k)%n_val_vir = n_val_tmp   ! -1 = all virtuals are valence

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

    subroutine build_embedded_fragment_space(frag_scf_info, nfrags, occupied_space, nmo_frag, P12, label)
        implicit none
        type(fragment_scf_info), intent(in) :: frag_scf_info(0:)
        integer, intent(in) :: nfrags
        logical, intent(in) :: occupied_space
        integer, allocatable, intent(out) :: nmo_frag(:)
        real(c_double), allocatable, intent(out) :: P12(:,:)
        character(len=*), intent(in) :: label

        integer, allocatable :: ao_map(:)
        integer :: k, i, nao_sm, nao_frag, nocc_k, nact_k, total_fragment_functions, offset, col_start

        nao_sm = size(frag_scf_info(0)%C_mo, 1)
        total_fragment_functions = 0
        allocate(nmo_frag(nfrags))

        do k = 1, nfrags
            nocc_k = int(frag_scf_info(k)%nocc)
            if (occupied_space) then
                nact_k = nocc_k
            else
                ! Respect n_val_vir: only the first n_val_vir virtual MOs are valence refs.
                ! -1 (default) means use all fragment virtual MOs.
                if (frag_scf_info(k)%n_val_vir >= 0) then
                    nact_k = int(frag_scf_info(k)%n_val_vir)
                else
                    nact_k = int(frag_scf_info(k)%nvir)
                end if
            end if
            nmo_frag(k) = nact_k
            total_fragment_functions = total_fragment_functions + nact_k
        end do

        allocate(P12(nao_sm, total_fragment_functions))
        P12 = 0.0d0

        offset = 0
        do k = 1, nfrags
            nao_frag = size(frag_scf_info(k)%C_mo, 1)
            nocc_k = int(frag_scf_info(k)%nocc)
            nact_k = nmo_frag(k)
            if (nact_k == 0) cycle
            if (occupied_space) then
                col_start = 1
            else
                col_start = nocc_k + 1
            end if

            call build_ao_map(frag_scf_info(0)%mol, frag_scf_info(k)%mol, ao_map)
            do i = 1, nao_frag
                P12(ao_map(i), offset+1:offset+nact_k) = &
                    frag_scf_info(k)%C_mo(i, col_start:col_start+nact_k-1)
            end do

            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nact_k, ' ', trim(label)
            offset = offset + nact_k
            deallocate(ao_map)
        end do
    end subroutine build_embedded_fragment_space

    subroutine build_iao_coefficients(C_block, S11, P12, CMO_IAO)
        implicit none
        real(c_double), intent(in) :: C_block(:,:), S11(:,:), P12(:,:)
        real(c_double), allocatable, intent(out) :: CMO_IAO(:,:,:)

        real(c_double), allocatable :: A(:,:,:), C3(:,:,:), S11_3(:,:,:)
        real(c_double), allocatable :: S12(:,:), S21(:,:), S22(:,:)
        real(c_double), allocatable :: S12_3(:,:,:), S21_3(:,:,:), S22_3(:,:,:), P12_3(:,:,:)
        real(c_double), allocatable :: temp(:,:)
        integer :: nao_sm, nref, nact

        nao_sm = size(S11, 1)
        nref = size(P12, 2)
        nact = size(C_block, 2)

        allocate(S12(nao_sm, nref), S21(nref, nao_sm), S22(nref, nref))
        call dgemm('N', 'N', nao_sm, nref, nao_sm, 1.0d0, S11, nao_sm, P12, nao_sm, 0.0d0, S12, nao_sm)
        call dgemm('T', 'N', nref, nao_sm, nao_sm, 1.0d0, P12, nao_sm, S11, nao_sm, 0.0d0, S21, nref)
        call dgemm('N', 'N', nref, nref, nao_sm, 1.0d0, S21, nref, P12, nao_sm, 0.0d0, S22, nref)

        allocate(A(nao_sm, nref, 1), C3(nao_sm, nact, 1), S11_3(nao_sm, nao_sm, 1))
        allocate(S12_3(nao_sm, nref, 1), S21_3(nref, nao_sm, 1), S22_3(nref, nref, 1), P12_3(nao_sm, nref, 1))

        C3(:,:,1) = C_block
        S11_3(:,:,1) = S11
        S12_3(:,:,1) = S12
        S21_3(:,:,1) = S21
        S22_3(:,:,1) = S22
        P12_3(:,:,1) = P12

        call MakeIAO_Simple_2013(A, C3, S11_3, S12_3, S21_3, S22_3, P12_3)

        allocate(temp(nao_sm, nact))
        allocate(CMO_IAO(nref, nact, 1))
        call dgemm('N', 'N', nao_sm, nact, nao_sm, 1.0d0, S11, nao_sm, C_block, nao_sm, 0.0d0, temp, nao_sm)
        call dgemm('T', 'N', nref, nact, nao_sm, 1.0d0, A(:,:,1), nao_sm, temp, nao_sm, 0.0d0, CMO_IAO(:,:,1), nref)

        deallocate(A, C3, S11_3, S12, S21, S22, S12_3, S21_3, S22_3, P12_3, temp)
    end subroutine build_iao_coefficients

    ! ================================================================
    ! build_ico_rotation: complementary space of the valence IAO projection
    !
    ! Given CMO_IAO (nref × nvir × 1) — virtual CMOs in the IAO basis —
    ! returns ICO_rot (nvir × nvir_hard), the columns of the right null
    ! space of CMO_IAO obtained via SVD of CMO_IAO^T.
    ! Reuses the SVD routine already present in linear_algebra.
    ! ================================================================
    subroutine build_ico_rotation(CMO_IAO, ICO_rot, nvir_hard)
        implicit none

        real(c_double), intent(in)               :: CMO_IAO(:,:,:)
        real(c_double), allocatable, intent(out) :: ICO_rot(:,:)
        integer,        intent(out)              :: nvir_hard

        real(c_double), allocatable :: CIAO_T(:,:,:), U(:,:,:), VT(:,:,:)
        real(c_double), allocatable :: svals(:)
        integer :: nref, nvir

        nref = size(CMO_IAO, 1)
        nvir = size(CMO_IAO, 2)
        nvir_hard = nvir - nref

        if (nvir_hard <= 0) then
            print '(A)', '  ICO: no complementary space (nref >= nvir), skipping hard virtuals'
            nvir_hard = 0
            return
        end if

        ! Transpose to (nvir × nref × 1) for SVD(A, svals, U, VT)
        ! After SVD: U is (nvir × nvir); null space = columns nref+1..nvir of U
        allocate(CIAO_T(nvir, nref, 1))
        CIAO_T(:,:,1) = transpose(CMO_IAO(:,:,1))
        allocate(U(nvir, nvir, 1), VT(nref, nref, 1), svals(nref))

        call SVD(CIAO_T, svals, U, VT)

        allocate(ICO_rot(nvir, nvir_hard))
        ICO_rot = U(:, nref+1:nvir, 1)

        deallocate(CIAO_T, U, VT, svals)
    end subroutine build_ico_rotation

    ! ================================================================
    ! build_hard_vir_fragment_space: P12 for hard virtual references
    !
    ! For each fragment k, uses fragment virtual MOs beyond the valence
    ! block: columns nocc_k + n_val_vir_k + 1 .. nocc_k + nvir_k.
    ! n_val_vir = -1 (default) → no hard virtual refs for that fragment.
    ! Reuses build_ao_map from this module.
    ! ================================================================
    subroutine build_hard_vir_fragment_space(frag_scf_info, nfrags, nmo_frag_hard, P12_hard, label)
        implicit none

        type(fragment_scf_info), intent(in)      :: frag_scf_info(0:)
        integer,                 intent(in)      :: nfrags
        integer, allocatable,    intent(out)     :: nmo_frag_hard(:)
        real(c_double), allocatable, intent(out) :: P12_hard(:,:)
        character(len=*),        intent(in)      :: label

        integer, allocatable :: ao_map(:)
        integer :: k, i, nao_sm, nao_frag, nocc_k, nvir_k, n_val_k, n_hard_k, total_hard, offset, col_start

        nao_sm     = size(frag_scf_info(0)%C_mo, 1)
        total_hard = 0
        allocate(nmo_frag_hard(nfrags))

        do k = 1, nfrags
            nocc_k = int(frag_scf_info(k)%nocc)
            nvir_k = int(frag_scf_info(k)%nvir)
            if (frag_scf_info(k)%n_val_vir >= 0) then
                n_val_k = int(frag_scf_info(k)%n_val_vir)
            else
                n_val_k = nvir_k   ! default: all virtuals are valence
            end if
            nmo_frag_hard(k) = max(0, nvir_k - n_val_k)
            total_hard = total_hard + nmo_frag_hard(k)
        end do

        allocate(P12_hard(nao_sm, total_hard))
        P12_hard = 0.0d0
        if (total_hard == 0) return

        offset = 0
        do k = 1, nfrags
            n_hard_k = nmo_frag_hard(k)
            if (n_hard_k == 0) cycle

            nocc_k   = int(frag_scf_info(k)%nocc)
            n_val_k  = int(frag_scf_info(k)%n_val_vir)   ! guaranteed >= 0 here
            nao_frag = size(frag_scf_info(k)%C_mo, 1)
            col_start = nocc_k + n_val_k + 1

            call build_ao_map(frag_scf_info(0)%mol, frag_scf_info(k)%mol, ao_map)
            do i = 1, nao_frag
                P12_hard(ao_map(i), offset+1:offset+n_hard_k) = &
                    frag_scf_info(k)%C_mo(i, col_start:col_start+n_hard_k-1)
            end do
            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', n_hard_k, ' '//trim(label)
            offset = offset + n_hard_k
            deallocate(ao_map)
        end do
    end subroutine build_hard_vir_fragment_space

    ! ================================================================
    ! run_rose_localization: fragment-based Pipek-Mezey localization
    !
    ! Builds the coupling matrix  A_k = C_frag_occ_k^T @ S12_k @ C_sm_occ
    ! for each fragment k, stacks them, then calls the existing ROSE
    ! Localization routine from make_ibo.
    ! ================================================================
    subroutine run_rose_localization(frag_scf_info, nfrags, nlo_frag_occ, nlo_frag_vir, c_lmo_occ, U_occ, c_lmo_vir, U_vir, ordering_mode, lmo_vir_hard)
        implicit none
        type(fragment_scf_info), intent(in) :: frag_scf_info(0:)
        integer, intent(in) :: nfrags
        integer, allocatable, intent(out) ::  nlo_frag_occ(:), nlo_frag_vir(:)
        real(c_double), allocatable, intent(out) :: c_lmo_occ(:,:), c_lmo_vir(:,:), U_occ(:,:), U_vir(:,:)
        ! lmo_vir_hard: localized hard virtual MOs in the IVO reference basis (nref_hard × nvir_hard).
        ! Only allocated when hard virtual references are present (any n_val_vir >= 0 and < nvir).
        real(c_double), allocatable, intent(out), optional :: lmo_vir_hard(:,:)
        character(len=*), intent(in), optional :: ordering_mode

        real(c_double), allocatable :: CMO_frag(:,:,:), LMO_frag(:,:,:), LMO_rot(:,:,:)
        real(c_double), allocatable :: P12(:,:)
        integer, allocatable :: nmo_frag(:), nlo_frag(:)
        real(c_double), allocatable :: frag_bias(:)
        ! Hard virtual workspace
        real(c_double), allocatable :: ICO_rot(:,:), C_ico(:,:), P12_hard(:,:)
        real(c_double), allocatable :: CMO_IVO(:,:,:), LMO_IVO(:,:,:), LMO_rot_hard(:,:,:)
        real(c_double), allocatable :: ico_energies(:)
        integer, allocatable :: nmo_frag_hard(:), nlo_frag_hard(:)
        integer :: nvir_hard, total_hard_refs
        logical :: has_hard_vir
        integer :: k, nao_sm, nocc_sm, nvir_sm
        character(len=32) :: active_ordering

        nao_sm  = size(frag_scf_info(0)%C_mo, 1)
        nocc_sm = int(frag_scf_info(0)%nocc)
        nvir_sm = int(frag_scf_info(0)%nvir)

        active_ordering = 'fragment'
        if (present(ordering_mode)) call normalize_rose_orbital_ordering(ordering_mode, active_ordering)

        print *, ''
        print *, '=========================================='
        print *, ' ROSE Fragment Localization'
        print *, '=========================================='
        print '(A,A)', ' Orbital ordering mode      : ', trim(active_ordering)

        ! ==================================================================
        !  OCCUPIED LOCALIZATION
        ! ==================================================================
        print *, ''
        print *, '--- Occupied space ---'
        print '(A,I0)', '  Supersystem occupied MOs : ', nocc_sm
        call build_embedded_fragment_space(frag_scf_info, nfrags, .true., nmo_frag, P12, 'occ reference functions')
        print '(A,I0)', '  Total fragment occ refs  : ', size(P12, 2)

        call build_iao_coefficients(frag_scf_info(0)%C_mo(:, 1:nocc_sm), frag_scf_info(0)%S, P12, CMO_frag)
        allocate(LMO_frag(size(CMO_frag, 1), size(CMO_frag, 2), 1))
        allocate(LMO_rot(nocc_sm, nocc_sm, 1))

        call Localization(CMO_frag, LMO_frag, nmo_frag, 2, LMO_rot)

        allocate(nlo_frag(nfrags))
        allocate(nlo_frag_occ(nfrags))
        allocate(frag_bias(nfrags))
        frag_bias = 0.0d0
        call sort_on_fragments(LMO_frag, nmo_frag, nlo_frag, frag_bias, LMO_rot)
        call apply_localization_ordering(CMO_frag, LMO_frag, LMO_rot, frag_scf_info(0)%mo_energies(1:nocc_sm), nlo_frag, active_ordering, 'Occupied')

        print *, ''
        print *, 'Occupied LMOs assigned to fragments:'
        do k = 1, nfrags
            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nlo_frag(k), ' LMOs'
        end do
        call print_orbital_populations(CMO_frag, LMO_frag, LMO_rot, frag_scf_info(0)%mo_energies(1:nocc_sm), nmo_frag, nfrags, 'OCC')

        ! Performa an C_mo -> C

        ! C_lmo_occ (nao x nocc) = C_mo(:,1:nocc) @ U_occ
        allocate(c_lmo_occ(nao_sm, nocc_sm))
        call dgemm('N', 'N', nao_sm, nocc_sm, nocc_sm, 1.0d0, frag_scf_info(0)%C_mo(:, 1:nocc_sm), nao_sm, &
                   LMO_rot(:,:,1), nocc_sm, 0.0d0, c_lmo_occ, nao_sm)

        U_occ = LMO_rot(:,:,1)
   
        nlo_frag_occ = nlo_frag

        deallocate(CMO_frag, LMO_frag, LMO_rot,nlo_frag, nmo_frag, frag_bias, P12)

        ! ==================================================================
        !  VALENCE VIRTUAL LOCALIZATION
        ! ==================================================================
        print *, ''
        print *, '--- Virtual space (valence) ---'
        print '(A,I0)', '  Supersystem virtual MOs  : ', nvir_sm
        call build_embedded_fragment_space(frag_scf_info, nfrags, .false., nmo_frag, P12, 'valence vir reference functions')
        print '(A,I0)', '  Total fragment vir refs  : ', size(P12, 2)

        call build_iao_coefficients(frag_scf_info(0)%C_mo(:, nocc_sm+1:nocc_sm+nvir_sm), frag_scf_info(0)%S, P12, CMO_frag)
        allocate(LMO_frag(size(CMO_frag, 1), size(CMO_frag, 2), 1))
        allocate(LMO_rot(nvir_sm, nvir_sm, 1))

        call Localization(CMO_frag, LMO_frag, nmo_frag, 2, LMO_rot)

        allocate(nlo_frag(nfrags))
        allocate(nlo_frag_vir(nfrags))
        allocate(frag_bias(nfrags))
        frag_bias = 0.0d0
        call sort_on_fragments(LMO_frag, nmo_frag, nlo_frag, frag_bias, LMO_rot)
        call apply_localization_ordering(CMO_frag, LMO_frag, LMO_rot, frag_scf_info(0)%mo_energies(nocc_sm+1:nocc_sm+nvir_sm), nlo_frag, active_ordering, 'Virtual')

        print *, ''
        print *, 'Valence virtual LMOs assigned to fragments:'
        do k = 1, nfrags
            print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nlo_frag(k), ' LMOs'
        end do
        call print_orbital_populations(CMO_frag, LMO_frag, LMO_rot, frag_scf_info(0)%mo_energies(nocc_sm+1:nocc_sm+nvir_sm), nmo_frag, nfrags, 'VIR_VAL')

    
        allocate(U_vir(size(LMO_rot, 1), size(LMO_rot, 2)))

        allocate(c_lmo_vir(nao_sm, nvir_sm))
        call dgemm('N', 'N', nao_sm, nvir_sm, nvir_sm, 1.0d0, frag_scf_info(0)%C_mo(:, nocc_sm+1:nocc_sm+nvir_sm), nao_sm, &
                   LMO_rot(:,:,1), nvir_sm, 0.0d0, c_lmo_vir, nao_sm)
        U_vir = LMO_rot(:,:,1)


      
        ! ==================================================================
        !  HARD VIRTUAL LOCALIZATION  (ICO / IVO step from original ROSE)
        !  Reuses: build_iao_coefficients, Localization, sort_on_fragments,
        !          recanonize — only build_ico_rotation is new.
        ! ==================================================================

        ! Check whether any fragment exposes hard virtual references
        has_hard_vir = .false.
        do k = 1, nfrags
            if (frag_scf_info(k)%n_val_vir >= 0 .and. &
                frag_scf_info(k)%n_val_vir < frag_scf_info(k)%nvir) then
                has_hard_vir = .true.
                exit
            end if
        end do

        if (has_hard_vir) then
            print *, ''
            print *, '--- Virtual space (hard) ---'

            ! 1. Build ICO: null space of valence IAO projection in virtual MO space.
            !    CMO_frag (nref_vir × nvir_sm) from above; ICO_rot (nvir_sm × nvir_hard).
            call build_ico_rotation(CMO_frag, ICO_rot, nvir_hard)

            if (nvir_hard > 0) then
                print '(A,I0)', '  Complementary (hard) virtual MOs : ', nvir_hard

                ! 2. ICO MOs in AO space: C_ico = C_vir @ ICO_rot
                nao_sm = size(frag_scf_info(0)%C_mo, 1)
                allocate(C_ico(nao_sm, nvir_hard))
                call dgemm('N', 'N', nao_sm, nvir_hard, nvir_sm, 1.0d0, &
                           frag_scf_info(0)%C_mo(:, nocc_sm+1), nao_sm, &
                           ICO_rot, nvir_sm, 0.0d0, C_ico, nao_sm)

                ! 3. Hard virtual reference functions from fragments
                call build_hard_vir_fragment_space(frag_scf_info, nfrags, nmo_frag_hard, P12_hard, 'hard vir reference functions')
                total_hard_refs = size(P12_hard, 2)
                print '(A,I0)', '  Total fragment hard vir refs     : ', total_hard_refs

                ! Original ROSE requires: nref_hard > nvir_hard_actual (more refs than virtuals).
                ! Cap nvir_hard_actual = min(nvir_hard, total_hard_refs - 1) for numerical stability.
                if (total_hard_refs >= 2) then
                    ! Truncate C_ico / ICO_rot to nvir_hard_actual columns (lowest-energy ICO MOs).
                    ! Here we just take the first columns; they come from SVD ordered by singular
                    ! value magnitude, so the complement near-zero columns are already ordered.
                    block
                        integer :: nvir_hard_actual
                        real(c_double), allocatable :: C_ico_trunc(:,:), ICO_rot_trunc(:,:)

                        nvir_hard_actual = min(nvir_hard, total_hard_refs - 1)
                        if (nvir_hard_actual < nvir_hard) then
                            print '(A,I0,A,I0,A)', '  Capping hard virtuals to ', nvir_hard_actual, &
                                ' (need nref > nvir; nref=', total_hard_refs, ')'
                            allocate(C_ico_trunc(nao_sm, nvir_hard_actual))
                            allocate(ICO_rot_trunc(nvir_sm, nvir_hard_actual))
                            C_ico_trunc   = C_ico(:, 1:nvir_hard_actual)
                            ICO_rot_trunc = ICO_rot(:, 1:nvir_hard_actual)
                        else
                            allocate(C_ico_trunc(nao_sm, nvir_hard_actual))
                            allocate(ICO_rot_trunc(nvir_sm, nvir_hard_actual))
                            C_ico_trunc   = C_ico
                            ICO_rot_trunc = ICO_rot
                        end if

                    ! 4. IVO construction: reuse build_iao_coefficients in the ICO subspace
                    call build_iao_coefficients(C_ico_trunc, frag_scf_info(0)%S, P12_hard, CMO_IVO)

                    ! 5. Localize hard virtuals (reuses Localization from make_ibo)
                    allocate(LMO_IVO(size(CMO_IVO, 1), size(CMO_IVO, 2), 1))
                    allocate(LMO_rot_hard(nvir_hard_actual, nvir_hard_actual, 1))
                    call Localization(CMO_IVO, LMO_IVO, nmo_frag_hard, 2, LMO_rot_hard)

                    ! 6. Sort by fragments (reuses sort_on_fragments from make_ibo)
                    allocate(nlo_frag_hard(nfrags))
                    frag_bias = 0.0d0
                    call sort_on_fragments(LMO_IVO, nmo_frag_hard, nlo_frag_hard, frag_bias, LMO_rot_hard)

                    ! 7. Recanonize using approximate ICO orbital energies:
                    !    eps_ico(i) = sum_j  ICO_rot_trunc(j,i)^2 * eps_vir(j)
                    allocate(ico_energies(nvir_hard_actual))
                    do k = 1, nvir_hard_actual
                        ico_energies(k) = sum(ICO_rot_trunc(:, k)**2 * &
                            frag_scf_info(0)%mo_energies(nocc_sm+1:nocc_sm+nvir_sm))
                    end do
                    call recanonize(LMO_rot_hard, ico_energies, ico_energies, nlo_frag_hard)

                    print *, ''
                    print *, 'Hard virtual LMOs assigned to fragments:'
                    do k = 1, nfrags
                        print '(A,I0,A,I0,A)', '  Fragment ', k, ': ', nlo_frag_hard(k), ' LMOs'
                    end do
                    call print_orbital_populations(CMO_IVO, LMO_IVO, LMO_rot_hard, ico_energies, nmo_frag_hard, nfrags, 'VIR_HARD')

                    if (present(lmo_vir_hard)) then
                        allocate(lmo_vir_hard(size(LMO_IVO, 1), size(LMO_IVO, 2)))
                        lmo_vir_hard = LMO_IVO(:,:,1)
                    end if

                    deallocate(LMO_IVO, LMO_rot_hard, nlo_frag_hard, ico_energies)
                    deallocate(CMO_IVO, C_ico_trunc, ICO_rot_trunc)
                    end block
                else
                    print '(A,I0,A)', '  WARNING: too few hard virtual refs (', total_hard_refs, '), skipping hard virtual localization'
                    if (present(lmo_vir_hard)) allocate(lmo_vir_hard(0, 0))
                end if

                deallocate(C_ico, P12_hard, nmo_frag_hard, ICO_rot)
            end if
        else
            if (present(lmo_vir_hard)) allocate(lmo_vir_hard(0, 0))
        end if

        nlo_frag_vir = nlo_frag

        deallocate(CMO_frag, LMO_frag, LMO_rot, nmo_frag, nlo_frag, frag_bias, P12)
    end subroutine run_rose_localization

    ! ================================================================
    ! print_orbital_populations: per-orbital fragment population table
    !
    ! For each LMO j, prints the population on each fragment:
    !   pop(k, j) = sum_i  LMO(offset_k+1 : offset_k+nmo_frag(k), j, 1)^2
    ! A well-localized orbital has ~1.0 on one fragment and ~0.0 elsewhere.
    ! ================================================================
    subroutine print_orbital_populations(CMO, LMO, LMO_rot, cmo_energies, nmo_frag, nfrags, label)
        implicit none
        real(c_double), intent(in) :: CMO(:,:,:)
        real(c_double), intent(in) :: LMO(:,:,:)
        real(c_double), intent(in) :: LMO_rot(:,:,:)
        real(c_double), intent(in) :: cmo_energies(:)
        integer, intent(in) :: nmo_frag(:)
        integer, intent(in) :: nfrags
        character(len=*), intent(in) :: label

        integer :: j, k, nact, imo, nmoi
        real(c_double) :: pop, total_pop
        real(c_double), allocatable :: lmo_energies(:)
        real(c_double), parameter :: HARTREE_TO_EV = 27.211386245988d0
        character(len=256) :: header_fmt

        nact = size(LMO, 2)
        call compute_lmo_energies(LMO_rot, cmo_energies, lmo_energies)

        ! Build format strings dynamically based on nfrags
        write(header_fmt, '(A,I0,A)') '(A6,A14,', nfrags, '(A12),A12)'

        print *, ''
        write(*, '(A,A,A)') '  Orbital populations (', trim(label), '):'
        print *, '   Reference orbitals:'
        write(*, header_fmt, advance='no') '   LMO'
        write(*, '(A14)', advance='no') ' Energy (eV)'
        do k = 1, nfrags
            write(*, '(A8,I0,A3)', advance='no') '  Frag ', k, '  '
        end do
        write(*, '(A12)') '     Total'

       ! Print each reference orbital
        do j = 1, nact
            write(*, '(I6)', advance='no') j
            write(*, '(F14.6)', advance='no') cmo_energies(j) * HARTREE_TO_EV
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

        print *, ''
        print *, '   Final localized orbitals:'
        write(*, header_fmt, advance='no') '   LMO'
        write(*, '(A14)', advance='no') ' Energy (eV)'
        do k = 1, nfrags
            write(*, '(A8,I0,A3)', advance='no') '  Frag ', k, '  '
        end do
        write(*, '(A12)') '     Total'

        ! Print each final localized orbital
        do j = 1, nact
            write(*, '(I6)', advance='no') j
            write(*, '(F14.6)', advance='no') lmo_energies(j) * HARTREE_TO_EV
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

        deallocate(lmo_energies)
    end subroutine print_orbital_populations


end module rose_interface_module
