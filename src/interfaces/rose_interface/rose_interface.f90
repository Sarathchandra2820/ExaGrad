module rose_interface_module

    use iso_c_binding
    use molecule_t,          only: molecule, fragment_scf_info
    use molecule_loader,     only: init_molecule, activate_molecule
    use one_eints,           only: build_hcore_overlap, build_s_inv, build_initial_guess
    use scf_module,          only: run_scf
    use fock_builder_module,  only: normalize_fock_method



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
        integer :: nao_frag, nocc_frag, nvir_frag, k
        character(len=64) :: method_env
        character(len=32) :: fock_method
        integer :: env_len, env_stat, unit_id, ios
        character(len=1024) :: frag_dir_name
        character(len=2048) :: frag_path



        ! ------------------------------------------------------------------
        ! Determine molecule directory from environment
        ! ------------------------------------------------------------------
        ! mol_dir = '.'
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
        allocate(frag_scf_info(nfrags))

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

    end subroutine perform_rhf_frag

    ! subroutine calculate_overlap(frag_scf_info, S11, S22, S12)
    !     implicit none 

    !     type(fragment_scf_info), intent(in) :: frag_scf_info(:)
    !     real(c_double), allocatable, intent(out) :: S11(:,:), S22(:,:), S12(:,:)

    !     !temps 
    !     real(c_double), allocatable :: temp1(:,:), temp2(:,:)

    !     nao = size(S11,1)


        


    ! end subroutine calculate_overlap


end module rose_interface_module
