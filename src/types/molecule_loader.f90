module molecule_loader
    use iso_c_binding
    use molecule_t, only: molecule

    implicit none

    type(molecule), save :: active_mol
    character(len=1024), save :: active_mol_dir = '.'

    integer :: natm = 0, nbas = 0, env_size = 0, nao = 0
    integer :: total_electrons = 0
    integer(c_int), allocatable :: atm(:), bas(:)
    real(c_double), allocatable :: env(:)
    integer, allocatable :: ao_loc(:)

    contains

    subroutine init_molecule(mol, dir)
        implicit none
        type(molecule), intent(inout) :: mol
        character(len=*), intent(in), optional :: dir
        character(len=1024) :: resolved_dir

        call resolve_molecule_dir(resolved_dir, dir)
        active_mol_dir = trim(resolved_dir)
        call read_molecule_from_dir(mol, trim(resolved_dir))
    end subroutine init_molecule

    subroutine activate_molecule(mol, dir)
        implicit none
        type(molecule), intent(in) :: mol
        character(len=*), intent(in), optional :: dir

        active_mol = mol
        if (present(dir)) then
            active_mol_dir = trim(dir)
        end if
        call sync_compatibility_views()
    end subroutine activate_molecule

    subroutine ensure_molecule_loaded()
        implicit none

        if (.not. allocated(active_mol%basis%env)) then
            call read_molecule_from_dir(active_mol, trim(active_mol_dir))
            call sync_compatibility_views()
        else if (.not. allocated(atm) .or. .not. allocated(bas) .or. .not. allocated(env) .or. .not. allocated(ao_loc)) then
            call sync_compatibility_views()
        end if
    end subroutine ensure_molecule_loaded

    subroutine resolve_molecule_dir(resolved_dir, requested_dir)
        implicit none
        character(len=*), intent(out) :: resolved_dir
        character(len=*), intent(in), optional :: requested_dir

        character(len=1024) :: env_dir
        integer :: env_len, env_stat

        resolved_dir = '.'
        if (present(requested_dir)) then
            if (len_trim(requested_dir) > 0) then
                resolved_dir = trim(requested_dir)
                return
            end if
        end if

        call get_environment_variable('EXAGRAD_MOL_DIR', env_dir, length=env_len, status=env_stat)
        if (env_stat == 0 .and. env_len > 0) then
            resolved_dir = trim(env_dir(1:env_len))
        end if
    end subroutine resolve_molecule_dir

    subroutine read_molecule_from_dir(mol, dir)
        implicit none
        type(molecule), intent(inout) :: mol
        character(len=*), intent(in) :: dir
        integer :: unit_id, i, l, nctr, nao_shl
        integer(c_int), allocatable :: atm_raw(:), bas_raw(:)
        integer :: env_size

        open(newunit=unit_id, file=trim(dir)//'/mol_info.txt', status='old')
        read(unit_id, *) mol%basis%natm
        read(unit_id, *) mol%basis%nbas
        read(unit_id, *) env_size
        read(unit_id, *) mol%basis%nao
        read(unit_id, *) mol%total_electrons
        close(unit_id)

        mol%natoms = int(mol%basis%natm)

        if (allocated(mol%basis%atm)) deallocate(mol%basis%atm)
        if (allocated(mol%basis%bas)) deallocate(mol%basis%bas)
        if (allocated(mol%basis%env)) deallocate(mol%basis%env)
        if (allocated(mol%basis%ao_loc)) deallocate(mol%basis%ao_loc)

        allocate(atm_raw(6 * int(mol%basis%natm)))
        allocate(bas_raw(8 * int(mol%basis%nbas)))
        allocate(mol%basis%env(env_size))

        open(newunit=unit_id, file=trim(dir)//'/ints/atm.txt', status='old')
        read(unit_id, *) atm_raw
        close(unit_id)

        open(newunit=unit_id, file=trim(dir)//'/ints/bas.txt', status='old')
        read(unit_id, *) bas_raw
        close(unit_id)

        open(newunit=unit_id, file=trim(dir)//'/ints/env.txt', status='old')
        read(unit_id, *) mol%basis%env
        close(unit_id)

        allocate(mol%basis%atm(6, int(mol%basis%natm)))
        allocate(mol%basis%bas(8, int(mol%basis%nbas)))
        mol%basis%atm = reshape(atm_raw, shape(mol%basis%atm))
        mol%basis%bas = reshape(bas_raw, shape(mol%basis%bas))
        deallocate(atm_raw, bas_raw)

        allocate(mol%basis%ao_loc(int(mol%basis%nbas) + 1))
        mol%basis%ao_loc(1) = 1_c_int
        do i = 1, int(mol%basis%nbas)
            l = int(mol%basis%bas(2, i))
            nctr = int(mol%basis%bas(4, i))
            nao_shl = (2*l + 1) * nctr
            mol%basis%ao_loc(i+1) = mol%basis%ao_loc(i) + int(nao_shl, c_int)
        end do

        print *, "Molecule initialized:"
        print *, "  Atoms:       ", mol%basis%natm
        print *, "  Shells:      ", mol%basis%nbas
        print *, "  AO functions:", mol%basis%nao
        print *, "  Electrons:   ", mol%total_electrons
    end subroutine read_molecule_from_dir

    subroutine sync_compatibility_views()
        implicit none

        natm = int(active_mol%basis%natm)
        nbas = int(active_mol%basis%nbas)
        nao = int(active_mol%basis%nao)
        total_electrons = active_mol%total_electrons
        env_size = 0
        if (allocated(active_mol%basis%env)) env_size = size(active_mol%basis%env)

        if (allocated(atm)) deallocate(atm)
        if (allocated(bas)) deallocate(bas)
        if (allocated(env)) deallocate(env)
        if (allocated(ao_loc)) deallocate(ao_loc)

        if (allocated(active_mol%basis%atm)) then
            allocate(atm(size(active_mol%basis%atm)))
            atm = reshape(active_mol%basis%atm, [size(active_mol%basis%atm)])
        end if

        if (allocated(active_mol%basis%bas)) then
            allocate(bas(size(active_mol%basis%bas)))
            bas = reshape(active_mol%basis%bas, [size(active_mol%basis%bas)])
        end if

        if (allocated(active_mol%basis%env)) then
            allocate(env(size(active_mol%basis%env)))
            env = active_mol%basis%env
        end if

        if (allocated(active_mol%basis%ao_loc)) then
            allocate(ao_loc(size(active_mol%basis%ao_loc)))
            ao_loc = int(active_mol%basis%ao_loc)
        end if
    end subroutine sync_compatibility_views

end module molecule_loader