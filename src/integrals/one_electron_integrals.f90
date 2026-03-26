module one_eints
    use, intrinsic :: iso_c_binding
    use libcint_interface
    use math_utils
    use molecule_t, only: molecule
    implicit none

    type(molecule), save :: active_mol

    integer :: natm = 0, nbas = 0, env_size = 0, nao = 0
    integer :: total_electrons = 0
    integer(c_int), allocatable :: atm(:), bas(:)
    real(c_double), allocatable :: env(:)
    integer, allocatable :: ao_loc(:)

contains

    subroutine init_molecule(mol)
        implicit none
        type(molecule), intent(inout), optional :: mol
        integer :: unit_id, i, l, nctr, nao_shl
        integer(c_int), allocatable :: atm_raw(:), bas_raw(:)

        open(newunit=unit_id, file='ints/mol_info.txt', status='old')
        read(unit_id, *) active_mol%basis%natm
        read(unit_id, *) active_mol%basis%nbas
        read(unit_id, *) env_size
        read(unit_id, *) active_mol%basis%nao
        read(unit_id, *) active_mol%total_electrons
        close(unit_id)

        active_mol%natoms = int(active_mol%basis%natm)

        if (allocated(atm)) deallocate(atm)
        if (allocated(bas)) deallocate(bas)
        if (allocated(env)) deallocate(env)
        if (allocated(ao_loc)) deallocate(ao_loc)

        if (allocated(active_mol%basis%atm)) deallocate(active_mol%basis%atm)
        if (allocated(active_mol%basis%bas)) deallocate(active_mol%basis%bas)
        if (allocated(active_mol%basis%env)) deallocate(active_mol%basis%env)
        if (allocated(active_mol%basis%ao_loc)) deallocate(active_mol%basis%ao_loc)

        allocate(atm_raw(6 * int(active_mol%basis%natm)))
        allocate(bas_raw(8 * int(active_mol%basis%nbas)))
        allocate(active_mol%basis%env(env_size))

        open(newunit=unit_id, file='ints/atm.txt', status='old')
        read(unit_id, *) atm_raw
        close(unit_id)

        open(newunit=unit_id, file='ints/bas.txt', status='old')
        read(unit_id, *) bas_raw
        close(unit_id)

        open(newunit=unit_id, file='ints/env.txt', status='old')
        read(unit_id, *) active_mol%basis%env
        close(unit_id)

        allocate(active_mol%basis%atm(6, int(active_mol%basis%natm)))
        allocate(active_mol%basis%bas(8, int(active_mol%basis%nbas)))
        active_mol%basis%atm = reshape(atm_raw, shape(active_mol%basis%atm))
        active_mol%basis%bas = reshape(bas_raw, shape(active_mol%basis%bas))
        deallocate(atm_raw, bas_raw)

        allocate(active_mol%basis%ao_loc(int(active_mol%basis%nbas) + 1))
        active_mol%basis%ao_loc(1) = 1_c_int
        do i = 1, int(active_mol%basis%nbas)
            l = int(active_mol%basis%bas(2, i))
            nctr = int(active_mol%basis%bas(4, i))
            nao_shl = (2*l + 1) * nctr
            active_mol%basis%ao_loc(i+1) = active_mol%basis%ao_loc(i) + int(nao_shl, c_int)
        end do

        call sync_compatibility_views()
        if (present(mol)) mol = active_mol

        print *, "Molecule initialized:"
        print *, "  Atoms:       ", active_mol%basis%natm
        print *, "  Shells:      ", active_mol%basis%nbas
        print *, "  AO functions:", active_mol%basis%nao
        print *, "  Electrons:   ", active_mol%total_electrons
    end subroutine init_molecule

    subroutine activate_molecule(mol)
        implicit none
        type(molecule), intent(in) :: mol

        active_mol = mol
        call sync_compatibility_views()
    end subroutine activate_molecule

    subroutine build_hcore_overlap(mol, S, Hcore)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double), allocatable, intent(out) :: S(:,:), Hcore(:,:)

        integer(c_int) :: shls(2), cdims(2), i, j, p
        integer :: dim_i, dim_j, max_shell_dim
        integer :: ai, aj, buf_idx
        integer :: nao_local, nbas_local
        integer :: i0, j0
        real(c_double), allocatable :: buf(:)

        call activate_molecule(mol)
        nao_local = int(mol%basis%nao)
        nbas_local = int(mol%basis%nbas)

        allocate(S(nao_local, nao_local))
        allocate(Hcore(nao_local, nao_local))
        S     = 0.0d0
        Hcore = 0.0d0

        max_shell_dim = maxval(ao_loc(2:nbas_local+1) - ao_loc(1:nbas_local))
        allocate(buf(max_shell_dim * max_shell_dim))
        do i = 1, nbas_local
            dim_i = ao_loc(i+1) - ao_loc(i)
            i0 = ao_loc(i)
            do j = 1, nbas_local
                dim_j = ao_loc(j+1) - ao_loc(j)
                j0 = ao_loc(j)
                cdims = [int(dim_i, c_int), int(dim_j, c_int)]

                shls(1) = i - 1
                shls(2) = j - 1

                buf(1:dim_i * dim_j) = 0.0d0
                p = cint1e_ovlp_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                do aj = 1, dim_j
                    do ai = 1, dim_i
                        buf_idx = (aj - 1) * dim_i + ai
                        S(i0 + ai - 1, j0 + aj - 1) = buf(buf_idx)
                    end do
                end do

                buf(1:dim_i * dim_j) = 0.0d0
                p = cint1e_kin_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                do aj = 1, dim_j
                    do ai = 1, dim_i
                        buf_idx = (aj - 1) * dim_i + ai
                        Hcore(i0 + ai - 1, j0 + aj - 1) = Hcore(i0 + ai - 1, j0 + aj - 1) + buf(buf_idx)
                    end do
                end do

                buf(1:dim_i * dim_j) = 0.0d0
                p = cint1e_nuc_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                do aj = 1, dim_j
                    do ai = 1, dim_i
                        buf_idx = (aj - 1) * dim_i + ai
                        Hcore(i0 + ai - 1, j0 + aj - 1) = Hcore(i0 + ai - 1, j0 + aj - 1) + buf(buf_idx)
                    end do
                end do

            end do
        end do
        deallocate(buf)

        print *, "1-Electron integrals (S, Hcore) built."
    end subroutine build_hcore_overlap

    subroutine build_s_inv(S, S_inv_sqrt)
        implicit none
        real(c_double), intent(in) :: S(:,:)
        real(c_double), allocatable, intent(out) :: S_inv_sqrt(:,:)

        real(c_double), allocatable :: evalS(:), evecS(:,:), work(:)
        real(c_double), allocatable :: tmp(:,:)
        integer(c_int) :: lwork, info
        real(c_double) :: tmp_work(1)
        integer :: n, i, j

        n = size(S, 1)

        allocate(S_inv_sqrt(n,n))
        allocate(evalS(n))
        allocate(evecS(n,n))
        evecS = S

        lwork = -1
        call dsyev('V', 'U', int(n, c_int), evecS, int(n, c_int), evalS, tmp_work, int(lwork, c_int), info)
        lwork = int(tmp_work(1))
        allocate(work(lwork))

        call dsyev('V', 'U', int(n, c_int), evecS, int(n, c_int), evalS, work, int(lwork, c_int), info)
        if (info /= 0) then
            print *, "Error in dsyev for S, info = ", info
            stop
        end if
        deallocate(work)

        allocate(tmp(n,n))
        do i = 1, n
            do j = 1, n
                if (evalS(i) > 1.0d-12) then
                    tmp(i,j) = (1.0d0 / sqrt(evalS(i))) * evecS(j,i)
                else
                    tmp(i,j) = 0.0d0
                end if
            end do
        end do

        call dgemm('N', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, evecS, int(n, c_int), tmp, int(n, c_int), &
                   0.0d0, S_inv_sqrt, int(n, c_int))

        deallocate(evalS)
        deallocate(evecS)
        deallocate(tmp)

        print *, "S^{-1/2} constructed."
    end subroutine build_s_inv

    subroutine build_initial_guess(mol, Hcore, X, P)
        implicit none
        type(molecule), intent(in) :: mol
        real(c_double), intent(in) :: Hcore(:,:), X(:,:)
        real(c_double), allocatable, intent(out) :: P(:,:)

        real(c_double), allocatable :: F_orth(:,:), C_orth(:,:), C(:,:)
        real(c_double), allocatable :: evalF(:), evecF(:,:), work(:)
        real(c_double) :: tmp_work(1)
        integer(c_int) :: lwork, info
        integer :: n, nocc

        n = size(Hcore, 1)

        allocate(P(n,n))
        allocate(F_orth(n,n))
        allocate(C(n,n))
        allocate(C_orth(n,n))
        call dgemm('N', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, Hcore, int(n, c_int), X, int(n, c_int), &
                   0.0d0, C_orth, int(n, c_int))
        call dgemm('T', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, X, int(n, c_int), C_orth, int(n, c_int), &
                   0.0d0, F_orth, int(n, c_int))

        allocate(evalF(n))
        allocate(evecF(n,n))
        evecF = F_orth

        lwork = -1
        call dsyev('V', 'U', int(n, c_int), evecF, int(n, c_int), evalF, tmp_work, int(lwork, c_int), info)
        lwork = int(tmp_work(1))
        allocate(work(lwork))

        call dsyev('V', 'U', int(n, c_int), evecF, int(n, c_int), evalF, work, int(lwork, c_int), info)
        if (info /= 0) then
            print *, "Error in dsyev for F_orth, info = ", info
            stop
        end if
        deallocate(work)

        call dgemm('N', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, X, int(n, c_int), evecF, int(n, c_int), &
                   0.0d0, C, int(n, c_int))

        nocc = mol%total_electrons / 2
        P = 0.0d0

        call dgemm('N', 'T', int(n, c_int), int(n, c_int), int(nocc, c_int), &
                   2.0d0, C(:,1:nocc), int(n, c_int), C(:,1:nocc), int(n, c_int), &
                   0.0d0, P, int(n, c_int))
        ! do j = 1, n
        !     do i = 1, n
        !         do occ = 1, nocc
        !             P(i,j) = P(i,j) + 2.0d0 * C(i,occ) * C(j,occ)
        !         end do
        !     end do
        ! end do

        deallocate(F_orth)
        deallocate(C_orth)
        deallocate(C)
        deallocate(evalF)
        deallocate(evecF)

        print *, "Initial guess Density Matrix (P) constructed."
    end subroutine build_initial_guess

    subroutine ensure_molecule_loaded()
        implicit none

        if (.not. allocated(active_mol%basis%env)) then
            call init_molecule(active_mol)
        else if (.not. allocated(atm) .or. .not. allocated(bas) .or. .not. allocated(ao_loc)) then
            call sync_compatibility_views()
        end if
    end subroutine ensure_molecule_loaded

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

end module one_eints
