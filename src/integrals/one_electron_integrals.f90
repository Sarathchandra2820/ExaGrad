module one_eints
    use, intrinsic :: iso_c_binding
    use libcint_interface
    use math_utils
    use molecule_t, only: molecule
    use molecule_loader, only: init_molecule, &
                               loader_activate_molecule => activate_molecule, &
                               loader_ensure_molecule_loaded => ensure_molecule_loaded, &
                               natm, nbas, nao, atm, bas, env, ao_loc
    implicit none

contains

    subroutine activate_molecule(mol)
        implicit none
        type(molecule), intent(in) :: mol

        call loader_activate_molecule(mol)
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

        call loader_ensure_molecule_loaded()
    end subroutine ensure_molecule_loaded

end module one_eints
