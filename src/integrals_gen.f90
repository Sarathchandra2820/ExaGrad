module one_eints
    use, intrinsic :: iso_c_binding
    use libcint_interface
    use math_utils
    implicit none

    ! Publicly accessible molecule and basis data
    integer :: natm, nbas, env_size, nao
    integer :: total_electrons
    integer(c_int), allocatable :: atm(:), bas(:)
    real(c_double), allocatable :: env(:)
    integer, allocatable :: ao_loc(:)

contains

    ! ------------------------------------------------------------------
    ! init_molecule: Read atm/bas/env from PySCF export and build ao_loc
    ! ------------------------------------------------------------------
    subroutine init_molecule()
        implicit none
        integer :: u, i, l, nctr, nao_shl

        open(newunit=u, file='ints/mol_info.txt', status='old')
        read(u, *) natm
        read(u, *) nbas
        read(u, *) env_size
        read(u, *) nao
        read(u, *) total_electrons
        close(u)

        allocate(atm(natm * 6))
        allocate(bas(nbas * 8))
        allocate(env(env_size))

        open(newunit=u, file='ints/atm.txt', status='old')
        do i = 1, size(atm)
            read(u, *) atm(i)
        end do
        close(u)

        open(newunit=u, file='ints/bas.txt', status='old')
        do i = 1, size(bas)
            read(u, *) bas(i)
        end do
        close(u)

        open(newunit=u, file='ints/env.txt', status='old')
        do i = 1, size(env)
            read(u, *) env(i)
        end do
        close(u)

        ! Build the shell -> AO index mapping
        allocate(ao_loc(nbas + 1))
        ao_loc(1) = 1
        do i = 1, nbas
            l    = bas((i-1)*8 + 2)   ! Angular momentum of shell i
            nctr = bas((i-1)*8 + 4)   ! Number of contracted functions
            nao_shl = (2*l + 1) * nctr
            ao_loc(i+1) = ao_loc(i) + nao_shl
        end do

        print *, "Molecule initialized:"
        print *, "  Atoms:       ", natm
        print *, "  Shells:      ", nbas
        print *, "  AO functions:", nao
        print *, "  Electrons:   ", total_electrons
    end subroutine init_molecule

    ! ------------------------------------------------------------------
    ! build_hcore_overlap: Evaluate S and Hcore via libcint 1e integrals
    ! ------------------------------------------------------------------
    subroutine build_hcore_overlap(S, Hcore)
        implicit none
        real(c_double), allocatable, intent(out) :: S(:,:), Hcore(:,:)

        integer(c_int) :: shls(2), cdims(1), i, j, ao_i, ao_j, p
        integer :: dim_i, dim_j
        real(c_double), allocatable :: buf(:,:)

        allocate(S(nao, nao))
        allocate(Hcore(nao, nao))
        S     = 0.0d0
        Hcore = 0.0d0
        cdims(1) = 0

        do i = 1, nbas
            dim_i = ao_loc(i+1) - ao_loc(i)
            do j = 1, nbas
                dim_j = ao_loc(j+1) - ao_loc(j)
                allocate(buf(dim_i, dim_j))

                shls(1) = i - 1   ! libcint uses 0-based shell indices
                shls(2) = j - 1

                ! Overlap S
                buf = 0.0d0
                p = cint1e_ovlp_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                S(ao_loc(i):ao_loc(i+1)-1, ao_loc(j):ao_loc(j+1)-1) = buf

                ! Kinetic Energy T
                buf = 0.0d0
                p = cint1e_kin_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                Hcore(ao_loc(i):ao_loc(i+1)-1, ao_loc(j):ao_loc(j+1)-1) = &
                    Hcore(ao_loc(i):ao_loc(i+1)-1, ao_loc(j):ao_loc(j+1)-1) + buf

                ! Nuclear Attraction V
                buf = 0.0d0
                p = cint1e_nuc_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
                Hcore(ao_loc(i):ao_loc(i+1)-1, ao_loc(j):ao_loc(j+1)-1) = &
                    Hcore(ao_loc(i):ao_loc(i+1)-1, ao_loc(j):ao_loc(j+1)-1) + buf

                deallocate(buf)
            end do
        end do

        print *, "1-Electron integrals (S, Hcore) built."
    end subroutine build_hcore_overlap

    ! ------------------------------------------------------------------
    ! build_s_inv: Computes S^{-1/2} from Overlap Matrix S
    ! ------------------------------------------------------------------
    subroutine build_s_inv(n, S, S_inv_sqrt)
        implicit none
        integer, intent(in) :: n
        real(c_double), intent(in) :: S(n,n)
        real(c_double), intent(out) :: S_inv_sqrt(n,n)

        real(c_double), allocatable :: evalS(:), evecS(:,:), work(:)
        real(c_double), allocatable :: tmp(:,:)
        integer(c_int) :: lwork, info, i, j
        real(c_double) :: tmp_work(1)

        allocate(evalS(n))
        allocate(evecS(n,n))
        evecS = S

        ! Query optimal workspace for dsyev
        lwork = -1
        call dsyev('V', 'U', int(n, c_int), evecS, int(n, c_int), evalS, tmp_work, int(lwork, c_int), info)
        lwork = int(tmp_work(1))
        allocate(work(lwork))

        ! Diagonalize S
        call dsyev('V', 'U', int(n, c_int), evecS, int(n, c_int), evalS, work, int(lwork, c_int), info)
        if (info /= 0) then
            print *, "Error in dsyev for S, info = ", info
            stop
        end if
        deallocate(work)

        ! Form s^{-1/2} * U^T = tmp
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

        ! Form S^{-1/2} = U * tmp
        ! S_inv_sqrt = evecS * tmp
        call dgemm('N', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, evecS, int(n, c_int), tmp, int(n, c_int), &
                   0.0d0, S_inv_sqrt, int(n, c_int))

        deallocate(evalS)
        deallocate(evecS)
        deallocate(tmp)
        
        print *, "S^{-1/2} constructed."
    end subroutine build_s_inv

    ! ------------------------------------------------------------------
    ! build_initial_guess: Construct initial density matrix P from Hcore
    ! ------------------------------------------------------------------
    subroutine build_initial_guess(n, nelec, Hcore, X, P)
        implicit none
        integer, intent(in) :: n, nelec
        real(c_double), intent(in) :: Hcore(n,n), X(n,n)
        real(c_double), intent(out) :: P(n,n)

        real(c_double), allocatable :: F_orth(:,:), C_orth(:,:), C(:,:)
        real(c_double), allocatable :: evalF(:), evecF(:,:), work(:)
        real(c_double) :: tmp_work(1)
        integer(c_int) :: lwork, info, i, j, nocc

        allocate(F_orth(n,n))
        allocate(C(n,n))

        ! F_orth = X^T * Hcore * X
        ! First do tmp = Hcore * X
        allocate(C_orth(n,n)) ! use C_orth as tmp space first
        call dgemm('N', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, Hcore, int(n, c_int), X, int(n, c_int), &
                   0.0d0, C_orth, int(n, c_int))
        ! Now F_orth = X^T * tmp
        call dgemm('T', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, X, int(n, c_int), C_orth, int(n, c_int), &
                   0.0d0, F_orth, int(n, c_int))

        ! Diagonalize F_orth
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

        ! The eigenvectors of F_orth are C_orth (in evecF)
        ! MO Coefficients C = X * C_orth
        call dgemm('N', 'N', int(n, c_int), int(n, c_int), int(n, c_int), &
                   1.0d0, X, int(n, c_int), evecF, int(n, c_int), &
                   0.0d0, C, int(n, c_int))

        ! Construct density matrix P
        ! P_ij = 2.0 * sum_{k=1}^{nocc} C_ik * C_jk
        nocc = nelec / 2
        P = 0.0d0
        do j = 1, n
            do i = 1, n
                do lwork = 1, nocc
                    P(i,j) = P(i,j) + 2.0d0 * C(i,lwork) * C(j,lwork)
                end do
            end do
        end do

        deallocate(F_orth)
        deallocate(C_orth)
        deallocate(C)
        deallocate(evalF)
        deallocate(evecF)
        
        print *, "Initial guess Density Matrix (P) constructed. E_core = ", sum(P * Hcore) / 2.0d0
        
    end subroutine build_initial_guess

end module one_eints
