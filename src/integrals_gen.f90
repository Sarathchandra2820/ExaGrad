module one_eints
    use, intrinsic :: iso_c_binding
    use libcint_interface
    implicit none

    ! Publicly accessible molecule and basis data
    integer :: natm, nbas, env_size, nao
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

end module one_eints
