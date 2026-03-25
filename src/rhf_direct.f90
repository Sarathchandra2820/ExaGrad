program rhf_direct
    use iso_c_binding
    use libcint_interface
    implicit none

    integer(c_int) :: u, i, j, shls(2), cdims(1)
    integer :: natm, nbas, env_size, nao
    
    integer(c_int), allocatable :: atm(:), bas(:)
    real(c_double), allocatable :: env(:)
    
    integer, allocatable :: ao_loc(:)
    integer :: l, nctr, nao_shl, dim_i, dim_j
    integer :: ao_i, ao_j, p, q
    
    real(c_double), allocatable :: S(:,:), Hcore(:,:), buf(:,:)
    
    print *, "=========================================="
    print *, " Fortran Direct RHF SCF powered by libcint"
    print *, "=========================================="

    ! ---------------------------------------------------------
    ! 1. Read the sizes from PySCF export
    ! ---------------------------------------------------------
    open(newunit=u, file='ints/mol_info.txt', status='old')
    read(u, *) natm
    read(u, *) nbas
    read(u, *) env_size
    read(u, *) nao
    close(u)
    
    ! ---------------------------------------------------------
    ! 2. Allocate and load the flat C-arrays
    ! PySCF exported flat arrays exactly formatted for CINT.
    ! ---------------------------------------------------------
    allocate(atm(natm * 6))     ! Each atm block has 6 elements
    allocate(bas(nbas * 8))     ! Each bas block has 8 elements
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
    
    print *, "Successfully loaded molecule and basis arrays!"
    print *, "  Atoms: ", natm
    print *, "  Shells: ", nbas
    print *, "  AO functions: ", nao

    ! ---------------------------------------------------------
    ! 3. Map Basis Shells to Atomic Orbital (AO) Indices
    ! ---------------------------------------------------------
    allocate(ao_loc(nbas + 1))
    ao_loc(1) = 1  ! Fortran 1-based indexing for arrays
    
    do i = 1, nbas
        ! In bas array: index 1 is L (angular momentum)
        ! index 3 is nctr (number of contracted functions)
        ! Fortran is 1-based, but bas is conceptually blocks of 8.
        l    = bas((i-1)*8 + 2) 
        nctr = bas((i-1)*8 + 4)
        
        nao_shl = (2*l + 1) * nctr
        ao_loc(i+1) = ao_loc(i) + nao_shl
    end do
    
    ! ---------------------------------------------------------
    ! 4. Evaluate 1-Electron Integrals (Overlap & Core Hamiltonian)
    ! ---------------------------------------------------------
    allocate(S(nao, nao))
    allocate(Hcore(nao, nao))
    S = 0.0d0
    Hcore = 0.0d0
    
    print *, "Generating 1-Electron Integrals (S and Hcore)..."
    
    ! libcint works on Shell Pairs (i, j)
    do i = 1, nbas
        dim_i = ao_loc(i+1) - ao_loc(i)
        
        do j = 1, nbas
            dim_j = ao_loc(j+1) - ao_loc(j)
            
            allocate(buf(dim_i, dim_j))
            buf = 0.0d0
            
            ! shls requires 0-based indices for CINT
            shls(1) = i - 1
            shls(2) = j - 1
            
            cdims(1) = 0
            ! 4.1 Compute Overlap (S)
            p = cint1e_ovlp_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
            do ao_i = 1, dim_i
                do ao_j = 1, dim_j
                    S(ao_loc(i)+ao_i-1, ao_loc(j)+ao_j-1) = buf(ao_i, ao_j)
                end do
            end do
            
            ! 4.2 Compute Kinetic Energy (T)
            buf = 0.0d0
            p = cint1e_kin_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
            do ao_i = 1, dim_i
                do ao_j = 1, dim_j
                    Hcore(ao_loc(i)+ao_i-1, ao_loc(j)+ao_j-1) = Hcore(ao_loc(i)+ao_i-1, ao_loc(j)+ao_j-1) + buf(ao_i, ao_j)
                end do
            end do
            
            ! 4.3 Compute Nuclear Attraction (V)
            buf = 0.0d0
            p = cint1e_nuc_sph(buf, cdims, shls, atm, natm, bas, nbas, env, c_null_ptr, c_null_ptr)
            do ao_i = 1, dim_i
                do ao_j = 1, dim_j
                    Hcore(ao_loc(i)+ao_i-1, ao_loc(j)+ao_j-1) = Hcore(ao_loc(i)+ao_i-1, ao_loc(j)+ao_j-1) + buf(ao_i, ao_j)
                end do
            end do
            
            deallocate(buf)
        end do
    end do
    
    print *, "1-Electron Integral phase complete!"
    print *, "Hcore(1,1) = ", Hcore(1,1)
    
end program rhf_direct
