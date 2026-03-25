program rhf_main
    use iso_c_binding
    use one_eints
    implicit none

    real(c_double), allocatable :: S(:,:), Hcore(:,:)

    print *, "=========================================="
    print *, " Fortran Direct RHF SCF powered by libcint"
    print *, "=========================================="

    ! Step 1: Read basis set and geometry from PySCF export
    call init_molecule()

    ! Step 2: Build the 1-electron integrals S and Hcore
    call build_hcore_overlap(S, Hcore)

    print *, "Hcore(1,1) = ", Hcore(1,1)
    print *, "S(1,1)     = ", S(1,1)

    ! Step 3: TODO - build S^{-1/2}, initial guess, Direct SCF loop
    


    
end program rhf_main
