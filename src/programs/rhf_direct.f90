program rhf_direct
    use iso_c_binding
    use molecule_t, only: molecule
    use one_eints, only: init_molecule, build_hcore_overlap
    implicit none

    type(molecule) :: mol
    real(c_double), allocatable :: S(:,:), Hcore(:,:)

    print *, "=========================================="
    print *, " Fortran Direct RHF SCF powered by libcint"
    print *, "=========================================="

    call init_molecule(mol)
    call build_hcore_overlap(mol, S, Hcore)

    print *, "1-Electron Integral phase complete!"
    print *, "Hcore(1,1) = ", Hcore(1,1)

end program rhf_direct
