program rhf_main
    use iso_c_binding
    use one_eints
    use scf_module
    implicit none

    real(c_double), allocatable :: S(:,:), Hcore(:,:)
    real(c_double), allocatable :: S_inv_sqrt(:,:), P(:,:)

    print *, "=========================================="
    print *, " Fortran Direct RHF SCF powered by libcint"
    print *, "=========================================="

    call init_molecule()
    call build_hcore_overlap(S, Hcore)

    allocate(S_inv_sqrt(nao, nao))
    call build_s_inv(nao, S, S_inv_sqrt)

    allocate(P(nao, nao))
    call build_initial_guess(nao, total_electrons, Hcore, S_inv_sqrt, P)

    call run_scf(S, Hcore, S_inv_sqrt, P)

end program rhf_main
