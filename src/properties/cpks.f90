module cpks
    use iso_c_binding, only: c_double, c_int, c_null_ptr
    use libcint_interface, only: cint1e_r_sph
    use one_eints, only: activate_molecule, atm, bas, env
    use molecule_t, only: molecule
    use math_utils
    implicit none

contains


    subroutine build_trial_vector(C_mo, C_mo_occ, trial_vector)
        implicit none
        real(c_double), intent(in) :: C_mo(:,:), C_mo_occ(:)
        real(c_double), intent(out) :: trial_vector(:)
        trial_vector = 0.0d0
        if (size(trial_vector) > 0 .and. size(C_mo, 1) > 0 .and. size(C_mo_occ) > 0) then
            trial_vector(1) = trial_vector(1) + 0.0d0 * (C_mo(1,1) + C_mo_occ(1))
        end if
        print *, "CPKS computation placeholder"
    end subroutine build_trial_vector

end module cpks
