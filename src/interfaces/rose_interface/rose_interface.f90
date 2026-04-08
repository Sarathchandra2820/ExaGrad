module rose_interface_module
    use iso_c_binding
    implicit none

contains

    ! -----------------------------------------------------------------------
    logical function env_is_true(name)
        implicit none
        character(len=*), intent(in) :: name
        character(len=32) :: env_val
        integer :: env_len, env_stat

        env_is_true = .false.
        call get_environment_variable(name, env_val, length=env_len, status=env_stat)
        if (env_stat /= 0) return

        env_val = adjustl(env_val(1:env_len))
        if (trim(env_val) == '1' .or. trim(env_val) == 'true' .or. trim(env_val) == 'TRUE') then
            env_is_true = .true.
        end if
    end function env_is_true

    subroutine localize_mo_spaces(S, C_mo, nocc, C_loc_occ, C_loc_vir, localized_occ, localized_vir)
        implicit none
        real(c_double), intent(in)  :: S(:,:), C_mo(:,:)
        integer, intent(in)         :: nocc
        real(c_double), intent(out) :: C_loc_occ(size(C_mo,1), nocc)
        real(c_double), intent(out) :: C_loc_vir(size(C_mo,1), size(C_mo,2)-nocc)
        logical, intent(out)        :: localized_occ, localized_vir

        integer :: nmo, nvir
        logical :: localize_all

        nmo = size(C_mo, 2)
        nvir = nmo - nocc

        C_loc_occ = C_mo(:,1:nocc)
        if (nvir > 0) C_loc_vir = C_mo(:,nocc+1:)
        localized_occ = .false.
        localized_vir = .false.

        if (size(S,1) /= size(C_mo,1)) then
            print *, '  ROSE interface skipped: inconsistent AO dimensions'
            return
        end if

        localize_all = env_is_true('EXAGRAD_LOCALIZE_ALL')
        if (localize_all .or. env_is_true('EXAGRAD_LOCALIZE_OCC')) then
            print *, '  ROSE interface scaffold active (occupied localization placeholder)'
            localized_occ = .true.
        end if

        if (nvir > 0 .and. (localize_all .or. env_is_true('EXAGRAD_LOCALIZE_VIR'))) then
            print *, '  ROSE interface scaffold active (virtual localization placeholder)'
            localized_vir = .true.
        end if
    end subroutine localize_mo_spaces

    ! -----------------------------------------------------------------------
    subroutine localize_occupied_orbitals(S, C_mo, nocc, C_loc_occ, localized)
        implicit none
        real(c_double), intent(in)  :: S(:,:), C_mo(:,:)
        integer, intent(in)         :: nocc
        real(c_double), intent(out) :: C_loc_occ(size(C_mo,1), nocc)
        logical, intent(out)        :: localized
        real(c_double), allocatable :: C_loc_vir(:,:)
        logical :: localized_vir

        allocate(C_loc_vir(size(C_mo,1), size(C_mo,2)-nocc))
        call localize_mo_spaces(S, C_mo, nocc, C_loc_occ, C_loc_vir, localized, localized_vir)
        deallocate(C_loc_vir)
    end subroutine localize_occupied_orbitals

end module rose_interface_module
