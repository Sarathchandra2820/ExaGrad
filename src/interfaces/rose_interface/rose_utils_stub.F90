! Minimal stub for the ROSE rose_utils module.
! Provides only the routines that are actually referenced at link time.
! make_RIBO_Fock is omitted to avoid pulling in rose_ao / rose_mo.
module rose_utils
    use rose_global
    implicit none
contains

    subroutine get_free_fileunit(unit)
        integer, intent(out) :: unit
        integer, parameter   :: LUN_MIN=100, LUN_MAX=1000
        logical :: opened
        integer :: lun
        unit = -1
        do lun = LUN_MIN, LUN_MAX
            inquire(unit=lun, opened=opened)
            if (.not. opened) then
                unit = lun
                return
            end if
        end do
    end subroutine get_free_fileunit

end module rose_utils
