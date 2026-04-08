! Minimal stub for the ROSE properties module.
! The Localization subroutine in make_ibo.F90 declares "use properties" but
! does not actually call any routine from it at runtime.  We provide this
! stub so that the module is available at compile time without pulling in
! the full ROSE I/O dependency chain (read_xyz, rose_mo, check, …).
module properties
    use rose_global
    implicit none
end module properties
