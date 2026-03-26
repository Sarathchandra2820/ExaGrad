module libcint_interface
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        integer(c_int) function cint1e_ovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(c, name="int1e_ovlp_sph")
            import :: c_int, c_double, c_ptr
            real(c_double), intent(out) :: out(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            type(c_ptr), value :: cache
        end function

        integer(c_int) function cint1e_nuc_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(c, name="int1e_nuc_sph")
            import :: c_int, c_double, c_ptr
            real(c_double), intent(out) :: out(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            type(c_ptr), value :: cache
        end function

        integer(c_int) function cint1e_kin_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(c, name="int1e_kin_sph")
            import :: c_int, c_double, c_ptr
            real(c_double), intent(out) :: out(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            type(c_ptr), value :: cache
        end function

        integer(c_int) function cint1e_r_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(c, name="int1e_r_sph")
            import :: c_int, c_double, c_ptr
            real(c_double), intent(out) :: out(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            type(c_ptr), value :: cache
        end function

        integer(c_int) function cint2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(c, name="int2e_sph")
            import :: c_int, c_double, c_ptr
            real(c_double), intent(out) :: out(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            type(c_ptr), value :: cache
        end function

        integer(c_int) function cint2c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(c, name="int2c2e_sph")
            import :: c_int, c_double, c_ptr
            real(c_double), intent(out) :: out(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            type(c_ptr), value :: cache
        end function

        integer(c_int) function cint3c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) bind(c, name="int3c2e_sph")
            import :: c_int, c_double, c_ptr
            real(c_double), intent(out) :: out(*)
            integer(c_int), intent(in) :: dims(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value :: opt
            type(c_ptr), value :: cache
        end function
    end interface

end module libcint_interface
