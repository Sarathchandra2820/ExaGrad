module math_utils
    use iso_c_binding
    implicit none

    interface

        function ddot(n, x, incx, y, incy) bind(c, name="ddot_")
            import :: c_double, c_int
            integer(c_int), intent(in) :: n, incx, incy
            real(c_double), intent(in) :: x(*), y(*)
            real(c_double) :: ddot
        end function ddot
        subroutine daxpy(n, da, dx, incx, dy, incy) bind(c, name="daxpy_")
            import :: c_double, c_int
            integer(c_int), intent(in) :: n, incx, incy
            real(c_double), intent(in) :: da
            real(c_double), intent(in) :: dx(*)
            real(c_double), intent(inout) :: dy(*)
        end subroutine daxpy
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) bind(c, name="dgemm_")
            import :: c_double, c_int, c_char
            character(kind=c_char), intent(in) :: transa, transb
            integer(c_int), intent(in) :: m, n, k, lda, ldb, ldc
            real(c_double), intent(in) :: alpha, beta
            real(c_double), intent(in) :: a(*), b(*)
            real(c_double), intent(inout) :: c(*)
        end subroutine dgemm

        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info) bind(c, name="dsyev_")
            import :: c_double, c_int, c_char
            character(kind=c_char), intent(in) :: jobz, uplo
            integer(c_int), intent(in) :: n, lda
            real(c_double), intent(inout) :: a(*)
            real(c_double), intent(out) :: w(*)
            real(c_double), intent(inout) :: work(*)
            integer(c_int), intent(in) :: lwork
            integer(c_int), intent(out) :: info
        end subroutine dsyev

        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info) bind(c, name="dgesv_")
            import :: c_double, c_int
            integer(c_int), intent(in)  :: n, nrhs, lda, ldb
            real(c_double), intent(inout) :: a(*)
            integer(c_int), intent(out)   :: ipiv(*)
            real(c_double), intent(inout) :: b(*)
            integer(c_int), intent(out)   :: info
        end subroutine dgesv

        subroutine dpotrf(uplo, n, a, lda, info) bind(c, name="dpotrf_")
            import :: c_double, c_int, c_char
            character(kind=c_char), intent(in) :: uplo
            integer(c_int), intent(in) :: n, lda
            real(c_double), intent(inout) :: a(*)
            integer(c_int), intent(out) :: info
        end subroutine dpotrf

        subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) bind(c, name="dtrsm_")
            import :: c_double, c_int, c_char
            character(kind=c_char), intent(in) :: side, uplo, transa, diag
            integer(c_int), intent(in) :: m, n, lda, ldb
            real(c_double), intent(in) :: alpha
            real(c_double), intent(in) :: a(*)
            real(c_double), intent(inout) :: b(*)
        end subroutine dtrsm
    end interface

contains

end module math_utils
