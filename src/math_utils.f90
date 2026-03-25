module math_utils
    use iso_c_binding
    implicit none

    ! BLAS and LAPACK interfaces
    interface
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
    end interface

contains

end module math_utils
