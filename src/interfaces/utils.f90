module utils
    use iso_c_binding, only: c_double
    implicit none
    

    contains 

    subroutine write_matrix_to_file(filename, data, status_msg)
        character(len=*), intent(in)  :: filename, status_msg
        real(c_double),   intent(in)  :: data(:,:)

        integer :: i, unit_id

        open(newunit=unit_id, file=filename, status=trim(status_msg), action='write', position='append')

        write(unit_id, '(A,I0,A,I0)') '# rows=', size(data, 1), '  cols=', size(data, 2)
        write(unit_id, '(A)') trim(status_msg)
        do i = 1, size(data, 1)
            write(unit_id, '(*(ES20.10))') data(i, :)
        end do
        close(unit_id)
    end subroutine write_matrix_to_file

end module utils

