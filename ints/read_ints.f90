module integral_read
    implicit none

    interface read_ints
        module procedure read_tensor, read_matrix, read_dipole
    end interface

contains

    subroutine read_tensor(filename, tensor, dim)
        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in) :: dim
        real, allocatable, intent(out) :: tensor(:,:,:,:)

        !---local variables ----

        integer :: u, ierr, i, j, k, l

        open(newunit=u, file=filename, status='old'. action='read', isostat=ierr)
        if(ierr /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if


        allocate(tensor(dim, dim, dim, dim))
        
        do i=1,dim
            do j=1,dim
                do k=1,dim 
                    do l=1,dim
                        read(u, *) tensor(i,j,k,l)
                    end do
                end do
            end do
        end do

        close(u)

    end subroutine read_tensor

    subroutine read_matrix(filename, matrix, dim)
        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in) :: dim
        real, allocatable, intent(out) :: matrix(:,:)

        !---local variables ----

        integer :: u, ierr, i, j

        open(newunit=u, file=filename, status='old'. action='read', isostat=ierr)
        if(ierr /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if

        allocate(matrix(dim, dim))
        
        do i=1,dim
            do j=1,dim
                read(u, *) matrix(i,j)
            end do
        end do

        close(u)

    end subroutine read_matrix

    subroutine read_dipole(filename, d_x, d_y, d_z, dim)
        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in) :: dim
        real, allocatable, intent(out) :: d_x(:,:), d_y(:,:), d_z(:,:)

        !---local variables ----

        integer :: u, ierr, i, j

        open(newunit=u, file=filename, status='old'. action='read', isostat=ierr)
        if(ierr /= 0) then
            print *, 'Error opening file: ', filename
            stop
        end if
        allocate(d_x(dim, dim), d_y(dim, dim), d_z(dim, dim))
        
        do i=1,dim
            do j=1,dim
                read(u, *) d_x(i,j)
            end do
        end do

        do i=1,dim
            do j=1,dim
                read(u, *) d_y(i,j)
            end do
        end do

        do i=1,dim
            do j=1,dim
                read(u, *) d_z(i,j)
            end do
        end do
        close(u)

    end subroutine read_dipole
    
end module integral_read








    