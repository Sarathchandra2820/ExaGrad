module molecule_t
    use iso_c_binding
    implicit none

    type :: basis_context
        integer(c_int) :: natm = 0
        integer(c_int) :: nbas = 0
        integer(c_int) :: nao = 0
        integer(c_int), allocatable :: atm(:,:)
        integer(c_int), allocatable :: bas(:,:)
        real(c_double), allocatable :: env(:)
        integer(c_int), allocatable :: ao_loc(:)
    end type basis_context

    type :: molecule
        integer :: natoms = 0
        integer :: total_electrons = 0
        real(c_double), allocatable :: coords(:,:)
        character(len=10), allocatable :: symbols(:)
        type(basis_context) :: basis
    end type molecule

end module molecule_t