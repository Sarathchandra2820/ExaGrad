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

    type :: fragment_scf_info
        type(molecule) :: mol
        real(c_double), allocatable :: S(:,:), C_mo(:,:)
        real(c_double), allocatable :: mo_energies(:)
        integer(c_int) :: nocc = 0
        integer(c_int) :: nvir = 0
    end type fragment_scf_info
        



end module molecule_t