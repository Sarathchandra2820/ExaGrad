module rose_global

! Module to hold global information like file units.

  implicit none

! Integer-type  file unit data
  integer, public :: luinp         = 5 ! Input file unit
  integer, public :: luout         = 6 ! Output file unit
 
  save

end module
