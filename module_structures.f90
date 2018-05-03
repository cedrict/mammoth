module structures

type compressedrowstorage                                   
   integer nr                             ! number of rows of (full) matrix
   integer nc                             ! number of columns of (full) matrix
   integer nz                             ! number of nonzeros
   integer NZ_SYMM
   integer,dimension(:),allocatable :: ia 
   integer,dimension(:),allocatable :: ja
   real(8),dimension(:),allocatable :: mat
   real(8),dimension(:),allocatable :: rhs
end type compressedrowstorage

type(compressedrowstorage) csrA

type chamber
   real(8) a,b,c
   real(8) xc,yc,zc
   real(8) emplacement
   real(8) T 
   integer ishape 
end type

type(chamber), dimension(:), allocatable :: mc

end module
