subroutine output_for_paraview (istep,np,nel,x,y,z,T,Tss,icon,eldens,qx,qy,qz,qznode)
implicit none
integer istep,np,nel
real(8), dimension(np)    :: x,y,z,T,Tss,qznode
real(8), dimension(nel)   :: eldens,qx,qy,qz
integer, dimension(8,nel) :: icon

integer i,iel
character(len=4) cistep

!=======================================

write(cistep,'(I4.4)') istep

open(unit=123,file='OUT/visu_'//cistep//'.vtu',status='replace',form='formatted')
write(123,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(123,*) '<UnstructuredGrid>'
write(123,*) '<Piece NumberOfPoints="',np,'" NumberOfCells="',nel,'">'
!.............................
write(123,*) '<PointData Scalars="scalars">'

write(123,*) '<DataArray type="Float32" Name="temperature" Format="ascii">'
do i=1,np
write(123,*) T(i)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32" Name="temperature-ssT" Format="ascii">'
do i=1,np
write(123,*) T(i)-Tss(i)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32" Name="heat flux qz" Format="ascii">'
do i=1,np
write(123,*) qznode(i)
end do
write(123,*) '</DataArray>'




write(123,*) '</PointData>'
write(123,*) '<Points>'
write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
do i=1,np
write(123,*) x(i),y(i),z(i)
end do
write(123,*) '</DataArray>'
write(123,*) '</Points>'
!.............................
write(123,*) '<CellData Scalars="scalars">'
write(123,*) '<DataArray type="Float32" Name="density" Format="ascii">'
do iel=1,nel
write(123,*) eldens(iel)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Name="heat flux" Format="ascii">'
do iel=1,nel
write(123,*) qx(iel),qy(iel),qz(iel)
end do
write(123,*) '</DataArray>'
write(123,*) '</CellData>'
!.............................
write(123,*) '<Cells>'
write(123,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
do iel=1,nel
write(123,*) icon(1,iel)-1,icon(2,iel)-1,icon(3,iel)-1,icon(4,iel)-1,&
             icon(5,iel)-1,icon(6,iel)-1,icon(7,iel)-1,icon(8,iel)-1
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
write(123,*) (iel*8,iel=1,nel)
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Int32" Name="types" Format="ascii">'
write(123,*) (12,iel=1,nel)
write(123,*) '</DataArray>'
write(123,*) '</Cells>'
write(123,*) '</Piece>'
write(123,*) '</UnstructuredGrid>'
write(123,*) '</VTKFile>'
close(123)

end subroutine
