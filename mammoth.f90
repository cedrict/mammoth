!==============================================!
!                                              !
! C. Thieulot ; May 2018                       !
!                                              !
!==============================================!
                                               !
program mammoth                                !
                                               !
use structures                                 !
                                               !
implicit none                                  !
                                               !
integer, parameter :: m=8                      ! number of nodes which constitute an element
integer, parameter :: ndofT=1                  ! number of dofs per node
integer, parameter :: nrhs=1
integer nnx                                    ! number of grid points in the x direction
integer nny                                    ! number of grid points in the y direction
integer nnz                                    ! number of grid points in the z direction
integer np                                     ! number of grid points
integer nelx                                   ! number of elements in the x direction
integer nely                                   ! number of elements in the y direction
integer nelz                                   ! number of elements in the z direction
integer nel                                    ! number of elements
integer NfemT                                  ! size of the FEM matrix \
integer nchamber                               ! number of magma chambers
integer itopo                                  ! type of surface topography
integer istep,nstep                            !
integer ipar(16)                               ! params for CG solver
integer inoode(m)                              !
integer, dimension(:,:), allocatable :: icon   ! connectivity array
integer, dimension(:), allocatable :: ipvt     ! work array needed by the solver 
integer, dimension(:), allocatable :: counting !
real(8), dimension(:), allocatable :: wk       ! cg solver array
real(8), dimension(:), allocatable :: sol      ! cg solver array
integer, dimension(:), allocatable :: iw       ! cg solver array
                                               !
integer i1,i2,i,j,k,iel,counter,iq,jq,kq       !
integer ik,jk,ikk,jkk,m1,m2,k1,k2,info         !
integer nz,j1,j2,ip,jp,nsees,maxits            !
integer output_frequency,inode,jnode           !
integer size_krylov_subspace,imc,idummy        !
integer init_temp,numarg,option_ID             !
                                               !  
real(8), parameter :: year=3600d0*24d0*365d0   !
real(8), parameter :: pi=4.d0*atan(1.d0)       ! 
real(8) Lx,Ly,Lz                               ! size of the numerical domain
real(8) density                                ! mass density $\rho$ of the material
real(8), dimension(:),  allocatable :: x,y,z   ! node coordinates arrays
real(8), dimension(:),  allocatable :: T,Tss   ! temperature field 
real(8), dimension(:),  allocatable :: rho     !  
real(8), dimension(:),  allocatable :: B       ! right hand side
real(8), dimension(:,:),allocatable :: A       ! FEM matrix
real(8), dimension(:),  allocatable :: qx,qy,qz!
real(8), dimension(:),  allocatable :: qznode  !
                                               !
real(8) rq,sq,tq,weightq                       ! local coordinate and weight of qpoint
real(8) xq,yq,zq                               ! global coordinate of qpoint
real(8) Ael(m*ndofT,m*ndofT)                   ! elemental FEM matrix
real(8) Bel(m*ndofT)                           ! elemental right hand side
real(8) N(m),dNdx(m),dNdy(m),dNdz(m)           !
real(8) dNdr(m),dNds(m),dNdt(m)                ! shape fcts and derivatives
real(8) M3D(m,m),F3D(m)                        !
real(8) jcob                                   ! determinant of jacobian matrix
real(8) jcb(3,3)                               ! jacobian matrix
real(8) jcbi(3,3)                              ! inverse of jacobian matrix
real(8) Aref,eps                               !
real(8) Bmat(3,8)                              ! B matrix
real(8) BmatT(8,3)                             ! B matrix transpose
real(8) Nvect(1,8),NvectT(8,1)                 !
real(8) temp(8)                                ! 
real(8) alpha,dt,dtwish                        !
real(8) heat_capacity                          !
real(8) heat_conductivity                      !
real(8) heat_production                        !
real(8) heatflux_bottom,heatflux_top           !
real(8) Kc(m,m)                                !
real(8) Tsurf,Tbot                             !
real(8) time,CFL_number                        !
real(8) surf,fixt,xi,zmax,dist                 !
real(8) dTdx,dTdy,dTdz,minqz,maxqz             !
real(8) rtol,atol,sstolerance                  !
real(8) sx,sy,sz,t1,t2                         !
real(8) fpar(16)                               ! params for CG solver
                                               !
logical,dimension(:),allocatable :: bc_fixT    ! prescribed b.c. array
logical,dimension(:),allocatable :: bottom_elt !
logical use_heatflux_bottom                    !
logical use_heatflux_top                       !
logical use_lapack                             ! 
logical debug                                  !
                                               !
character(len=4) cistep                        !
character(len=255) arg                         !
character(len=16) inputfile                    !
                                               !
external cg                                    ! Conjugate Gradient Solver
                                               !
!==============================================!

write(*,'(a)') '//+++++++++++++++++++++++++++++++++++++++++++++++//'
write(*,'(a)') '//+++++++++++++++++ MAMMOTH +++++++++++++++++++++//'
write(*,'(a)') '//+++++++++++++++++++++++++++++++++++++++++++++++//'

numarg = command_argument_count()

if (numarg/=1) stop 'You forgot the input file. Duh.'
option_ID=1
call getarg(option_ID,arg)
read(arg,*) inputfile

print *,'input file: ',trim(inputfile)

!==============================================!

open(unit=111,file=trim(inputfile),action='read')
open(unit=1000,file="OUT/heat_flux_statistics.dat")
open(unit=1001,file="OUT/heat_flux_statistics_qz_surface.dat")
open(unit=1002,file="OUT/temperature_statistics.dat")
open(unit=1003,file="OUT/xi.dat")
open(unit=1004,file="OUT/dt.dat")

!==============================================!


read(111,*) Lx,Ly,Lz
read(111,*) nelx,nely,nelz

nnx=nelx+1
nny=nely+1
nnz=nelz+1
np=nnx*nny*nnz
nel=nelx*nely*nelz
NfemT=np*ndofT
sx=Lx/nelx
sy=Ly/nely
sz=Lz/nelz

print *,'Lx=',Lx,'m'
print *,'Ly=',Ly,'m'
print *,'Lz=',Lz,'m'
print *,'nelx=',nelx
print *,'nely=',nely
print *,'nelz=',nelz
print *,'nnx=',nnx
print *,'nny=',nny
print *,'nnz=',nnz
print *,'np=',np
print *,'NfemT=',NfemT
print *,'sx=',sx,'m'
print *,'sy=',sy,'m'
print *,'sz=',sz,'m'

!==============================================!

read(111,*) density
read(111,*) heat_capacity
read(111,*) heat_conductivity
read(111,*) heat_production

print *,'density=',density
print *,'heat_capacity',heat_capacity
print *,'heat_conductivity',heat_conductivity
print *,'heat_production',heat_production

read(111,*) idummy
use_heatflux_bottom=(idummy==1)

if (use_heatflux_bottom) then 
read(111,*) heatflux_bottom
print *,'heatflux_bottom=',heatflux_bottom
else
read(111,*) Tbot
print *,'Tbot=',Tbot
end if

use_heatflux_top=.false.

read(111,*) Tsurf
print *,'Tsurf=',Tsurf

read(111,*) init_temp
print *,'init_temp=',init_temp

read(111,*) dtwish
print *,'dtwish=',dtwish,'yr'
dtwish=dtwish*year

!==============================================!
! numerical parameters 
!==============================================!

read(111,*) nstep
print *,'nstep=',nstep

read(111,*) output_frequency
print *,'output_frequency=',output_frequency

read(111,*) sstolerance
print *,'sstolerance=',sstolerance

read(111,*) itopo
print *,'itopo=',itopo

!==============================================!
! magma chambers layout/properties
!==============================================!

read(111,*) nchamber
print *, 'nchamber=',nchamber

allocate(mc(nchamber))

read(111,*) mc(1:nchamber)%ishape
read(111,*) mc(1:nchamber)%xc
read(111,*) mc(1:nchamber)%yc
read(111,*) mc(1:nchamber)%zc
read(111,*) mc(1:nchamber)%a
read(111,*) mc(1:nchamber)%b
read(111,*) mc(1:nchamber)%c
read(111,*) mc(1:nchamber)%T

mc(1:nchamber)%xc=mc(1:nchamber)%xc*Lx
mc(1:nchamber)%yc=mc(1:nchamber)%yc*Ly
mc(1:nchamber)%zc=mc(1:nchamber)%zc*Lz

print *,'chambers ishape:',mc(:)%ishape
print *,'chambers xc:',mc(:)%xc
print *,'chambers yc:',mc(:)%yc
print *,'chambers zc:',mc(:)%zc

!==============================================!
! numerical stuff / do not change
!==============================================!

eps=1.d-10
alpha=0.5
use_lapack=.false.
debug=.false.
CFL_number=0.1 ! not more than 0.5

print *,'CFL_number=',CFL_number

!==============================================!
!===[allocate memory]==========================!
!==============================================!

allocate(x(np))
allocate(y(np))
allocate(z(np))
allocate(T(np))
allocate(Tss(np)) 
allocate(icon(m,nel))
allocate(bc_fixT(NfemT))
allocate(rho(nel))
allocate(qx(nel))
allocate(qy(nel))
allocate(qz(nel))
allocate(qznode(np))
allocate(counting(np))
allocate(bottom_elt(nel)) ; bottom_elt=.false.
if (use_lapack) then
allocate(A(NfemT,NfemT))
allocate(B(NfemT))
end if

!==============================================!
!===[grid points setup]========================!
!==============================================!

counter=0
do i=0,nelx
do j=0,nely
do k=0,nelz
   counter=counter+1
   x(counter)=dble(i)*Lx/dble(nelx)
   y(counter)=dble(j)*Ly/dble(nely)
   z(counter)=dble(k)*Lz/dble(nelz)
end do
end do
end do

!==============================================!
!===[connectivity]=============================!
!==============================================!

counter=0
do i=1,nelx
do j=1,nely
do k=1,nelz
counter=counter+1
icon(1,counter)=nny*nnz*(i-1)+nnz*(j-1)+k
icon(2,counter)=nny*nnz*(i  )+nnz*(j-1)+k
icon(3,counter)=nny*nnz*(i  )+nnz*(j  )+k
icon(4,counter)=nny*nnz*(i-1)+nnz*(j  )+k
icon(5,counter)=nny*nnz*(i-1)+nnz*(j-1)+k+1
icon(6,counter)=nny*nnz*(i  )+nnz*(j-1)+k+1
icon(7,counter)=nny*nnz*(i  )+nnz*(j  )+k+1
icon(8,counter)=nny*nnz*(i-1)+nnz*(j  )+k+1
if (k==1) bottom_elt(counter)=.true.
end do
end do
end do

!==============================================!
!=====[matrix setup in CSR format]=============!
!==============================================!

csrA%nr=np
csrA%nc=np

csrA%nz=8*8                                  & ! 8 corners with 8 neighbours
       +(nnx-2)*(nny-2)*(nnz-2)*27           & ! all the inside nodes with 27 neighbours
       +(4*(nnx-2)+4*(nny-2)+4*(nnz-2))*12   & ! the edge nodes with 12 neighbours  
       +2*(nnx-2)*(nny-2)*18                 & ! 2 faces
       +2*(nnx-2)*(nnz-2)*18                 & ! 2 faces
       +2*(nny-2)*(nnz-2)*18                   ! 2 faces
      
csrA%nz=(csrA%nz-csrA%nr)/2+csrA%nr

!write(*,'(a)')       'CSR matrix format                   ||' 
print *,'csrA%nr=',csrA%nr
print *,'csrA%nz=',csrA%nz

allocate(csrA%ia(csrA%nr+1))
allocate(csrA%ja(csrA%nz))  
allocate(csrA%mat(csrA%nz)) 
allocate(csrA%rhs(csrA%nr)) 

   counter=0
   nz=0
   csrA%ia(1)=1
   do i1=1,nnx
   do j1=1,nny
   do k1=1,nnz
      ip=nny*nnz*(i1-1)+(j1-1)*nnz + k1 ! node number
      nsees=0
      do i2=-1,1 ! exploring neighbouring nodes
      do j2=-1,1 ! exploring neighbouring nodes
      do k2=-1,1 ! exploring neighbouring nodes
         i=i1+i2
         j=j1+j2
         k=k1+k2
         if (i>=1 .and. i<= nnx .and. j>=1 .and. j<=nny .and. k>=1 .and. k<=nnz) then ! if node exists
            jp=nny*nnz*(i-1)+(j-1)*nnz + k ! node number
if (jp>=ip) then
            nz=nz+1
            csrA%ja(nz)=jp
            nsees=nsees+1
            counter=counter+1
end if
         end if
      end do
      end do
      end do
      csrA%ia(ip+1)=csrA%ia(ip)+nsees
   end do
   end do
   end do

!print *,counter,csrA%nz

!==============================================!
!=====[compute dt]=============================!
!==============================================!

dt=min(CFL_number*(min(sx,sy,sz))**2/(heat_conductivity/density/heat_capacity),dtwish)

write(1004,*) istep,dt/year,dtwish/year,time/year ; call flush(1004)

print *,'dt=',dt/year,'yr'

!==============================================!
!=====[initial temperature]====================!
!==============================================!

! generate steady state temperature field

do i=1,np
   T(i)=-heatflux_bottom/heat_conductivity*(z(i)-Lz)+Tsurf
end do

Tss=T

! add to it magma chambers

do imc=1,nchamber
   select case(mc(imc)%ishape)
   case(1) ! spheroid  
      do i=1,np
      if ( (x(i)-mc(imc)%xc)**2/mc(imc)%a**2 + &
           (y(i)-mc(imc)%yc)**2/mc(imc)%b**2 + &
           (z(i)-mc(imc)%zc)**2/mc(imc)%c**2 <= 1.d0) T(i)=mc(imc)%T
      end do
   case(2) ! cuboid  
      do i=1,np
      if ( abs(x(i)-mc(imc)%xc)<mc(imc)%a/2 .and. &
           abs(y(i)-mc(imc)%yc)<mc(imc)%b/2 .and. &
           abs(z(i)-mc(imc)%zc)<mc(imc)%c/2 ) T(i)=mc(imc)%T
      end do
   case(3) ! pancake  
      do i=1,np
      if ( (x(i)-mc(imc)%xc)**2/mc(imc)%a**2 + &
           (y(i)-mc(imc)%yc)**2/mc(imc)%b**2 <= 1.d0 .and. &
           abs(z(i)-mc(imc)%zc)<mc(imc)%c/2 ) T(i)=mc(imc)%T
      end do
   case default
      stop 'wrong chamber type'
   end select
end do

!==============================================!
!=====[implement topgraphy]====================!
!==============================================!

! deform free surface

counter=0
do i=1,nnx
do j=1,nny
do k=1,nnz
   counter=counter+1
   if (k==nnz) then
      select case(itopo)
      case(1) 
         if (abs(x(counter)-Lx/2)<20d3) &
         z(counter)=z(counter)+3.d3*(1+cos((x(counter)-Lx/2.)/20d3*pi))
      case(2)
         if (abs(y(counter)-Ly/2)<20d3) &
         z(counter)=z(counter)+1.d3*(1+cos((y(counter)-Lx/2.)/20d3*pi))
      case(3)
         dist=sqrt((x(counter)-Lx/2.)**2+(y(counter)-Ly/2.)**2)
         if (dist<20.d3) &
         z(counter)=z(counter)+2.5d3*(1+cos(dist/20d3*pi))

      case default
         stop 'itopo parameter unknown'
      end select
      write(999,*) x(counter),y(counter),z(counter)-Lz
   end if 
end do
end do
end do

! adapt mesh underneath

counter=0
do i=1,nnx
do j=1,nny
do k=1,nnz
   counter=counter+1
   
   jp=nny*nnz*(i-1)+nnz*(j-1)+nnz
   zmax=z(jp)

   if (k<nnz .and. k>1) then
   z(counter)=zmax/nelz*(k-1)
   end if 
end do
end do
end do






!**************************************************************************************************
!****** time stepping loop ************************************************************************
!**************************************************************************************************

time=0

do istep=1,nstep 

write(cistep,'(I4.4)') istep

write(*,'(a)') '//+++++++++++++++++++++++++++++++++++++++++++++++//'
write(*,'(a)') '//+++++++++++++ istep='//cistep//'+++++++++++++++++++++++//'
write(*,'(a)') '//+++++++++++++++++++++++++++++++++++++++++++++++//'

!==============================================!
!=====[define bc]==============================!
!==============================================!

bc_fixT=.false.

ip=0
do i=1,nnx
do j=1,nny
do k=1,nnz
   ip=ip+1 
   if (.not.use_heatflux_bottom) then
      if (k==1) then
         bc_fixT(ip)=.true. ; T(ip)=Tbot
      endif
   end if
   if (.not.use_heatflux_top) then
      if (k==nnz) then
         bc_fixT(ip)=.true. ; T(ip)=Tsurf
      endif
   endif
end do
end do
end do

!==============================================!
!=====[build FE matrix]========================!
!==============================================!

call cpu_time(t1)

if (use_lapack) then 
   A=0.d0
   B=0.d0
else
   csrA%mat=0
   csrA%rhs=0
end if

do iel=1,nel

   Ael=0.d0
   Bel=0.d0

   inoode(1:m)=icon(1:m,iel)

   temp(1:m)=T(inoode)

   do iq=-1,1,2
   do jq=-1,1,2
   do kq=-1,1,2

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      tq=kq/sqrt(3.d0)
      weightq=1.d0*1.d0*1.d0

      N(1)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0-tq)
      N(2)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0-tq)
      N(3)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0-tq)
      N(4)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0-tq)
      N(5)=0.125d0*(1.d0-rq)*(1.d0-sq)*(1.d0+tq)
      N(6)=0.125d0*(1.d0+rq)*(1.d0-sq)*(1.d0+tq)
      N(7)=0.125d0*(1.d0+rq)*(1.d0+sq)*(1.d0+tq)
      N(8)=0.125d0*(1.d0-rq)*(1.d0+sq)*(1.d0+tq)

      dNdr(1)= - 0.125d0*(1.d0-sq)*(1.d0-tq)    
      dNdr(2)= + 0.125d0*(1.d0-sq)*(1.d0-tq)    
      dNdr(3)= + 0.125d0*(1.d0+sq)*(1.d0-tq)    
      dNdr(4)= - 0.125d0*(1.d0+sq)*(1.d0-tq)    
      dNdr(5)= - 0.125d0*(1.d0-sq)*(1.d0+tq)    
      dNdr(6)= + 0.125d0*(1.d0-sq)*(1.d0+tq)    
      dNdr(7)= + 0.125d0*(1.d0+sq)*(1.d0+tq)    
      dNdr(8)= - 0.125d0*(1.d0+sq)*(1.d0+tq)    

      dNds(1)= - 0.125d0*(1.d0-rq)*(1.d0-tq)    
      dNds(2)= - 0.125d0*(1.d0+rq)*(1.d0-tq)    
      dNds(3)= + 0.125d0*(1.d0+rq)*(1.d0-tq)  
      dNds(4)= + 0.125d0*(1.d0-rq)*(1.d0-tq)
      dNds(5)= - 0.125d0*(1.d0-rq)*(1.d0+tq)    
      dNds(6)= - 0.125d0*(1.d0+rq)*(1.d0+tq)    
      dNds(7)= + 0.125d0*(1.d0+rq)*(1.d0+tq)    
      dNds(8)= + 0.125d0*(1.d0-rq)*(1.d0+tq)    

      dNdt(1)= - 0.125d0*(1.d0-rq)*(1.d0-sq)    
      dNdt(2)= - 0.125d0*(1.d0+rq)*(1.d0-sq)    
      dNdt(3)= - 0.125d0*(1.d0+rq)*(1.d0+sq)    
      dNdt(4)= - 0.125d0*(1.d0-rq)*(1.d0+sq)    
      dNdt(5)= + 0.125d0*(1.d0-rq)*(1.d0-sq)    
      dNdt(6)= + 0.125d0*(1.d0+rq)*(1.d0-sq)    
      dNdt(7)= + 0.125d0*(1.d0+rq)*(1.d0+sq)  
      dNdt(8)= + 0.125d0*(1.d0-rq)*(1.d0+sq) 

      jcb=0.d0    
      do k=1,8    
      jcb(1,1)=jcb(1,1)+dNdr(k)*x(inoode(k))
      jcb(1,2)=jcb(1,2)+dNdr(k)*y(inoode(k))
      jcb(1,3)=jcb(1,3)+dNdr(k)*z(inoode(k))
      jcb(2,1)=jcb(2,1)+dNds(k)*x(inoode(k))
      jcb(2,2)=jcb(2,2)+dNds(k)*y(inoode(k))
      jcb(2,3)=jcb(2,3)+dNds(k)*z(inoode(k))
      jcb(3,1)=jcb(3,1)+dNdt(k)*x(inoode(k))
      jcb(3,2)=jcb(3,2)+dNdt(k)*y(inoode(k))
      jcb(3,3)=jcb(3,3)+dNdt(k)*z(inoode(k))
      enddo    
    
      jcob=jcb(1,1)*jcb(2,2)*jcb(3,3) &    
          +jcb(1,2)*jcb(2,3)*jcb(3,1) &    
          +jcb(2,1)*jcb(3,2)*jcb(1,3) &    
          -jcb(1,3)*jcb(2,2)*jcb(3,1) &    
          -jcb(1,2)*jcb(2,1)*jcb(3,3) &    
          -jcb(2,3)*jcb(3,2)*jcb(1,1)    

      jcbi(1,1)=(jcb(2,2)*jcb(3,3)-jcb(2,3)*jcb(3,2))/jcob    
      jcbi(2,1)=(jcb(2,3)*jcb(3,1)-jcb(2,1)*jcb(3,3))/jcob    
      jcbi(3,1)=(jcb(2,1)*jcb(3,2)-jcb(2,2)*jcb(3,1))/jcob  
      jcbi(1,2)=(jcb(1,3)*jcb(3,2)-jcb(1,2)*jcb(3,3))/jcob
      jcbi(2,2)=(jcb(1,1)*jcb(3,3)-jcb(1,3)*jcb(3,1))/jcob    
      jcbi(3,2)=(jcb(1,2)*jcb(3,1)-jcb(1,1)*jcb(3,2))/jcob    
      jcbi(1,3)=(jcb(1,2)*jcb(2,3)-jcb(1,3)*jcb(2,2))/jcob    
      jcbi(2,3)=(jcb(1,3)*jcb(2,1)-jcb(1,1)*jcb(2,3))/jcob    
      jcbi(3,3)=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))/jcob    

      xq=0.d0
      yq=0.d0
      zq=0.d0
      do k=1,m
         xq=xq+N(k)*x(icon(k,iel))
         yq=yq+N(k)*y(icon(k,iel))
         zq=zq+N(k)*z(icon(k,iel))
         dNdx(k)=jcbi(1,1)*dNdr(k)&    
                +jcbi(1,2)*dNds(k)&    
                +jcbi(1,3)*dNdt(k)    
         dNdy(k)=jcbi(2,1)*dNdr(k)&    
                +jcbi(2,2)*dNds(k)&    
                +jcbi(2,3)*dNdt(k)    
         dNdz(k)=jcbi(3,1)*dNdr(k)&    
                +jcbi(3,2)*dNds(k)&    
                +jcbi(3,3)*dNdt(k)  
      end do

      Nvect(1,:)=N(:)  
      NvectT(:,1)=N(:) 

      Bmat(1,:)=dNdx 
      Bmat(2,:)=dNdy
      Bmat(3,:)=dNdz

      BmatT=transpose(Bmat)     

      Kc=matmul(BmatT,Bmat)*heat_conductivity*weightq*jcob 

      M3D=matmul(NvectT,Nvect)*density*heat_capacity*weightq*jcob 

      F3D=N(:)*jcob*weightq*heat_production

      Ael=Ael+(M3D+Kc*alpha*dt)   

      Bel=Bel+matmul(M3D-Kc*(1.d0-alpha)*dt,temp(1:8))+F3D*dt    

   end do
   end do
   end do

   !==============================
   !=====[impose bc conditions]===
   !==============================

   do i=1,m
      inode=icon(i,iel)
      if (bc_fixT(inode)) then
      fixt=T(inode)
      Aref=Ael(i,i)
      do j=1,m
         Bel(j)=Bel(j)-Ael(j,i)*fixt
         Ael(i,j)=0.d0
         Ael(j,i)=0.d0
      enddo
      Ael(i,i)=Aref
      Bel(i)=Aref*fixt
      endif
   enddo

   !=====================
   !=====[impose flux]===
   !=====================

   if (use_heatflux_bottom .and. bottom_elt(iel)) then
      inoode(1)=icon(1,iel)  
      inoode(2)=icon(2,iel)  
      inoode(3)=icon(3,iel)  
      inoode(4)=icon(4,iel)  
      surf=(x(inoode(3))-x(inoode(1)))*(y(inoode(3))-y(inoode(1)))
      if (.not. bc_fixT(inoode(1))) Bel(1)=Bel(1)+surf*heatflux_bottom*0.25d0*dt
      if (.not. bc_fixT(inoode(2))) Bel(2)=Bel(2)+surf*heatflux_bottom*0.25d0*dt
      if (.not. bc_fixT(inoode(3))) Bel(3)=Bel(3)+surf*heatflux_bottom*0.25d0*dt
      if (.not. bc_fixT(inoode(4))) Bel(4)=Bel(4)+surf*heatflux_bottom*0.25d0*dt
   end if

   !=====================
   !=====[assemble]======
   !=====================

   if (use_lapack) then
      do k1=1,m
         ik=icon(k1,iel)
         do i1=1,ndofT
            ikk=ndofT*(k1-1)+i1
            m1=ndofT*(ik-1)+i1
            do k2=1,m
               jk=icon(k2,iel)
               do i2=1,ndofT
                  jkk=ndofT*(k2-1)+i2
                  m2=ndofT*(jk-1)+i2
                  A(m1,m2)=A(m1,m2)+Ael(ikk,jkk)
               end do
            end do
            B(m1)=B(m1)+Bel(ikk)
         end do
      end do
   else
      do k1=1,m
      inode=icon(k1,iel)
      do k2=1,m
         jnode=icon(k2,iel)
         do k=csrA%ia(inode),csrA%ia(inode+1)-1
            if (csrA%ja(k)==jnode) then
               csrA%mat(k)=csrA%mat(k)+Ael(k1,k2)
               exit
            end if
         end do
      end do
      csrA%rhs(inode)=csrA%rhs(inode)+Bel(k1)
      end do
   end if

end do

call cpu_time(t2) ; write(*,'(a,f8.3,a)') 'make matrix: ',t2-t1,'s'

!==============================================!
!=====[solve system]===========================!
!==============================================!

call cpu_time(t1)

if (use_lapack) then
   if (debug) then
   print *,'A :',minval(A),maxval(A)
   print *,'B :',minval(B),maxval(B)
   end if
   allocate(ipvt(NfemT))
   call DGESV( NfemT, NRHS, A, NfemT, ipvt, B, NfemT, INFO )
   deallocate(ipvt)
   T=B
else
   if (debug) then
   print *,'A :',minval(csrA%mat),maxval(csrA%mat)
   print *,'B :',minval(csrA%rhs),maxval(csrA%rhs)
   end if
   maxits=1000
   size_krylov_subspace=100
   rtol=1.d-12!1.d-8
   atol=1.d50 ! I do not want to use atol
   allocate(wk(csrA%nr*10))
   allocate(iw(csrA%nr*3))
   allocate(sol(csrA%nr))
   ipar(2)=0 ! no prec
   ipar(3)=1 ! stopping criterion
   ipar(4)=csrA%nr*10
   ipar(5)=100 ! size krylov subspace
   ipar(6)=1000 ! maxits
   fpar(1)=1.d-28 ! rtol
   fpar(2)=1.d50 ! atol
   sol=0
   call runrc(csrA%nr,csrA%nz,csrA%rhs,sol, &
              ipar,fpar,wk,T,csrA%mat,csrA%ja,csrA%ia,cg)
   T=sol
   deallocate(wk)
   deallocate(iw)
   deallocate(sol)
end if

call cpu_time(t2) ; write(*,'(a,f8.3,a)') 'solve system: ',t2-t1,'s'

write(*,'(a,2f10.3,a)') 'T (m/M):',minval(T),maxval(T),' C'

time=time+dt

write(*,'(a,f12.2,a)') 'time=',time/year,' yr'

write(1002,*) time/year,minval(T),maxval(T),sum(T)/NfemT,istep ; call flush(1002)

!==============================================!
!=====[compute heat flux]======================!
!==============================================!

do iel=1,nel

   inoode(1:m)=icon(1:m,iel)

   rq=0
   sq=0
   tq=0

   dNdr(1)= - 0.125d0*(1.d0-sq)*(1.d0-tq)    
   dNdr(2)= + 0.125d0*(1.d0-sq)*(1.d0-tq)    
   dNdr(3)= + 0.125d0*(1.d0+sq)*(1.d0-tq)    
   dNdr(4)= - 0.125d0*(1.d0+sq)*(1.d0-tq)    
   dNdr(5)= - 0.125d0*(1.d0-sq)*(1.d0+tq)    
   dNdr(6)= + 0.125d0*(1.d0-sq)*(1.d0+tq)    
   dNdr(7)= + 0.125d0*(1.d0+sq)*(1.d0+tq)    
   dNdr(8)= - 0.125d0*(1.d0+sq)*(1.d0+tq)    

   dNds(1)= - 0.125d0*(1.d0-rq)*(1.d0-tq)    
   dNds(2)= - 0.125d0*(1.d0+rq)*(1.d0-tq)    
   dNds(3)= + 0.125d0*(1.d0+rq)*(1.d0-tq)  
   dNds(4)= + 0.125d0*(1.d0-rq)*(1.d0-tq)
   dNds(5)= - 0.125d0*(1.d0-rq)*(1.d0+tq)    
   dNds(6)= - 0.125d0*(1.d0+rq)*(1.d0+tq)    
   dNds(7)= + 0.125d0*(1.d0+rq)*(1.d0+tq)    
   dNds(8)= + 0.125d0*(1.d0-rq)*(1.d0+tq)    

   dNdt(1)= - 0.125d0*(1.d0-rq)*(1.d0-sq)    
   dNdt(2)= - 0.125d0*(1.d0+rq)*(1.d0-sq)    
   dNdt(3)= - 0.125d0*(1.d0+rq)*(1.d0+sq)    
   dNdt(4)= - 0.125d0*(1.d0-rq)*(1.d0+sq)    
   dNdt(5)= + 0.125d0*(1.d0-rq)*(1.d0-sq)    
   dNdt(6)= + 0.125d0*(1.d0+rq)*(1.d0-sq)    
   dNdt(7)= + 0.125d0*(1.d0+rq)*(1.d0+sq)  
   dNdt(8)= + 0.125d0*(1.d0-rq)*(1.d0+sq) 

   jcb=0.d0    
   do k=1,8    
      jcb(1,1)=jcb(1,1)+dNdr(k)*x(inoode(k))
      jcb(1,2)=jcb(1,2)+dNdr(k)*y(inoode(k))
      jcb(1,3)=jcb(1,3)+dNdr(k)*z(inoode(k))
      jcb(2,1)=jcb(2,1)+dNds(k)*x(inoode(k))
      jcb(2,2)=jcb(2,2)+dNds(k)*y(inoode(k))
      jcb(2,3)=jcb(2,3)+dNds(k)*z(inoode(k))
      jcb(3,1)=jcb(3,1)+dNdt(k)*x(inoode(k))
      jcb(3,2)=jcb(3,2)+dNdt(k)*y(inoode(k))
      jcb(3,3)=jcb(3,3)+dNdt(k)*z(inoode(k))
   enddo    
    
   jcob=jcb(1,1)*jcb(2,2)*jcb(3,3) &    
       +jcb(1,2)*jcb(2,3)*jcb(3,1) &    
       +jcb(2,1)*jcb(3,2)*jcb(1,3) &    
       -jcb(1,3)*jcb(2,2)*jcb(3,1) &    
       -jcb(1,2)*jcb(2,1)*jcb(3,3) &    
       -jcb(2,3)*jcb(3,2)*jcb(1,1)    

   jcbi(1,1)=(jcb(2,2)*jcb(3,3)-jcb(2,3)*jcb(3,2))/jcob    
   jcbi(2,1)=(jcb(2,3)*jcb(3,1)-jcb(2,1)*jcb(3,3))/jcob    
   jcbi(3,1)=(jcb(2,1)*jcb(3,2)-jcb(2,2)*jcb(3,1))/jcob  
   jcbi(1,2)=(jcb(1,3)*jcb(3,2)-jcb(1,2)*jcb(3,3))/jcob
   jcbi(2,2)=(jcb(1,1)*jcb(3,3)-jcb(1,3)*jcb(3,1))/jcob    
   jcbi(3,2)=(jcb(1,2)*jcb(3,1)-jcb(1,1)*jcb(3,2))/jcob    
   jcbi(1,3)=(jcb(1,2)*jcb(2,3)-jcb(1,3)*jcb(2,2))/jcob    
   jcbi(2,3)=(jcb(1,3)*jcb(2,1)-jcb(1,1)*jcb(2,3))/jcob    
   jcbi(3,3)=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))/jcob    

   do k=1,m
      dNdx(k)=jcbi(1,1)*dNdr(k)&    
             +jcbi(1,2)*dNds(k)&    
             +jcbi(1,3)*dNdt(k)    
      dNdy(k)=jcbi(2,1)*dNdr(k)&    
             +jcbi(2,2)*dNds(k)&    
             +jcbi(2,3)*dNdt(k)    
      dNdz(k)=jcbi(3,1)*dNdr(k)&    
             +jcbi(3,2)*dNds(k)&    
             +jcbi(3,3)*dNdt(k)  
   end do

   dTdx = sum(dNdx(1:m)*T(inoode(1:m)))
   dTdy = sum(dNdy(1:m)*T(inoode(1:m)))
   dTdz = sum(dNdz(1:m)*T(inoode(1:m)))

   qx(iel)=heat_conductivity*dTdx
   qy(iel)=heat_conductivity*dTdy
   qz(iel)=heat_conductivity*dTdz

end do

qznode=0
counting=0
do iel=1,nel
   qznode(icon(1,iel))=qznode(icon(1,iel))+qz(iel)
   qznode(icon(2,iel))=qznode(icon(2,iel))+qz(iel)
   qznode(icon(3,iel))=qznode(icon(3,iel))+qz(iel)
   qznode(icon(4,iel))=qznode(icon(4,iel))+qz(iel)
   qznode(icon(5,iel))=qznode(icon(5,iel))+qz(iel)
   qznode(icon(6,iel))=qznode(icon(6,iel))+qz(iel)
   qznode(icon(7,iel))=qznode(icon(7,iel))+qz(iel)
   qznode(icon(8,iel))=qznode(icon(8,iel))+qz(iel)
   counting(icon(1,iel))=counting(icon(1,iel))+1
   counting(icon(2,iel))=counting(icon(2,iel))+1
   counting(icon(3,iel))=counting(icon(3,iel))+1
   counting(icon(4,iel))=counting(icon(4,iel))+1
   counting(icon(5,iel))=counting(icon(5,iel))+1
   counting(icon(6,iel))=counting(icon(6,iel))+1
   counting(icon(7,iel))=counting(icon(7,iel))+1
   counting(icon(8,iel))=counting(icon(8,iel))+1
end do
qznode=qznode/counting

write(1000,*) time/year,minval(qx),maxval(qx),&
                        minval(qy),maxval(qy),&
                        minval(qz),maxval(qz),&
                        minval(qznode),maxval(qznode) ; call flush(1000)

minqz=+1d50
maxqz=-1d50
counter=0
do i=1,nelx
do j=1,nely
do k=1,nelz
   counter=counter+1
   if (k==nelz) then
      minqz=min(minqz,abs(qz(counter)))
      maxqz=max(maxqz,abs(qz(counter)))
   end if
end do
end do
end do
write(1001,*) time/year,minqz,maxqz ; call flush(1001)

write(*,'(a,2f10.5)') 'qx (m/M)',minval(qx),maxval(qx)
write(*,'(a,2f10.5)') 'qy (m/M)',minval(qy),maxval(qy)
write(*,'(a,2f10.5)') 'qz (m/M)',minval(qz),maxval(qz)

!==============================================!

call cpu_time(t1)

!if (istep==1 .or. mod(istep,output_frequency)==0) then
!   open(unit=123,file='OUT/solution_'//cistep//'.dat',status='replace')
!   do i=1,np
!   write(123,'(4f20.10)') x(i),y(i),z(i),T(i)
!   end do
!   close(123)
!end if

if (istep==1 .or. mod(istep,output_frequency)==0) &
call output_for_paraview (istep,np,nel,x,y,z,T,Tss,icon,rho,qx,qy,qz,qznode)

!==============================================!
! output T profiles at the surface

if (istep==1 .or. mod(istep,output_frequency)==0) then
open(unit=2000,file='OUT/temperature_z_profile_middle_'//cistep//'.dat',status='replace')
open(unit=2001,file='OUT/heat_flux_x_profile_'//cistep//'.dat',status='replace')
open(unit=2002,file='OUT/heat_flux_y_profile_'//cistep//'.dat',status='replace')
open(unit=2003,file='OUT/temperature_z_profile_all_'//cistep//'.dat',status='replace')
do ip=1,np
   if (abs(x(ip)-Lx/2.)<eps .and. abs(y(ip)-Ly/2.)<eps) then
      write(2000,*) z(ip),T(ip)
   end if
   if (abs(y(ip)-Ly/2.)<eps .and. abs(z(ip)-Lz)<eps) then
      write(2001,*) x(ip),abs(qznode(ip))
   end if
   if (abs(x(ip)-Lx/2.)<eps .and. abs(z(ip)-Lz)<eps) then
      write(2002,*) y(ip),abs(qznode(ip))
   end if
   write(2003,*) z(ip),T(ip)
end do
close(2000)
close(2001)
close(2002)
close(2003)
end if

call cpu_time(t2) ; write(*,'(a,f8.3,a)') 'output: ',t2-t1,'s'

!==============================================!
! assess whether model is back to steady state temperature field

xi=maxval(abs(T-Tss))

write(1003,*) istep,xi,sstolerance ; call flush(1003)

write(*,'(a,f7.2,a,f7.2)') 'xi (C/yr)=',xi,' sstolerance=',sstolerance

if (xi<sstolerance) stop 'steady state reached'


end do ! time stepping loop

!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************

end program





