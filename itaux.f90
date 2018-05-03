subroutine runrc(n,nz,rhs,sol,ipar,fpar,wk,guess,a,ja,ia,solver)

implicit none

integer n,nz,ipar(16),ia(n+1),ja(nz)
real(8) fpar(16),rhs(n),sol(n),guess(n),wk(*),a(nz)
external solver
character(len=42) shift

!c-----------------------------------------------------------------------
!c     the actual tester. It starts the iterative linear system solvers
!c     with a initial guess suppied by the user.
!c
!c     The structure {au, jau, ju} is assumed to have the output from
!c     the ILU* routines in ilut.f.
!c
!c-----------------------------------------------------------------------
!c     local variables
!c
      integer i, iou, its
      real(8) res, dnrm2
!c     real dtime, dt(2), time
!c     external dtime
      external dnrm2
      save its,res


shift = '                                       ||'

!c
!c     ipar(2) can be 0, 1, 2, please don't use 3
!c
!      if (ipar(2).gt.2) then
!         print *, 'I can not do both left and right preconditioning.'
!         return
!      endif
!c
!c     normal execution
!c
      its = 0
      res = 0.0D0

      do i = 1, n
         sol(i) = guess(i)
      enddo

      iou = 700
      ipar(1) = 0
!c     time = dtime(dt)
 10   call solver(n,rhs,sol,ipar,fpar,wk)
!c
!c     output the residuals
!c
      if (ipar(7).ne.its) then
         write (iou, *) its, real(res) ; call flush(iou)
         its = ipar(7)
      endif
      res = fpar(5)
!c
      if (ipar(1).eq.1) then
         !call amux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         call spmv_symm(n,nz, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
         goto 10

         !if (its<25) goto 10

      !else if (ipar(1).eq.2) then
      !   call atmux(n, wk(ipar(8)), wk(ipar(9)), a, ja, ia)
      !   goto 10
      !else if (ipar(1).eq.3 .or. ipar(1).eq.5) then
      !   call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
      !   goto 10
      !else if (ipar(1).eq.4 .or. ipar(1).eq.6) then
      !   call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
      !   goto 10
      else if (ipar(1).le.0) then
         if (ipar(1).eq.0) then
            !write(*,'(a)') shift//'Iterative solver converged'
         else if (ipar(1).eq.-1) then
            print *, 'Iterative solver has iterated too many times.'
         else if (ipar(1).eq.-2) then
            print *, 'Iterative solver was not given enough work space.'
            print *, 'The work space should at least have ', ipar(4),' elements.'
         else if (ipar(1).eq.-3) then
            print *, 'Iterative sovler is facing a break-down.'
         else
            print *, 'Iterative solver terminated. code =', ipar(1)
         endif
      endif
!c     time = dtime(dt)

      write (6,'(a,i6,a,es9.3)') 'inner: #its:',ipar(7),' res:',res 

      write(702,*) ipar(7)

      write (701, *) '# retrun code =', ipar(1),'convergence rate =', fpar(7)
!c     write (iou, *) '# total execution time (sec)', time

!*****check the error***********************

      !call amux(n,sol,wk,a,ja,ia)
      call spmv_symm(n,nz,sol,wk,a,ja,ia)
      do i = 1, n
         wk(n+i) = sol(i) -1.0D0
         wk(i) = wk(i) - rhs(i)
      enddo
      write (701, *) '# the actual residual norm is', dnrm2(n,wk,1)
      write (701, *) '# the error norm is', dnrm2(n,wk(1+n),1)

call flush(700)
call flush(701)
call flush(702)

      return
      end

!*****************************************
!function distdot(n,x,ix,y,iy)
!integer n, ix, iy
!real*8 distdot, x(*), y(*), ddot
!external ddot
!distdot = ddot(n,x,ix,y,iy)
!return
!end
!*****************************************
function distdot(n,x,ix,y,iy)
integer n, ix, iy
real(8) distdot, x(*), y(*)
distdot = dot_product(x(1:n),y(1:n))
return
end
!*****************************************



