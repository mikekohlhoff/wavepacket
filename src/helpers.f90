module helpers
  use settings, only : dp
  implicit none
  
  contains
  
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine expmat(z,n)
      implicit none
!     -----------------------------------------------------------------
!     Exponentiates a complex symmetric matrix z to give exp(z).
!     -----------------------------------------------------------------
      integer :: n, lwork, ierr, i, j, k
      complex(dp) :: z(n, n), zeta(n), v(n, n), work(2*n), s
      real(dp) :: rwork(2*n), zr, zi, exp, cos, sin
      
      lwork = 2*n
      call zgeev ('N','V',n,z,n,zeta,v,n,v,n,work,lwork,rwork,ierr)
      if (ierr .ne. 0) stop 'expmat 1'
      do j = 1,n
         s = (0.d0,0.d0)
         do i = 1,n
            s = s+v(i,j)*v(i,j)
         enddo
         s = 1.d0/sqrt(s)
         do i = 1,n
            v(i,j) = s*v(i,j)
         enddo
      enddo
      do j = 1,n
         do i = 1,j
            z(i,j) = (0.d0,0.d0)
         enddo
      enddo
      do k = 1,n
         zr = real(zeta(k), kind=dp)
         zi = aimag(zeta(k))
         zeta(k) = exp(zr)*cmplx(cos(zi),sin(zi), kind=dp)
         do j = 1,n
            s = v(j,k)*zeta(k)
            do i = 1,j
               z(i,j) = z(i,j)+v(i,k)*s
            enddo
         enddo
      enddo
      do j = 1,n
         do i = 1,j-1
            z(j,i) = z(i,j)
         enddo
      enddo
      return
    end subroutine expmat
    
    
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function pnorm (psi)
      use settings, only : nr, nt
      implicit none
!     -----------------------------------------------------------------
!     Wavepacket left on grid
!     -----------------------------------------------------------------
      complex(dp) :: psi(nr, nt)
      real(dp) :: dble
      real(dp) :: pnorm
      integer :: i, j
      pnorm = 0.d0
      do i = 1,nr
         do j=1,nt
            pnorm = pnorm+(dble(psi(i,j))**2+aimag(psi(i,j))**2)!*pi!*a(i)**2
         enddo
      enddo
      return
    end function pnorm

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function pnorm2 (psi,dist,r,cost,rmid)
      use settings, only : nr, nt
      implicit none
!     -----------------------------------------------------------------
!     Wavepacket left inside area bound by absorbing potential and surface
!     -----------------------------------------------------------------
      real(dp) :: rmid, dist
      complex(dp) :: psi(nr, nt), dble
      real(dp) :: r(nr),cost(nt)
      real(dp) :: d, pnorm2, z
      integer :: i, j
      d = -dist
      pnorm2 = 0.d0
      do i = 1,nr
         do j = 1,nt
            z = r(i)*cost(j)
            if ((z.gt.d).and.(r(i).lt.rmid)) then
!            if (r(i).lt.rmid) then
               pnorm2 = pnorm2+(aimag(psi(i,j))**2+dble(psi(i,j))**2)
            endif
         enddo
      enddo
      return
    end function pnorm2


end module helpers
