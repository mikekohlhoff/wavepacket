module surface
  use settings, only : dp
  implicit none
  
  real(dp), parameter :: pi = dacos(-1.d0)
  
  ! chulkov potential parameters for metal surface of interest
  real(dp), parameter :: xhar = 27.2113845d0
  
  ! Jellium surface parameters for aluminium surface (see paper by Jennings and Jones)
  real(dp), parameter :: z0=  0.7d0
  real(dp), parameter :: v0 = -0.574d0
  real(dp), parameter :: beta = 1.25d0
  real(dp), parameter :: ajel = -1.d0-4.d0*v0/beta
  real(dp), parameter :: bjel = -2.d0*(v0/ajel)
  
  !-----Cu100---------------------      
  real(dp), parameter :: as = 3.415d0
  real(dp), parameter :: a10 = -11.480d0/xhar
  real(dp), parameter :: a1 = 6.1d0/xhar
  real(dp), parameter :: a2 = 3.782d0/xhar
  real(dp), parameter :: bs = 2.539d0
      
  !-----Cu111--------------
  !      as = 3.94d0
  !      a10 = -11.895d0/xhar
  !      a1 = 5.14d0/xhar
  !      a2 = 4.3279d0/xhar
  !      bs = 2.9416d0
  !------------------------
  real(dp), parameter :: z1 = 5.d0*pi/(4.d0*bs)
  real(dp), parameter :: a20 = a2-a10-a1
  real(dp), parameter :: a3 = -a20-a2/dsqrt(2.d0)
  real(dp), parameter :: aa = bs*a2*dsin(bs*z1)/a3
  real(dp), parameter :: xlam = 2.d0*aa
  real(dp), parameter :: zim = z1-1.d0/aa*dlog(-aa/(2.d0*a3))
  
  contains
  
    real(dp) function dpotsurf(r,z,d)
      implicit none
      real(dp) :: r, z, d
      real(dp) :: zz, vee, vep, xl
      zz = z+d
!     electron-image electron attraction
      if (zz.gt.z0) then
        vee = 0.25d0*(1.d0-dexp(-beta*(zz-z0)))/(zz-z0)**2 - 0.25d0*beta*dexp(-beta*(zz-z0))/(zz-z0)
      else
        vee = -v0/(ajel*dexp(bjel*(zz-z0))+1.d0)**2 * ajel*bjel*dexp(bjel*(zz-z0))
      endif
      if (zz.gt.0) then
        xl = dsqrt(r**2 + 4.d0*d**2 + 4.d0*d*z)
        vep = (-4.d0*d-2.d0*z)/xl**3
      else
        vep = 0.d0
      endif
      dpotsurf = vee+vep
      return
    end function dpotsurf

    double precision function potsurf(r,z,d)
      implicit none
      real(dp) :: r, z, d
      real(dp) :: zz, y2, vee, xl, dabs, vep
!     -------------------------------------------------------------------------
!     compute the surface potential at distance d from surface (surface is at z=-d)
!     -------------------------------------------------------------------------
      zz = z+d
      y2 = (r**2-z**2)
!     electron-image electron attraction for Chulkov pseudopotentials
!      pi = dacos(-1.d0)
!      if (zz.le.0.d0) then
!         vee = a10 + a1*dcos(2.d0*pi/as*(zz))
!      elseif ((zz.gt.0.d0).and.(zz.lt.z1)) then
!         vee = -a20 + a2*dcos(bs*(zz))
!      elseif ((zz.ge.z1).and.(zz.le.zim)) then
!         vee = a3*dexp(-aa*(zz-z1))
!      elseif (zz.gt.zim) then
!         vee = (dexp(-xlam*(zz-zim))-1.d0)/(4.d0*(zz-zim))
!      endif
!      
!     electron-image electron attraction
!     uncomment for Chulkov
      if (zz.gt.z0) then
       vee = 0.25d0*(-1.d0+dexp(-beta*(zz-z0)))/(zz-z0)
      else
       vee = v0/(ajel*dexp(bjel*(zz-z0))+1.d0)
      endif
!     electron-image proton repulsion
      xl = dsqrt(y2+(d+dabs(zz))**2)
!      if (zz.gt.0.d0) then
         vep = +1.d0/xl
!      else
!         vep = +1.d0/dsqrt(d**2+y2)
!      endif
      potsurf = vee+vep
      return
    end function potsurf
 

  
end module surface
