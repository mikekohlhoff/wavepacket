      program cwdvr_wavepacket
c     calculates the time evolution of H atom Rydberg state
c     electronic wavefunction as it approaches the metal surface       
c     The wf is propagated on Coulomb Wave DVR (radial) and Legendre
c     DVR (angular)
c     Several options are possible:
c     1) calculate using a constant velocity for the incidence
c     2) calculate using the mean-field approximation: the positive ion
c        travels under a mean-potential described by the electronic
c        wavefunction
c     3) can described incidence at Jellium aluminium surface, or at
c        Cu(100) or Cu(111) surface (using Chulkov pseudopotential) 
c        just comment or uncomment lines in potsurf subroutine
c     4) can easily uncomment out pseudo-potential for xenon calcs
c        (search pseudo and the terms should appear)
c     see Eric So thesis
      implicit double precision (a-h,o-z)
      parameter (mm = 3 )      ! the number of states to run calcs for  !!!!!!!!!!!!!!!!!!REMEMBER TO CHANGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      parameter (nmax = 2000)  ! the upper limit of the number of CWDVR radial points
      dimension n0(mm)          
      dimension k0(mm)         
      dimension r1(nmax)
      common /options/ ntest,meanfield,mplot  !options for running program
      common /mass/ xmass      ! mass of ion core
      common /array/ nr,nt     ! number of radial and angular points in propgation
      common /param/ pi,dk,z,zacc,ef,dt,finish !mix of parameters
      common /wave/ xlmin,xlmax,dx,eta1 !some parameters for CWDVR
      common /optics/ eminr,rabs,rmid   !parameters for absorbing potential
      common /teq0/ dist0,vel0          !initial conditions at t=0
      common /state/ ni                 !n of interest 
      common /interp/ b,nrb,ntb         !initial diagonalisation grid parameters 
      common /jellium/ z0,v0,beta,ajel,bjel !Jellium surface potential parameters
      common /outputs/ nwstep,npstep,nfstep 
      common /chulpot/ a10,a1,as,z1,a20,a2,bs,zim,a3,xlam,aa
      pi = dacos(-1.d0)
c
cccc  OPTIONS --------------------------------
      ntest = 0    ! ntest=1 stops calculation after initial diagonalisation, use this option if you want to find the state number first!
                               ! ntest=0 diagonalises and the propagates the wf      
      meanfield = 1! meanfield = 1 for meanfield calc, meanfield=0 for constant velocity calc
      mplot = 1    ! output wavefunction for plotting? (note: output files can get big! don't just plot out wf willy nilly!)
c---------------------------------------------
cccc  Parameters that need modifying
c
c     state parameters
      ni = 3       ! the principal quantum number of interest
c     the number of the state from the initial diagonalisation required,
c     if you don't know, run with ntest = 1 to find the state number first
      n0(1) = 4     !state 1
      n0(2) = 5     !state 1
      n0(3) = 6     !state 1
c      
c     the parabolic quantum number of the state of interest, required
c     for finding out the energy shift from infinite distance and thus
c     the change of velocity from v_infinite
      k0(1) = -2    !state 1
      k0(2) =  0    !state 2
      k0(3) = +2    !state 3
c
c     initial conditions
      dist0 = 6.d0*ni**2  ! initial distance to start propagation (where the initial diagonalisation is carried out)
      vel0 = -3.d-4       ! velocity at inifinite distance or the velocity for constant vel calc
      ef = 0.d0           ! applied electric field +ve for ion-extraction field, -ve for electron-extraction
      xmass = 1836.d0     ! mass of ion core
c      
c     initial diagonalisation grid parameters (regularised Laguerre)
      b = 0.3d0    ! radial scaling parameter
      nrb = 40     ! no. of radial points for initial diagonalisation
      ntb = 20     ! no. of angular points for initial diagonalisation
c      
c     absorbing potential parameters
      rmid = dist0 !absorbing boundary position on radial grid, AND where flux detector plane sits, it is recommended to set it to dist0
      rabs = 30.d0 !width of absorbing boundary, needs to be wide enough to absorb the lowest energy components
c
c     CWDVR parameters (the propagation grid)
      nt = 20      !no. of cwdvr angular points, can be greater than ntb if desired
      dk = 3.0d0   !coulomb wave parameter: inc. dk-> 1.more points, 2.smaller sep further out, 3.more even distribution at larger dist
      z = 50.d0    !inc.  z-> smaller the first grid point                
      rmax = rmid+rabs  !maximum radial point of the grid
c
c     parameters for initial CWDVR grid point search
      zacc = 1.d-8 ! newton-raphson accuracy
      xmin = 8.d-3 ! Lower bound of 1st zero
      nx = 100000  ! Number of grid points to scan for zero, may need to increase this if using high dk parameter    
c
c     propagation parameter
      dt = 1.0d0   ! timestep
      finish = 1e-2! stop when population less than 'finish'
c
c     outputting parameters
      nwstep = 100 ! number of time steps (+1)  between each output
      nfstep = 1   ! number of time steps (+1) between evaluation of meanfield
      npstep = 500 ! if outputting wf, the time steps between succesive outputs
      npt = nt     ! no. of angular points in outputting the wavefunction, setting more than nt will give interpolated values
c      
c     END OF PARAMETERS THAT REQUIRE "ROUTINE" MODIFICATION
c---------------------------------------------------------------------------------------------------------------------------
c------------------------------------------------------------------
c     absorbing potential parameters
      pminr = 2.d0*pi/rabs
      eminr = 0.5d0*pminr**2
c
c     CWDVR parameters and parameters for scanning CWDVR points
      eta1 = -z/dk
      xmax = rmax*dk      
      dx = (xmax-xmin)/nx 
      xlmin = 0
      xlmax = 0
c
c     Jellium surface parameters for aluminium surface (see paper by Jennings and Jones)
      z0 = 0.7d0
      v0 = -0.574d0
      beta = 1.25d0
      ajel = -1.d0-4.d0*v0/beta
      bjel = -2.d0*(v0/ajel)
c
c     chulkov potential parameters for metal surface of interest
      xhar = 27.2113845d0
c
c-----Cu100---------------------      
      as = 3.415d0
      a10 = -11.480d0/xhar
      a1 = 6.1d0/xhar
      a2 = 3.782d0/xhar
      bs = 2.539d0
c      
c-----Cu111--------------
c      as = 3.94d0
c      a10 = -11.895d0/xhar
c      a1 = 5.14d0/xhar
c      a2 = 4.3279d0/xhar
c      bs = 2.9416d0
c------------------------
      z1 = 5.d0*pi/(4.d0*bs)
      a20 = a2-a10-a1
      a3 = -a20-a2/dsqrt(2.d0)
      aa = bs*a2*dsin(bs*z1)/a3
      xlam = 2.d0*aa
      zim = z1-1.d0/aa*dlog(-aa/(2.d0*a3))
c
c..........................CALCULATIONS....BEGIN....................................
c
c     get CWDVR zeros roughly, by scanning r and finding when the sign of the coulomb function changes
      call couzero(r1,nr,nmax,xmin,nx)
      if (meanfield.eq.1) then
         print*,'Calculation with Mean-Field Approx.'
      else
         print*,'Calculation with Constant Velocity Approx'
      endif
      print*,'P A R A M E T E R S :'
      print*,'====================='
      write(6,'(A,i5.2,/)')'hydrogen n = ',ni
      write(6,'(A,/,A,/,A,i5.2,4X,A,f5.2,4X,A,i5.2,/)')
     &'initial diagonalisation with lag DVR:',
     &'-------------------------------------',
     &'nr = ',nrb,'b = ',b,'nt = ',ntb
c
      write(6,'(A,/,A,/,A,f5.2,4X,A,f5.2,/)')
     &'cwdvr grid parameters:',
     &'----------------------',
     &'dk = ',dk,'z = ',z
c
      write(6,'(A,/,A,/,A,i5.2,4X,A,i5.2,4X,A,f5.2,4X,A,i5.2,4X,A,f5.2,
     & 4X,A,f5.2/)')
     &'wpp:',
     &'----',
     &'number of cwdvr pts = ',nr,'angular points nt = ',nt,'dt = ',dt,
     & 'time steps between outputs=',npstep,
     &'velocity at infinite distance= ',vel0,'E-field=',ef
c      
      call main(r1,nmax,npt,n0,k0,mm)
      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine main(r1,nmax,npt,n0,k0,mm)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      common /param/ pi,dk,z,zacc,ef,dt,finish
      common /wave/ xlmin,xlmax,dx,eta1
      common /optics/ eminr,rabs,rmid
      common /teq0/ dist0,vel0
      common /options/ ntest,meanfield,mplot  
      common /interp/ b,nrb,ntb
      common /outputs/ nwstep,npstep,nfstep 
      dimension n0(mm),r1(nmax),r(nr),a(nr),wt(nr)
      dimension k0(mm),velacc(mm)
      dimension btg(nt,nt),cost(nt)
      dimension ptg(nt,npt),pcost(npt)
      dimension hwav(nrb*ntb,mm),basis(nr,nrb)
      double complex psi(nr,nt),etr(nr,nr,nt),ev2(nr,nt)
      character*17 filen1
      character*18 filen2
      character*15 filen3
      character*8 filen4
      character*13 filen6
      character*17 filen5
      double complex forcez
      dimension zx(nr,nt)
      dimension voptic(nr)
c     setup
      t0 = dsecnd()           !clear timer
      do i = 1,nr             !fill r with the approximate zeroes calculated previously
         r(i) = r1(i)
      enddo
c     Use Newton-Raphson method (differential required) to find exact zeroes and radial weights
      call rtnewt(r,a,wt)
c 
c     work out operators, exponentiate Kinetic operator, work out flux operator
      call xoperats(r,etr,cost,btg,rmid,ntest,
     +             mrflag,voptic)
c
c     work out angular FBR function at npt plotting points
      if (mplot.eq.1) then
         call anginterp(npt,ptg,pcost)
      endif
c
c     work out electron z position for all r , cos(theta)
      do i = 1,nr
         do j = 1,nt
            zx(i,j) = r(i)*cost(j)
         enddo
      enddo
c
c     set over states to propagate
      do istat = 1,mm
         nw = n0(istat)
c        get intial wavefunction
         if (istat.eq.1) then
            call init(btg,ef,psi,hwav,b,r,nr,nt,a,
     &           nrb,n0,basis,ntb,velacc,k0,mm) !diagonalise in regularised-laguerre DVR
c           hwav contains all the wanted mm eigenvectors from initial diagonalisation            
         else
c           no need to digaonlise again just project next state in hwav onto CWDVR and propagate 
            call wavfunc(btg,psi,hwav,b,nr,nt,a,nrb,basis,ntb,istat,mm)
         endif
c         
c     open files for output
         write(filen1,'(A,i4.4)')'forward_flux.',nw
         write(filen2,'(A,i4.4)')'backward_flux.',nw
         write(filen4,'(A,i4.4)')'pop.',nw
         open(100+nw,file=filen1,status='unknown')
         open(200+nw,file=filen2,status='unknown')
         open(400+nw,file=filen4,status='unknown')
         if (mplot.eq.1) then
            write(filen5,'(A,i4.4)')'wavefunction.',nw
            open(500+nw,file=filen5,status='unknown')
            write(500+nw,*),nr,npt,dt*npstep,rmid
         endif
         if (meanfield.eq.1) then
            write(filen3,'(A,i4.4)')'totalforce.',nw
            write(filen6,'(A,i4.4)')'velocity.',nw
            open(300+nw,file=filen3,status='unknown')
            open(600+nw,file=filen6,status='unknown')
         endif
c     
         if (istat.eq.1) then
            t1 = dsecnd()
            write (6,61) (t1-t0)/60.d0
  61        format(/1x,'Setup completed in:',f9.2,' minutes'//1x) 
            print*,'Time....................Population...........
     +........Dist'
         endif
c         
c        propagation
         if (meanfield.eq.1) then
            vel = +velacc(istat) !include energy shift acceleration
         else
            vel = vel0 !assume velocity at infinite distance
         endif
         t0 = t1
         ptold = 10.d0
         istep = -1
         dist = dist0
         ionstep = 0
         dedz = 0.d0
         dft = nfstep*dt
         p0 = 1.d0
1        istep = istep+1
         time = istep*dt
         if (istep .gt. 0) then
c        work out time-depedent potential energy operator
            call potop (ev2,dt,r,cost,time,dist,ef,zx)
c        propagate wavefunction with split operator
            call split (psi,etr,ev2,btg)
         endif
c
c        work out population on grid at time t
         pt = pnorm(psi)
c
c        output
         writstep = mod(istep,nwstep)
         if (writstep.eq.0.d0) then
c           output flux
            call flux2(psi,cost,mrflag,voptic,fg2,fl2,vel)
            write(100+nw,*)dist,fl2
            write(200+nw,*)dist,fg2
c           work out wavefunction population inside area bound by absorbing boundary and surface at time t
c            pt2 = pnorm2(psi,dist,r,cost,rmid)
            write(400+nw,*) dist,pt
            write (6,*) time,pt,dist
         endif
c        plot wavefunction
         if (mplot.eq.1) then
            plotstep = mod(istep,npstep)
            if (plotstep.eq.0.d0) then
               call plot(psi,r,wt,nw,npt,btg,ptg,pcost,dist)
            endif
         endif
c
          if (meanfield.ne.1) then
             dist = dist + vel*dt
             goto 426
          endif
c        work out force on ion-core and propgate with velocity verlet
         fstep = mod(istep,nfstep)
         if (fstep.eq.0.d0) then
            call denergy2(psi,r,dist,ef,zx,dedz)
            if (writstep.eq.0.d0) then
                write(600+nw,*) dist,vel,time
                write(300+nw,*) dist,dedz
            endif
            call fion (psi,dist,dedz,dft,vel,zx,ionstep,dedzo)
            ionstep = 1
         endif
c
         dist = dist + vel*dt -0.5d0*dedz*dt**2
c         
c        Stop if ion goes backwards too far!
426      if (dist.gt.(dist0+50.d0)) stop 'reversed too far!'
c
c        WARN if wavepacket grows
         if (pt-ptold.gt.1.d-8) print*, 'growing!'
         ptold = pt
c
c        stop when distance from surface =0
         if (dist.lt.0.d0) go to 847
         if (pt .gt. finish) go to 1
847      t1 = dsecnd()
         write (6,63) (t1-t0)/60.d0
  63     format(/1x,'Propagation took:',f24.2,' minutes')
         close(100+nw)
         close(200+nw)
         close(300+nw)
         close(400+nw)
         close(500+nw)
         close(600+nw)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plot(psi,r,wt,nw,npt,btg,ptg,pcost,dist)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      double complex psi(nr,nt),phi(nr,nt),pxi(nr,npt)
      dimension r(nr),wt(nr),btg(nt,nt)
      dimension pcost(npt),ptg(nt,npt)
c     coeffs: angular DVR-> angular FBR conversion
      r0 = 0.d0
      r1 = 1.d0
      nrs = 2*nr
      call dgemm('n','t',nrs,nt,nt,r1,psi,nrs,btg,nt,r0,phi,nrs)
c     coeff x FBR basis functions (evaluated at npt angular points)
      r0 = 0.d0
      r1 = 1.d0
      nrs = 2*nr
      call dgemm('n','n',nrs,npt,nt,r1,phi,nrs,ptg,nt,r0,pxi,nrs)
      do i = 1,nr
         do j = 1,npt
            si = dble((pxi(i,j))*dconjg(pxi(i,j)))/(wt(i)) !*pi/(wt(i))
            zz = r(i)*pcost(j)
            ro = r(i)*dsqrt(1.d0-(pcost(j)**2))
            if ((i.eq.1).and.(j.eq.1)) then
               write(500+nw,*) ro,zz,si,dist
            else
               write(500+nw,*) ro,zz,si,0 !saves space, only need dist once
            endif
         enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine expmih (etr,tr,r,cent,dt,voptic)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      common /optics/ eminr,rabs,rmid
      double complex etr,vopt
      dimension etr(nr,nr,nt)
      dimension tr(nr,nr),r(nr)
      dimension cent(nt)
      dimension voptic(nr)
c     -----------------------------------------------------------------
c     Sets up the exponentiated Hamiltonian factors
c     needed for split operator propagation.
c     -----------------------------------------------------------------
c
c     exp(-iTdt)
      do k = 1,nt
         do j = 1,nr
            voptic(j) = 0.d0       
            do i = 1,nr
               etr(i,j,k) = dcmplx(0.d0,-dt)*tr(i,j)
            enddo
            rmrmid = r(j)-rmid
            call optpot (rmrmid,rabs,eminr,vopt)
            if (rmrmid .gt. 0.d0) then
                voptic(j) = dimag(vopt)
            endif
            vopt = tr(j,j)+(cent(k)/r(j)**2)+vopt!+pseudopot(k-1,r(j))
            etr(j,j,k) = dcmplx(0.d0,-dt)*vopt
         enddo
         call expmat(etr(1,1,k),nr)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine expmat(z,n)
      implicit double precision (a-h,o-z)
      double complex z,zeta,v,work,s
c     -----------------------------------------------------------------
c     Exponentiates a complex symmetric matrix z to give exp(z).
c     -----------------------------------------------------------------
      dimension z(n,n),zeta(n)
      dimension v(n,n),work(2*n),rwork(2*n)
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
         zr = dreal(zeta(k))
         zi = dimag(zeta(k))
         zeta(k) = exp(zr)*dcmplx(cos(zi),sin(zi))
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
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine optpot (rmrmid,rabs,eminr,vopt)
      implicit double precision (a-h,o-z)
      double complex vopt
c     -----------------------------------------------------------------
c     Transmission-free negative imaginary absorbing potential.
c     -----------------------------------------------------------------
      parameter (c = 2.62206d0)
      vopt = (0.d0,0.d0)
      if (rmrmid .gt. 0.d0) then
         x = c*rmrmid/rabs
         y = 4.d0/(c-x)**2+4.d0/(c+x)**2-8.d0/c**2  
         vopt = dcmplx(0.d0,-eminr*y)
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xoperats(r,etr,cost,btg,rmid,ntest,
     +         mrflag,voptic)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      common /param/ pi,dk,z,zacc,ef,dt,finish
      dimension r(nr),tr(nr,nr)
      dimension btg(nt,nt),cost(nt),cent(nt)
      double complex etr(nr,nr,nt)
      dimension voptic(nr)
c 
c     radial dvr kinetic operator
c
      call coutr(r,nr,dk,z,tr)
c
c     angular grid points,transform matrix,weights,centrifugal
c
      call legdvr(nt,cost,btg,cent)
c
      if (ntest.ne.1) then
c        exponentiate kinetic energy operator
         print*,'exponentiating KE operator'
         call expmih(etr,tr,r,cent,dt,voptic)
c   
c        work out absorbing boundary limit
         do ii = 1,nr
            if (r(ii).ge.rmid) then
               mrflag = ii
               goto 555
            endif
         enddo
555   endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine couzero(r,nr,nmax,xmin,nx)
      implicit double precision (a-h,o-z)
c     -----------------------------------------------------------------
c     approx location of coulomb function, F_0(eta,x), zeros for DVR grid.
c     -----------------------------------------------------------------
      dimension r(nmax)
      common /wave/ xlmin,xlmax,dx,eta1
c
c     First find approximate location of zeros by stepping
      nr = 0
      ifail = 0
      call coulfg(xmin,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
      if (ifail.ne.0) stop 'couzero | coulfg 1'
      fo = f 
c
      do i = 1,nx+1
         x = xmin+i*dx 
         ifail = 0
         call coulfg(x,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
         if (ifail.ne.0) stop 'couzero | coulfg 2'
         if (f*fo.lt.0.d0) then
            nr = nr+ 1   
            if (nr.gt.nmax) stop 'couzero | nr > nmax'
            r(nr) = x - 0.5d0*dx
         endif
         fo = f
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rtnewt(r,a,wt)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      common /param/ pi,dk,z,zacc,ef,dt,finish
      dimension r(nr),a(nr),wt(nr)
      common /wave/ xlmin,xlmax,dx,eta1
c
c     Newton-Raphson method to find roots more accurately. 
      do i = 1,nr
         x1 = r(i) - 0.5d0*dx 
         x2 = r(i) + 0.5d0*dx 
         xz = r(i)
         j=0
123      j=j+1
         ifail = 0
         call coulfg(xz,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
         if (ifail.ne.0) stop 'couzero | coulfg 3'
         dz = f/fp 
         xz = xz - dz
         if ((x1-xz)*(xz-x2).lt.0.d0)
     +        pause 'jumped out of bracket'
         if (abs(dz).gt.zacc) go to 123
c         print*, 'newtr converged after',j,'iterations'
         call coulfg(xz,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
         if (ifail.ne.0) stop 'couzero | coulfg 3'
         a(i) = fp*sqrt(dk)
         wt(i) = pi/(dk*(fp**2))
         r(i) = xz/dk
      enddo
c     test if seperation of grid points is ok
      slimit = pi/dsqrt(2.d0*0.574d0)
      xmaxdif = 0.d0
      ncount = 0
      do i = 1,nr-1 
         dif = dabs(r(i+1)-r(i))
         if (dif.gt.xmaxdif) xmaxdif = dif
         if (dif.gt.slimit) then
            ncount = ncount+1
           !write(6,'(A,i9.2,A,i9.2,A)')'point',i,'/',nr,'too far apart'
         endif 
      enddo
      write(6,'(A,f9.4)')'max spacing',xmaxdif
      if (ncount.ne.0) then
         !write(6,'(A,f9.2)') 'sep need to be lt',slimit
c         pause 'need to inc. dk paramter!!!!!!!!'
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP, 
     *                  MODE1,KFN,IFAIL)
      implicit double precision (a-h,o-z)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. BARNETT           MANCHESTER  MARCH   1981                    C
C                                                                      C
C  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           C
C                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           C
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           C
C  THIS VERSION WRITTEN UP       IN    CPC XX (1982) YYY-ZZZ           C
C                                                                      C
C  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), C
C   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    C
C   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   C
C   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     C
C   THE DIRAC EQUATION ,ALSO SPHERICAL + CYLINDRICAL BESSEL EQUATIONS  C
C                                                                      C
C  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    C
C  STARTING ARRAY ELEMENT IS M1 = MAX0(  INT(XLMIN+ACCUR),0) + 1       C
C      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 C
C                                                                      C
C  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     C
C            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    C
C            = 3      F               CALL TO AT LEAST LENGTH (1)      C
C  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            C
C            = 1 SPHERICAL   BESSEL                                    C
C            = 2 CYLINDRICAL BESSEL                                    C
C  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          C
C                                                                      C
C  PRECISION&  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    C
C   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     C
C   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   C
C   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   C
C   FOR MANTISSAS OF 56 + 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) C
C   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DIMENSION    FC(1),GC(1),FCP(1),GCP(1)
      LOGICAL      ETANE0,XLTURN
      COMMON       /STEED/ PACCQ,NFP,NPQ,IEXP,M1
C***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
C***  COULFG HAS CALLS TO&  SQRT, ABS, AMOD ,  INT, SIGN, FLOAT, MIN1 
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0d0, 1.0d0, 2.0d0, 1.0d2, 2.0d4/
      DATA HALF,TM30 / 0.5d0, 1.0d-30 / 
c      DATA RT2EPI /0.79788 45608 02865/
C *** THIS CONSTANT IS  DSQRT(TWO/PI)&  USE Q0 FOR IBM REAL*16& D0 FOR
C ***  REAL*8 + CDC DOUBLE P&  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
C
c                       ACCUR = R1MACH(4)
      pi=dacos(-1.d0)
      accur = 1.d-20
      RT2EPI=dsqrt(2.d0/pi)
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      MODE  = 1
      IF(MODE1 .EQ. 2 .OR. MODE1 .EQ. 3 ) MODE = MODE1
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO
      ACC   = ACCUR 
      ACC4  = ACC*TEN2*TEN2
      ACCH  = dSQRT(ACC)
C ***    TEST RANGE OF XX, EXIT IF.LE. SQRT(ACCUR) OR IF NEGATIVE
C
      IF(XX .LE. ACCH)                          GO TO 100
      X     = XX
      XLM   = XLMIN 
      IF(KFN .EQ. 2)  XLM = XLM - HALF
      IF(XLM .LE. -ONE .OR. XLMAX .LT. XLMIN)   GO TO 105
      E2MM1 = ETA*ETA + XLM*XLM + XLM
      XLTURN= X*(X - TWO*ETA) .LT. XLM*XLM + XLM
      DELL  = XLMAX - XLMIN + ACC
      IF( ABS(MOD(DELL,ONE)) .GT. ACC) WRITE(6,2040)XLMAX,XLMIN,DELL
      LXTRA =   INT(DELL)
      XLL   = XLM +  FLOAT(LXTRA)
C ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
C ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
      M1  = MAX0(  INT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
C
C ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
    2 EK  = ETA / PK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 =  PK + ONE
C ***   TEST ENSURES B1 .NE. ZERO FOR NEGATIVE ETA; FIXUP IS EXACT.
             IF( ABS(ETA*X + PK*PK1) .GT. ACC)  GO TO 3
             FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
             PK   =  TWO + PK 
      GO TO 2
    3 D   =  ONE/((PK + PK1)*(XI + EK/PK1))
      DF  = -FCL*(ONE + EK*EK)*D
            IF(FCL .NE. ONE )  FCL = -ONE
            IF(D   .LT. ZERO)  FCL = -FCL
      F   =  F  + DF
C
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
C
      P     = ONE
    4 PK    = PK1
        PK1 = PK1 + ONE
        EK  = ETA / PK
        TK  = (PK + PK1)*(XI + EK/PK1)
        D   =  TK - D*(ONE + EK*EK)
              IF( ABS(D) .GT. ACCH)             GO TO 5
              WRITE (6,1000) D,DF,ACCH,PK,EK,ETA,X
              P = P  +   ONE
              IF( P .GT. TWO )                  GO TO 110
    5 D     = ONE/D 
              IF (D .LT. ZERO) FCL = -FCL
        DF  = DF*(D*TK - ONE) 
        F   = F  + DF
              IF(PK .GT. PX)                    GO TO 110
      IF( ABS(DF) .GE.  ABS(F)*ACC)             GO TO 4
                  NFP = PK - XLL - 1
      IF(LXTRA .EQ. 0)                          GO TO 7
C
C *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
C
      FCL = FCL*TM30
      FPL = FCL*F
      IF(MODE .EQ. 1) FCP(L1) = FPL
                      FC (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO 6  LP = 1,LXTRA
         IF(ETANE0) EL = ETA/XL
         IF(ETANE0) RL =  dSQRT(ONE + EL*EL)
         SL    =  EL  + XL*XI 
         L     =  L1  - LP
         FCL1  = (FCL *SL + FPL)/RL
         FPL   =  FCL1*SL - FCL *RL
         FCL   =  FCL1
         FC(L) =  FCL
         IF(MODE .EQ. 1) FCP(L)  = FPL
         IF(MODE .NE. 3 .AND. ETANE0) GC(L+1) = RL
    6 XL = XL - ONE 
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
C ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM 
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM 
C
    7 IF( XLTURN ) CALL JWKB(X,ETA,MAX(XLM,ZERO),FJWKB,GJWKB,IEXP)
      IF( IEXP .GT. 1 .OR. GJWKB .GT. ONE/(ACCH*TEN2))  GO TO 9
          XLTURN = .FALSE.
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X - ETA)
      BI =  TWO
      DR =  BR/(BR*BR + BI*BI)
      DI = -BI/(BR*BR + BI*BI)
      DP = -XI*(AR*DI + AI*DR)
      DQ =  XI*(AR*DR - AI*DI)
    8 P     = P  + DP
         Q  = Q  + DQ
         PK = PK + TWO
         AR = AR + PK
         AI = AI + WI
         BI = BI + TWO
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         C  = ONE/(D*D + DI*DI)
         DR =  C*D
         DI = -C*DI 
         A  = BR*DR - BI*DI - ONE
         B  = BI*DR + BR*DI
         C  = DP*A  - DQ*B
         DQ = DP*B  + DQ*A
         DP = C
         IF(PK .GT. TA)                         GO TO 120
      IF( ABS(DP)+ ABS(DQ).GE.( ABS(P)+ ABS(Q))*ACC)   GO TO 8
                      NPQ   = PK/TWO
                      PACCQ = HALF*ACC/MIN( ABS(Q),ONE)
                      IF( ABS(P) .GT.  ABS(Q)) PACCQ = PACCQ* ABS(P)
C
C *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM 
C
      GAM = (F - P)/Q
            IF(Q .LE. ACC4* ABS(P))             GO TO 130
      W   = ONE/ dSQRT((F - P)*GAM + Q)
            GO TO 10
C *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 250 + XLTURN = .TRUE.
    9 W   = FJWKB
      GAM = GJWKB*W 
      P   = F
      Q   = ONE
C
C *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C
   10                     ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .EQ. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .EQ. 2) BETA  =  dSQRT(XI)*RT2EPI
      FCM  =  SIGN(W,FCL)*BETA
           FC(M1)  = FCM
                      IF(MODE .EQ. 3)           GO TO 11
           IF(.NOT. XLTURN)   GCL =  FCM*GAM
           IF(      XLTURN)   GCL =  GJWKB*BETA
           IF( KFN .NE. 0 )   GCL = -GCL
           GC(M1)  = GCL
           GPL =  GCL*(P - Q/GAM) - ALPHA*GCL
                      IF(MODE .EQ. 2)           GO TO 11
           GCP(M1) = GPL
           FCP(M1) = FCM*(F - ALPHA)
   11 IF(LXTRA .EQ. 0 ) RETURN
C *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
C *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
         W    = BETA*W/ ABS(FCL)
         MAXL = L1 - 1
      DO 12 L = M1,MAXL
                      IF(MODE .EQ. 3)           GO TO 12
                      XL = XL + ONE
         IF(ETANE0)   EL = ETA/XL
         IF(ETANE0)   RL = GC(L+1)
                      SL = EL + XL*XI
         GCL1     = ((SL - ALPHA)*GCL - GPL)/RL
         GPL      =   RL*GCL -  (SL + ALPHA)*GCL1 
         GCL      = GCL1
         GC(L+1)  = GCL1
                      IF(MODE .EQ. 2)           GO TO 12
         GCP(L+1) = GPL
         FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
   12 FC(L+1)     = W* FC(L+1)
      RETURN
 1000 FORMAT(/' CF1 ACCURACY LOSS& D,DF,ACCH,K,ETA/K,ETA,X = ',1P7E9.2/)
C
C ***    ERROR MESSAGES
C
  100 IFAIL = -1
      WRITE(6,2000) XX,ACCH
 2000 FORMAT(' FOR XX = ',1PE12.3,' TRY SMALL-X  SOLUTIONS',
     *' OR X NEGATIVE'/ ,' SQUARE ROOT ACCURACY PARAMETER =  ',E12.3/)
      RETURN
  105 IFAIL = -2
      WRITE (6,2005) XLMAX,XLMIN,XLM
 2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES&XLMAX,XLMIN,XLM = ',
     *1P3E15.6/)
      RETURN
  110 IFAIL =  1
      WRITE (6,2010) ABORT,F ,DF,PK,PX,ACC
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',/
     *' F,DF,PK,PX,ACCUR =  ',1P5E12.3//)
      RETURN
  120 IFAIL =  2
      WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',/
     *' P,Q,DP,DQ,ACCUR =  ',1P4E17.7,E12.3//)
      RETURN
  130 IFAIL =  3
      WRITE (6,2030) P,Q,ACC,DELL,LXTRA,M1
 2030 FORMAT(' FINAL Q.LE. ABS(P)*ACC*10**4 , P,Q,ACC = ',1P3E12.3,4X,
     *' DELL,LXTRA,M1 = ',E12.3,2I5 /)
      RETURN
 2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P3E20.10/)
      END 
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
      REAL      XX,ETA1,XL,FJWKB,GJWKB, ZERO
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0
C *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C *** CALLS AMAX1,SQRT,ALOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0E0, 0.5E0, 1.0E0, 6.0E0, 10.0E0 /
      DATA  DZERO, RL35, ALOGE  /0.0E0, 35.0E0, 0.43429 45 E0 /
      X     = XX
      ETA   = ETA1
      GH2   = X*(ETA + ETA - X)
      XLL1  = MAX(XL*XL + XL,DZERO)
      IF(GH2 + XLL1 .LE. ZERO) RETURN
       HLL  = XLL1 + SIX/RL35 
       HL   = SQRT(HLL)
       SL   = ETA/HL + HL/X
       RL2  = ONE + ETA*ETA/HLL
       GH   = SQRT(GH2 + HLL)/X
       PHI  = X*GH - HALF*( HL*LOG((GH + SL)**2/RL2) - LOG(GH) )
          IF(ETA .NE. ZERO) PHI = PHI - ETA*ATAN2(X*GH,X - ETA)
      PHI10 = -PHI*ALOGE
      IEXP  =  INT(PHI10)
      IF(IEXP .GT. 250) GJWKB = TEN**(PHI10 - FLOAT(IEXP))
      IF(IEXP .LE. 250) GJWKB = EXP(-PHI)
      IF(IEXP .LE. 250) IEXP  = 0
      FJWKB = HALF/(GH*GJWKB) 
      RETURN
      END 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coutr(r,n,dk,z,t)
      implicit double precision (a-h,o-z)
c     -----------------------------------------------------------------
c     Coulomb dvr kinetic energy
c     -----------------------------------------------------------------
      dimension r(n),t(n,n)
      te = dk**2
      do i = 1,n
         ri = r(i)
         do j = 1,n
            rj = r(j) 
            t(i,j) = 1.d0/(ri-rj)**2 
         enddo 
         t(i,i) =(te + 2.d0*z/ri)/6.d0 
      enddo
c     add diag coulomb potential and for case of XE, bardsley pseudo-potential
c      ad = 4.044d0
c      aq = 14.235d0
c      d = 1.d0
      do i = 1,n
         t(i,i) = t(i,i)-1.d0/r(i) 
c     +                            -ad*0.5d0/((r(i)**2+d**2)**2)
c     +                           -aq*0.5d0/((r(i)**2+d**2)**3)
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine symevp (a,n,d,ierr)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Uses LAPACK DSYEV to diagonalise a real symmetric matrix.
c     ----------------------------------------------------------------- 
c
      dimension a(n,n),d(n)
      dimension work(34*n)
c
      lwork = 34*n 
      call dsyev ('V','U',n,a,n,d,work,lwork,ierr) 
c
c     normalise eigenvectors
      do j= 1,n
         s = 0.d0
         do i = 1,n
            s = s+a(i,j)**2
         enddo
         s = dsqrt(s)
         do i = 1,n
            a(i,j) = a(i,j)/s
         enddo
      enddo
      return
      end
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine legdvrw(n,x,btg,twt,cent)
      implicit double precision (a-h,o-z)
c     -----------------------------------------------------------------
c     Legendre polynomial discrete variable representation.
c     -----------------------------------------------------------------
      dimension x(n),d(n),btg(n,n),work(2*(n-1)),twt(n),cent(n)
c     compute jacobi matrix
      do i = 1,n
         x(i) = 0.d0
         d(i) = (dble(i)/(dsqrt(4.d0*(i**2)-1.d0)))
      enddo
c     compute zeros (eigenvalues)
      ierr = 0
      call dstev('V',n,x,d,btg,n,work,ierr)
      if (ierr.ne.0) stop 'lagrt | dgtsv'
c     work out legendre weights
      ovnorm = dsqrt(2.d0)
      do i = 1,n
         if (btg(1,i).lt.0.d0) then
          do j = 1,n
             btg(j,i) = btg(j,i)*(-1.d0)
          enddo
         endif
         twt(i) = (btg(1,i)*ovnorm)**2
      enddo
c     work out centrifugal term
      j = 0
      do i = 1,n
         cent(i) = 0.5d0*j*(j+1)
         j = j+1
      enddo
      return
      end
      subroutine legdvr(n,x,btg,cent)
      implicit double precision (a-h,o-z)
c     -----------------------------------------------------------------
c     Legendre polynomial discrete variable representation.
c     -----------------------------------------------------------------
      dimension x(n),d(n),btg(n,n),work(2*(n-1)),cent(n)
c     compute jacobi matrix
      do i = 1,n
         x(i) = 0.d0
         d(i) = (dble(i)/(dsqrt(4.d0*(i**2)-1.d0)))
      enddo
c     compute zeros (eigenvalues)
      ierr = 0
      call dstev('V',n,x,d,btg,n,work,ierr)
      if (ierr.ne.0) stop 'lagrt | dgtsv'
c     work out legendre weights
      ovnorm = dsqrt(2.d0)
      do i = 1,n
         if (btg(1,i).lt.0.d0) then
          do j = 1,n
             btg(j,i) = btg(j,i)*(-1.d0)
          enddo
         endif
      enddo
c     work out centrifugal term
      j = 0
      do i = 1,n
         cent(i) = 0.5d0*j*(j+1)
         j = j+1
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function gaus(a,pi,a0,width)
      implicit double precision(a-h,o-z)
      anorm=dsqrt(1.d0/(width*dsqrt(pi)))
      ex=dexp(-((a-a0)**2)/(2.d0*(width**2)))
      gaus=anorm*ex
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potop (ev2,dt,r,cost,time,dist,ef,zx)
      implicit double precision (a-h,o-z)
c     -------------------------------------------------------------------------
c     compute the time-dependent potential energy operator
c     -------------------------------------------------------------------------
      double complex ev2
      common /array/ nr,nt
      dimension ev2(nr,nt),r(nr),cost(nt)
      dimension zx(nr,nt)
c     exp(-iVdt/2)
      dt2 = 0.5d0*dt
      do i = 1,nr
         do j = 1,nt
            zz = zx(i,j)
            zd = zz + dist
            if (zd.lt.0.d0) then
               v = potsurf(r(i),zz,dist)
            else
               v = potsurf(r(i),zz,dist)
     +             +ef*zd
            endif
            ar = -v*dt2
            ev2(i,j) = dcmplx(dcos(ar),dsin(ar))
         enddo
      enddo
      return
      end
 
      double precision function dpotsurf(r,z,d)
      implicit double precision (a-h,o-z)
      common /jellium/ z0,v0,beta,ajel,bjel
      zz = z+d
c     electron-image electron attraction
      if (zz.gt.z0) then
        vee = 0.25d0*(1.d0-dexp(-beta*(zz-z0)))/(zz-z0)**2
     +       - 0.25d0*beta*dexp(-beta*(zz-z0))/(zz-z0)
      else
        vee = -v0/(ajel*dexp(bjel*(zz-z0))+1.d0)**2
     +       * ajel*bjel*dexp(bjel*(zz-z0))
      endif
      if (zz.gt.0) then
        xl = dsqrt(r**2 + 4.d0*d**2 + 4.d0*d*z)
        vep = (-4.d0*d-2.d0*z)/xl**3
      else
        vep = 0.d0
      endif
      dpotsurf = vee+vep
      return
      end

      double precision function potsurf(r,z,d)
      implicit double precision (a-h,o-z)
      common /jellium/ z0,v0,beta,ajel,bjel
      common /chulpot/ a10,a1,as,z1,a20,a2,bs,zim,a3,xlam,aa
c     -------------------------------------------------------------------------
c     compute the surface potential at distance d from surface (surface is at z=-d)
cc     -------------------------------------------------------------------------
      zz = z+d
      y2 = (r**2-z**2)
c     electron-image electron attraction for Chulkov pseudopotentials
c      pi = dacos(-1.d0)
c      if (zz.le.0.d0) then
c         vee = a10 + a1*dcos(2.d0*pi/as*(zz))
c      elseif ((zz.gt.0.d0).and.(zz.lt.z1)) then
c         vee = -a20 + a2*dcos(bs*(zz))
c      elseif ((zz.ge.z1).and.(zz.le.zim)) then
c         vee = a3*dexp(-aa*(zz-z1))
c      elseif (zz.gt.zim) then
c         vee = (dexp(-xlam*(zz-zim))-1.d0)/(4.d0*(zz-zim))
c      endif
c      
cc    electron-image electron attraction
      if (zz.gt.z0) then
       vee = 0.25d0*(-1.d0+dexp(-beta*(zz-z0)))/(zz-z0)
      else
       vee = v0/(ajel*dexp(bjel*(zz-z0))+1.d0)
      endif
cc    electron-image proton repulsion
      xl = dsqrt(y2+(d+dabs(zz))**2)
c      if (zz.gt.0.d0) then
         vep = +1.d0/xl
c      else
c         vep = +1.d0/dsqrt(d**2+y2)
c      endif
      potsurf = vee+vep
      return
      end
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine split (psi,etr,ev2,btg)
      implicit double precision (a-h,o-z)
      double complex psi,etr,ev2
      double complex phi,c0,c1
c     -----------------------------------------------------------------
c     Evolves the wavepacket through a time step dt
c     using the symmetric split operator method.
c     -----------------------------------------------------------------
c     common block
      common /array/ nr,nt
c     input arrays
      dimension psi(nr*nt)
      dimension etr(nr*nr*nt)
      dimension ev2(nr*nt)
      dimension btg(nt*nt)
c     local array
      dimension phi(nr*nt)
c     multiplication by exp(-iVdt/2)
      n = nr*nt
      do ij = 1,n
         psi(ij) = ev2(ij)*psi(ij)
      enddo
c     angular grid to basis transformation
      r0 = 0.d0
      r1 = 1.d0
      nrs = 2*nr
      call dgemm ('n','t',nrs,nt,nt,r1,psi,nrs,btg,nt,r0,phi,nrs)
c
c     multiplication by exp(-iTdt)
c
      c0 = (0.d0,0.d0)
      c1 = (1.d0,0.d0)
      ns = 1
C$OMP PARALLEL DO
C$OMP& SHARED(nt,nr,c0,c1,etr,phi,psi)
C$OMP& PRIVATE(k,ke,kp)
      do k = 1,nt
         ke = (k-1)*nr*nr + 1
         kp = (k-1)*nr + 1
         call zgemv
     +        ('n',nr,nr,c1,etr(ke),nr,phi(kp),1,c0,psi(kp),1)
c         call zgemm
c     +   ('n','n',nr,ns,nr,c1,etr(ke),nr,phi(kp),nr,c0,psi(kp),nr)
      enddo
C$OMP END PARALLEL DO
c
c     angular basis to grid transformation
c
      nrs = 2*nr
      call dgemm ('n','n',nrs,nt,nt,r1,psi,nrs,btg,nt,r0,phi,nrs)
c
c     multiplication by exp(-iVdt/2)
c
      do ij = 1,n
         psi(ij) = ev2(ij)*phi(ij)
      enddo
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fion (psi,dist,dedz,dft,vel,zx,ionstep,dedzo)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      common /mass/ xmass
      common /teq0/ dist0,vel0
      double complex psi
      dimension psi(nr*nt)
      dimension zx(nr*nt)
      v0 = vel
      if (ionstep.ne.0) then
        v1 = v0 + 0.5d0*dedzo*dft/xmass
      else
        v1 = v0
      endif
      v5 = v1 + 0.5d0*dedz*dft/xmass
      vel = v5
c
      n = nr*nt
      do i = 1,n
         arg = (-vel+v0)*zx(i)
         psi(i) = psi(i)*dcmplx(dcos(arg),dsin(arg))
      enddo
c
      dedzo = dedz
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine denergy2(psi,r,dist,ef,zx,dedz)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      double complex psi
      dimension psi(nr,nt)
      dimension r(nr)
      dimension zx(nr,nt)
      dedz = 0.d0

      vint = 0.d0
      do i = 1,nt
         do j = 1,nr
            zz = zx(j,i)
            zd = zz + dist
            if (zd.lt.0.d0) then
               ev = dpotsurf(r(j),zz,dist)
            else
               ev = dpotsurf(r(j),zz,dist)
     +                     + ef
            endif
            vint =  vint + dble(dconjg(psi(j,i))*ev*psi(j,i))
        enddo
      enddo
      dedz =  -(-ef +0.25d0/(dist**2) + vint) 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function pnorm (psi)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      double complex psi
      dimension psi(nr,nt)
c     -----------------------------------------------------------------
c     Wavepacket left on grid
c     -----------------------------------------------------------------
      pnorm = 0.d0
      do i = 1,nr
         do j=1,nt
            pnorm = pnorm+(dble(psi(i,j))**2+dimag(psi(i,j))**2)!*pi!*a(i)**2
         enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function pnorm2 (psi,dist,r,cost,rmid)
      implicit double precision (a-h,o-z)
      common /array/ nr,nt
      double complex psi
      dimension r(nr),cost(nt)
      dimension psi(nr,nt)
c     -----------------------------------------------------------------
c     Wavepacket left inside area bound by absorbing potential and surface
c     -----------------------------------------------------------------
      d = -dist
      pnorm2 = 0.d0
      do i = 1,nr
         do j = 1,nt
            z = r(i)*cost(j)
            if ((z.gt.d).and.(r(i).lt.rmid)) then
c            if (r(i).lt.rmid) then
               pnorm2 = pnorm2+
     &          (dimag(psi(i,j))**2+dble(psi(i,j))**2)
            endif
         enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function en(xn,xm,xk,f)
      implicit double precision (a-h,o-z)
      a=-1.d0/(2.d0*(xn**2))
      b=3.d0*f*xn*xk/2.d0
      c=-(f**2)*(xn**4)*(17.d0*(xn**2)-3.d0*(xk**2)
     +   -9.d0*(xm**2)+19.d0)/16.d0
      d=(3.d0/32.d0)*(xn**7)*xk*(23.d0*(xn**2)-(xk**2)
     +   +11.d0*(xm**2)+39.d0)*(f**3)
      e=-(xn**10)*(f**4)*(5487.d0*(xn**4)+35182.d0*(xn**2)
     +  -1134.d0*(xm**2)*(xk**2)+1806.d0*(xn**2)*(xk**2))
     +   /1024.d0
      g=-(xn**10)*(f**4)*(-3402.d0*(xn**2)*(xm**2)+
     +  147.d0*(xk**4)-549.d0*(xm**4)+5754.d0*(xk**2)-
     +  8622.d0*(xm**2)+16211.d0)/1024.d0
      h=3.d0*(xn**13)*xk*(f**5)*(10563.d0*(xn**4)+90708.d0*
     +  (xn**2)+220.d0*(xm**2)*(xk**2)+98.d0*(xn**2)*
     +  (xk**2))/1024.d0
      p=3.d0*(xn**13)*xk*(f**5)*(772.d0*(xn**2)*(xm**2)-
     +  21.d0*(xk**4)+725.d0*(xm**4)+780.d0*(xk**2)+
     +  830.d0*(xm**2)+59293.d0)/1024.d0
      en=a+b+c+d+e+g+h+p
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function pseudopot(l,r)
      implicit double precision (a-h,o-z)
c     -------------------------------------------------------------------------- 
c     Bardsley Pseudopotential for Xenon
c     ---------------------------------------------------------------------- 
      if (l.eq.0) then
c       a=9.102d0
c       b=0.511d0
           a = -1.1798376d0
           b = 0.01d0
      endif
      if (l.eq.1) then
c       a=2.567d0
c       b=0.224d0
           a = -0.89700294d0
           b = 0.01d0
      endif
      if (l.eq.2) then
c       a=-0.468d0
c       b=0.190d0
           a = -0.4990644d0
           b = 0.01d0
      endif
      if (l.eq.3) then
         a = -0.0082312d0
         b = 0.01d0
      endif
      if (l.ge.4) then
         a = 0.d0
         b = 0.01d0
      endif
      pseudopot = a*exp(-b*(r**2))
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init(pbtg,ef,psi,hwav,b,x,nx,ntcw,a,
     &               nr,n0,basis,nt,velacc,k0,mm)
c--------------------------------------------------------------------------
c     intitial diagonalisation on regularised Laguerre DVR to find
c     eigenstates to propagate      
c     NB x(nx) are the cwdvr grid points and r(nr) are the laguerre points
c     cost(nt) are the ang points for initial diagonalisation and cost(ntcw) and those for propagation
c--------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /state/ ni
      common /options/ ntest,meanfield,mplot 
      common /teq0/ dist0,vel0
      dimension btg(nt,nt),cost(nt),cent(nt)
      dimension pbtg(ntcw,ntcw),ptg(nt,ntcw)
      dimension r(nr),hr(nr,nr)
      dimension h(nr*nt,nr*nt),e(nr*nt)
      dimension vt(nt,nt),ht(nt,nt,nr),phi(nx,nt)
      dimension p(nx,ntcw)
      dimension trans(nx),x(nx),a(nx),basis(nx,nr)
      dimension n0(mm), hwav(nr*nt,mm)
      double complex psi(nx,ntcw)
      dimension k0(mm),velacc(mm)
c
      n = nr*nt
c     work out nt legendre ang functions/points for initial diagonalisation
      call legdvr(nt,cost,btg,cent)
c
c     work out laguerre grid points and kinetic energy operator for regularised laguerre dvr
      call lagdvr(nr,b,r,hr)
c
c     work out L_(N)(x) with x at the nx cwdvr grid points
      call lagtrans(x,nx,trans,nr,b)
c
c     work out f_(i)(x) for i=1,...,nr laguerre basis functions at x(j) (j=1,...,nx) cwdvr points
c     BAYE REGULARISED LAGUERRE: f_(i)(x)=(-1)^(i)*x_(i)^(-0.5)*x*L_(N)(x)*e^(-x/2)/(x-x_(i))
      do j = 1,nx
         do i = 1,nr
            basis(j,i) = ((-1.d0)**i)*trans(j)
     +                   /(dsqrt(r(i)/b)*(x(j)-r(i))/b)
         enddo
      enddo
c
c      add the diagonal (wrt radial coord) coulomb potential and for XE the bardsley pseudo-potential
c      ad = 4.044d0
c      aq = 14.235d0
c      d = 1.d0
      do i = 1,nr
          hr(i,i) = hr(i,i)-1.d0/r(i)
c     +                      -ad*0.5d0/((r(i)**2+d**2)**2)
c     +                      -aq*0.5d0/((r(i)**2+d**2)**3)
      enddo
c     add the diagonal (wrt angular and radial coord) centrifugal potential
      k = 1
      do i = 1,nr
         do j = 1,nt
            h(k,k) = cent(j)/r(i)**2 !+pseudopot(j-1,r(i)) !l-dependent pseudo-potential for xe
            k = k + 1
         enddo
      enddo
c     angular potential similarity transform (angular DVR => FBR)
      r0 = 0.d0
      r1 = 1.d0
      do i = 1,nr
         do j = 1,nt
            do k = 1,nt
               ht(j,k,i) = 0.d0
            enddo
            zz = r(i)*cost(j)
            zd = zz + dist0
            if (zd.lt.0.d0) then
               ht(j,j,i) = potsurf(r(i),zz,dist0)
            else
               ht(j,j,i) = ef*zd+potsurf(r(i),zz,dist0)
            endif
         enddo
         call dgemm ('n','t',nt,nt,nt,r1,ht(1,1,i),nt,btg,nt,r0,vt,nt)
         call dgemm ('n','n',nt,nt,nt,r1,btg,nt,vt,nt,r0,ht(1,1,i),nt)
      enddo
c     Build total Hamiltonian
      do i = 1,n
         ir = (i-1)/nt + 1
         it = mod(i-1,nt) + 1
         do j = 1,n
            jr = (j-1)/nt + 1
            jt = mod(j-1,nt) + 1
            if (i.ne.j) h(i,j) = 0.d0
            if (it.eq.jt) then
                h(i,j) = h(i,j) + hr(ir,jr)
            endif
            if (ir.eq.jr) then
               h(i,j) = h(i,j) + ht(it,jt,ir)
            endif
         enddo
      enddo
c     Diagonalise Hamiltonian
c     print*,'diagonalising intial Hamiltonian'
      ierr = 0
      call symevp (h,n,e,ierr)
      if (ierr.ne.0) stop 'hdvr | dgeev 2'
      if (ntest.eq.1) then
         enn = -0.5d0/dble(ni)**2 + 0.25d0/dist0
         do i = 1,n
            ei = e(i)!*4.35974417d-18/1.60218d-19
c            if ((ei.ge.-0.34d0).and.(ei.le.-0.18d0)) then
            diff = dabs((enn-ei)/enn)
            if (diff.le.0.2d0) then
               ndd = int(dist0)
               write(ndd,*) ei,i,diff,enn
               print*,ei,i,diff,enn
            endif
         enddo
         stop
      endif
c     store selection of wavefunctions (eigenvectors)
      do istat = 1,mm
         nwav = n0(istat)
         do j = 1,n 
            hwav(j,istat) = h(j,nwav)
         enddo
c    work out energy shift from infinite distance
         denn = (e(nwav)-ef*dist0-0.25d0/dist0)
     &     -en(dble(ni),0.d0,dble(k0(istat)),ef)
c    work out corresponding velocity change from infinite distance
         velacc(istat) = -dsqrt(vel0**2-denn*2.d0/1836.d0)
      enddo
c      
c     Tranform Eigenvector to CWDVR basis:
c     INTERPOLATE RADIAL POINTS
      do k = 1,nt        !loop over angles
         do j = 1,nx     ! loop over interpolation points (cwdvr)
            jc = nt*(j-1)+k
            phi(j,k) = 0.0d0
            do i = 1,nr  !loop over basis set (laguerre)
               ic = nt*(i-1)+k
               phi(j,k) = phi(j,k)
     +                    +(hwav(ic,1)*basis(j,i)/dsqrt(b))/a(j)
c   coeff for state ij with angle k and r=j on CWDVR is: 
c   (sum of coeff for state ij over the laguerre basis i=1,nr) x CWDVR weight for r=j 
            enddo
         enddo
      enddo
c     interpolate angular functions
      do i = 1,nt
         do j = 1,ntcw
              ptg(i,j) = pbtg(i,j)
         enddo
      enddo
c     coeff of phi (in FBR) x [FBR basis functions(evaluated at ntcw angular points)*gaussian weight]
c     ---> angular DVR with ntcw ang points
      r0 = 0.d0
      r1 = 1.d0
      call dgemm('n','n',nx,ntcw,nt,r1,phi,nx,ptg,nt,r0,p,nx)
c
c normalise wavefunction
      aintegral = 0.d0
      do i = 1,nx 
         do j = 1,ntcw       
            aintegral = aintegral+p(i,j)**2
         enddo
      enddo
      sa = dsqrt(aintegral)
      do i = 1,nx 
         do j = 1,ntcw
            p(i,j)=p(i,j)/sa 
         enddo
      enddo
c
      do i = 1,nx
         do j = 1,ntcw
            psi(i,j) = dcmplx(p(i,j),0.d0)
         enddo
      enddo
      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wavfunc(pbtg,psi,hwav,b,nx,ntcw,a,nr,basis,nt,istat,mm)
      implicit double precision (a-h,o-z)
c--------------------------------------------------------------------------
c     NB x(nx) are the cwdvr grid points and r(nr) are the laguerre points
c     cost(nt) are the ang points for initial diagonalisation and cost(ntcw) and those for propagation
c--------------------------------------------------------------------------
      dimension phi(nx,nt),hwav(nr*nt,mm),basis(nx,nr)
      dimension a(nx),ptg(nt,ntcw),pbtg(ntcw,ntcw)
      dimension p(nx,ntcw)
      double complex psi(nx,ntcw)
c
      n = nr*nt
      do k = 1,nt        !loop over angles
         do j = 1,nx     ! loop over interpolation points (cwdvr)
            jc = nt*(j-1)+k
            phi(j,k) = 0.0d0
            do i = 1,nr  !loop over basis set (laguerre)
               ic = nt*(i-1)+k
               phi(j,k) = phi(j,k)
     +                    +(hwav(ic,istat)*basis(j,i)/dsqrt(b))/a(j)
c   coeff for state ij with angle k and r=j on CWDVR is: 
c   (sum of coeff for state ij over the laguerre basis i=1,nr) x CWDVR weight for r=j 
            enddo
         enddo
      enddo
c     interpolate angular functions
      do i = 1,nt
         do j = 1,ntcw
              ptg(i,j) = pbtg(i,j)
         enddo
      enddo
c     coeff of phi (in FBR) x [FBR basis functions(evaluated at ntcw angular points)*gaussian weight]
c     ---> angular DVR with ntcw ang points
      r0 = 0.d0
      r1 = 1.d0
      call dgemm('n','n',nx,ntcw,nt,r1,phi,nx,ptg,nt,r0,p,nx)
c
c normalise wavefunction
      aintegral = 0.d0
      do i = 1,nx 
         do j = 1,ntcw       
            aintegral = aintegral+p(i,j)**2
         enddo
      enddo
      sa = dsqrt(aintegral)
      do i = 1,nx 
         do j = 1,ntcw
            p(i,j)=p(i,j)/sa 
         enddo
      enddo
c
      do i = 1,nx
         do j = 1,ntcw
            psi(i,j) = dcmplx(p(i,j),0.d0)
         enddo
      enddo
      return 
      end

      subroutine lagdvr(n,b,r,t)
      implicit double precision (a-h,o-z)
      dimension r(n),d(n),t(n,n),z(n,n),work(2*n-2)!,wt(n)
c     -------------------------------------------------------------------------
c     compute the zeroes of Ln(x) to obtain n DVR grid points
c     ans also compute the kinetic energy operator matrix for reg-laguerre DVR
c     -------------------------------------------------------------------------

      do i = 1,n
         r(i) = (2*i-1)
         d(i) = -i
      enddo
c     compute zeros (eigenvalues)
      ierr = 0
      call dstev('V',n,r,d,z,n,work,ierr)
      if (ierr.ne.0) stop 'lagrt | dgtsv'
c
cc      do i=1,n
cc         rtgwt=z(1,i)
cc         wt(i)=rtgwt*sqrt(b)/(dexp(-r(i)/(2.d0)))
cc      enddo
c     computer kinetic energy operator matrix
      obs = 1.d0/(b*b)
      do i = 1,n
         do j = 1,n
            if (i.ne.j) then
                if (mod(i-j,2).eq.0) then
                   t(i,j) = ((r(i)+r(j))/((dsqrt(r(j)*r(i)))*
     +                      ((r(i)-r(j))**2)))*obs *0.5d0
                else
                   t(i,j)=((r(i)+r(j))/((dsqrt(r(j)*r(i)))*
     +                    ((r(i)-r(j))**2)))*obs *(-0.5d0)

                endif
            else
                t(i,i)=((4.d0+(4.d0*n+2.d0)*r(i)-r(i)**2)
     +                 /(12.d0*r(i)**2))*obs *0.5d0
            endif
         enddo
      enddo
      do i = 1,n
         r(i) = r(i)*b
      enddo
      write(6,'(A,i9.2)')'initial diagonalisation nr = ',n
      write(6,'(A,f9.2)')'initial diagonalisation grid rmax = ',r(n)
      return
      end

      subroutine lagtrans(x,nx,trans,nr,b)
      implicit double precision (a-h,o-z)
      dimension x(nx),ts(0:nr),trans(nx)
c     ---------------------------------------------------------------------------
c     computer value of Baye reg-laguerre basis functions at CWDVR grid points
c     ---------------------------------------------------------------------------
      do i = 1,nx
         ts(0) = 1.d0
         ts(1) = 1.d0-x(i)/b
         do j = 1,nr-1
            xj = dble(j)
            ts(j+1) = (2.d0*xj+1.d0-x(i)/b)*ts(j)/
     +                   (xj+1.d0)-(xj)*ts(j-1)/(xj+1.d0)
         enddo
         trans(i) = (x(i)/b)*ts(nr)*dexp(-x(i)*0.5d0/b)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine flux2(psi,cost,mrflag,voptic,fg2,fl2,vel)
      implicit double precision (a-h,o-z)
      double complex psi
      common /array/ nr,nt
      dimension psi(nr,nt),cost(nt)
      dimension voptic(nr)
      fg2 = 0.d0
      fl2 = 0.d0
      do k = 1,nt
         dp = 0.d0
         do i = mrflag,nr
            if (voptic(i).ne.0.d0) then
               dp = dp + 
     &          2.d0*dble(dconjg(psi(i,k))*voptic(i)*psi(i,k))
     &         /vel
            endif
         enddo
         if (cost(k).gt.0.d0) then
            fg2 = fg2+dp
         else
            fl2 = fl2+dp
         endif
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine anginterp(npt,ptg,pcost)
      implicit double precision (a-h,o-z)
c     --------------------------------------------------------------------------
c     work out angular FBR function at npt angular plotting points
c     --------------------------------------------------------------------------
      common /array/ nr,nt
      dimension pbtg(npt,npt),pcost(npt),pcent(npt),ptwt(npt)
      dimension ptg(nt,npt)
      call legdvrw(npt,pcost,pbtg,ptwt,pcent)
      do i = 1,nt
         do j = 1,npt
              ptg(i,j) = pbtg(i,j)/sqrt(ptwt(j))
         enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

