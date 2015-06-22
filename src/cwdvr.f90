      program cwdvr_wavepacket
!     calculates the time evolution of H atom Rydberg state
!     electronic wavefunction as it approaches the metal surface       
!     The wf is propagated on Coulomb Wave DVR (radial) and Legendre
!     DVR (angular)
!     Several options are possible:
!     1) calculate using a constant velocity for the incidence
!     2) calculate using the mean-field approximation: the positive ion
!        travels under a mean-potential described by the electronic
!        wavefunction
!     3) can described incidence at Jellium aluminium surface, or at
!        Cu(100) or Cu(111) surface (using Chulkov pseudopotential) 
!        just comment or uncomment lines in potsurf subroutine
!     4) can easily uncomment out pseudo-potential for xenon calcs
!        (search pseudo and the terms should appear)
!     see Eric So thesis
      use settings, only : loadSettings
      implicit none
      call loadSettings()
      call main()
      
      stop
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine main()
      use settings, only : v0, MEANFIELD, MPLOT, n0, nStates, field, dist0, b, min_pop, &
              timestep, nt, nr, ni, dk, npt, z, rmid, npstep, nfstep, nwstep, nrb, ntb, dp, TEST, TESTWF
      use helpers, only : pnorm
      implicit none
      integer, parameter :: nmax = 2000  ! the upper limit of the number of CWDVR radial points
      real(dp) :: r1(nmax)
      
      real(dp), allocatable :: r(:),a(:),wt(:)
      real(dp) :: velacc(nStates)
      real(dp) :: btg(nt,nt),cost(nt)
      real(dp) :: ptg(nt,npt),pcost(npt)
      real(dp) :: hwav(nrb*ntb,nStates)
      real(dp), dimension(:, :), allocatable :: basis
      complex(dp), allocatable :: psi(:, :),etr(:, :, :),ev2(:, :)
      character(LEN=17) filen1
      character(LEN=18) filen2
      character(LEN=15) filen3
      character(LEN=8) filen4
      character(LEN=13) filen6
      character(LEN=17) filen5
      real(dp), allocatable :: zx(:, :), voptic(:)
      integer :: t0, t1, clkrate
      real(dp) :: elapt
      integer :: i, j, mrflag, istat, nw, istep, ionstep, writstep, plotstep, fstep
      real(dp) :: vel, ptold, dist, dedz, dft, p0, time, pt
      real(dp) :: fg2, fl2, dedzo
      
      
!     pre stuff that used not to be in main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     get CWDVR zeros roughly, by scanning r and finding when the sign of the coulomb function changes
      call couzero(r1,nr,nmax)
!     we now know the number of zero crossings, and can allocate arrays
      allocate(r(nr))
      allocate(a(nr))
      allocate(wt(nr))
      allocate(basis(nr,nrb))
      allocate(psi(nr,nt))
      allocate(etr(nr,nr,nt))
      allocate(ev2(nr,nt))
      allocate(zx(nr,nt))
      allocate(voptic(nr))

      if (MEANFIELD .eqv. .TRUE.) then
         print*,'Calculation with Mean-Field Approx.'
      else
         print*,'Calculation with Constant Velocity Approx'
      endif
      print*,'---------------------------------------'
      print*,'Parameters:'
      print*,'====================='
      write(6,'(A,i5.2,/)')'hydrogen n = ',ni
      write(6,'(A,/,A,/,A,i5.2,/,A,f5.2,/,A,i5.2,/)')'initial diagonalisation with lag DVR:',&
                                        '-------------------------------------',&
                                        'nr = ',nrb,'b = ',b,'nt = ',ntb

      write(6,'(A,/,A,/,A,f5.2,/,A,f5.2,/)')'cwdvr grid parameters:','----------------------','dk = ',dk,'z = ',z

      write(6,'(A,/,A,/,A,i5.2,/,A,i5.2,/,A,f5.2,/,A,i5.2,/,A,f5.2,/,A,f5.2/)') & 
                                        'wpp:','----','number of cwdvr pts = ',nr,'angular points nt = ', & 
                                        nt,'dt = ',timestep, 'time steps between outputs=',npstep, &
                                        'velocity at infinite distance= ',v0,'E-field=',field
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     setup
      
      call system_clock(t0, clkrate)           !clear timer
      do i = 1,nr             !fill r with the approximate zeroes calculated previously
         r(i) = r1(i)
      enddo
!     Use Newton-Raphson method (differential required) to find exact zeroes and radial weights
      call rtnewt(r,a,wt)
! 
!     work out operators, exponentiate Kinetic operator, work out flux operator
      call xoperats(r,etr,cost,btg,mrflag,voptic)
!
!     work out angular FBR function at npt plotting points

      if (MPLOT .eqv. .TRUE.) then
         call anginterp(npt,ptg,pcost)
      endif
!
!     work out electron z position for all r , cos(theta)
      do i = 1,nr
         do j = 1,nt
            zx(i,j) = r(i)*cost(j)
         enddo
      enddo
!
!     set over states to propagate
      do istat = 1,nStates
         nw = n0(istat)
!        get intial wavefunction
         if (istat.eq.1) then
            write(6, *)'calling init'
            call init(btg,psi,hwav,r,nr,nt,a,nrb,basis,ntb,velacc) !diagonalise in regularised-laguerre DVR
!           hwav contains all the wanted mm eigenvectors from initial diagonalisation            
         else
!           no need to digaonlise again just project next state in hwav onto CWDVR and propagate 
            call wavfunc(btg,psi,hwav,nr,nt,a,nrb,basis,ntb,istat,nStates)
         endif
!         
!     open files for output
         write(filen1,'(A,i2.2)')'forward_flux.',nw
         write(filen2,'(A,i2.2)')'backward_flux.',nw
         write(filen4,'(A,i2.2)')'pop.',nw
         open(100+nw,file=filen1,status='unknown')
         open(200+nw,file=filen2,status='unknown')
         open(400+nw,file=filen4,status='unknown')
         if (MPLOT .eqv. .TRUE.) then
            write(filen5,'(A,i2.2)')'wavefunction.',nw
            open(500+nw,file=filen5,status='unknown')
            write(500+nw,*) nr,npt,timestep*npstep,rmid
         endif
         if (MEANFIELD .eqv. .TRUE.) then
            write(filen3,'(A,i2.2)')'totalforce.',nw
            write(filen6,'(A,i2.2)')'velocity.',nw
            open(300+nw,file=filen3,status='unknown')
            open(600+nw,file=filen6,status='unknown')
         endif
!     
         if (istat.eq.1) then
            call system_clock(t1, clkrate)
            write (6,61) (t1-t0)/clkrate/60.d0
  61        format(/1x,'Setup completed in:',f9.2,' minutes'//1x) 
            print*,'Time....................Population...................Dist'
         endif
!         
!        propagation
         if (MEANFIELD .eqv. .TRUE.) then
            vel = +velacc(istat) !include energy shift acceleration
         else
            vel = v0 !assume velocity at infinite distance
         endif
         t0 = t1
         ptold = 10.d0
         istep = -1
         dist = dist0
         ionstep = 0
         dedz = 0.d0
         dft = nfstep*timestep
         p0 = 1.d0
1        istep = istep+1
         time = istep*timestep
         if (istep .gt. 0) then
!        work out time-depedent potential energy operator
            call potop (ev2,r,dist,zx)
!        propagate wavefunction with split operator
            call split (psi,etr,ev2,btg)
        endif
!
!        work out population on grid at time t
         pt = pnorm(psi)

         
!
!        output
         writstep = mod(istep,nwstep)
         if (writstep.eq.0.d0) then
!           output flux
            call flux2(psi,cost,mrflag,voptic,fg2,fl2,vel)
            write(100+nw,*)dist,fl2
            write(200+nw,*)dist,fg2
!           work out wavefunction population inside area bound by absorbing boundary and surface at time t
!            pt2 = pnorm2(psi,dist,r,cost,rmid)
            write(400+nw,*) dist,pt
            write (6,*) time,pt,dist
         endif
!        plot wavefunction
         if (MPLOT .eqv. .TRUE.) then
            plotstep = mod(istep,npstep)
            if (plotstep.eq.0.d0) then
               call plot(psi,r,wt,nw,npt,btg,ptg,pcost,dist)
            endif
         endif
!
          if (MEANFIELD .eqv. .FALSE.) then
             dist = dist + vel*timestep
             goto 426
          endif
!        work out force on ion-core and propgate with velocity verlet
         fstep = mod(istep,nfstep)
         if (fstep.eq.0.d0) then
            call denergy2(psi,r,dist,zx,dedz)
            if (writstep.eq.0.d0) then
                write(600+nw,*) dist,vel,time
                write(300+nw,*) dist,dedz
            endif
            call fion (psi,dedz,dft,vel,zx,ionstep,dedzo)
            ionstep = 1
         endif
!
         dist = dist + vel*timestep -0.5d0*dedz*timestep**2
         
!         
!        Stop if ion goes backwards too far!
426      if (dist > (dist0 + 50.d0)) stop 'reversed too far!'
!
!        WARN if wavepacket grows
         if (pt - ptold > 1.d-8) print*, 'growing!'
         ptold = pt
!
!        stop when distance from surface =0
         if (dist < 0.d0) go to 847
         if ((pt > min_pop) .AND. .NOT.TESTWF)  go to 1
847      call system_clock(t1, clkrate)
         write (6,61) (t1-t0)/clkrate/60.d0
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plot(psi,r,wt,nw,npt,btg,ptg,pcost,dist)
      use settings, only : nr, nt, dp
      implicit none
      integer :: nw, npt, nrs, i, j
      real(dp) :: dist
      complex(dp) psi(nr,nt),phi(nr,nt),pxi(nr,npt)
      real(dp) :: r(nr),wt(nr),btg(nt,nt)
      real(dp) :: pcost(npt),ptg(nt,npt)
      real(dp) :: r0, r1, si, zz, ro, dble
!     coeffs: angular DVR-> angular FBR conversion
      r0 = 0.d0
      r1 = 1.d0
      nrs = 2*nr
      call dgemm('n','t',nrs,nt,nt,r1,psi,nrs,btg,nt,r0,phi,nrs)
!     coeff x FBR basis functions (evaluated at npt angular points)
      r0 = 0.d0
      r1 = 1.d0
      nrs = 2*nr
      call dgemm('n','n',nrs,npt,nt,r1,phi,nrs,ptg,nt,r0,pxi,nrs)
      do i = 1,nr
         do j = 1,npt
            si = dble((pxi(i,j))*conjg(pxi(i,j)))/(wt(i)) !*pi/(wt(i))
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine expmih (etr,tr,r,cent,voptic)
      use settings, only : timestep, nr, nt, rmid, rabs, eminr, dp
      use helpers, only : expmat
      implicit none
!     -----------------------------------------------------------------
!     Sets up the exponentiated Hamiltonian factors
!     needed for split operator propagation.
!     -----------------------------------------------------------------
      complex(dp) :: etr(nr, nr, nt), vopt
      real(dp) :: tr(nr,nr),r(nr)
      real(dp) :: cent(nt)
      real(dp) :: voptic(nr)
      integer :: i, j, k
      real(dp) :: rmrmid

!     exp(-iTdt)
      do k = 1,nt
         do j = 1,nr
            voptic(j) = 0.d0       
            do i = 1,nr
               etr(i,j,k) = cmplx(0.d0, -timestep, kind=dp)*tr(i,j)
            enddo
            rmrmid = r(j)-rmid
            call optpot (rmrmid,rabs,eminr,vopt)
            if (rmrmid .gt. 0.d0) then
                voptic(j) = aimag(vopt)
            endif
            vopt = tr(j,j)+(cent(k)/r(j)**2)+vopt!+pseudopot(k-1,r(j))
            etr(j,j,k) = cmplx(0.d0, -timestep, kind=dp)*vopt
         enddo
         call expmat(etr(1,1,k),nr)
      enddo
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine optpot (rmrmid,rabs,eminr,vopt)
      use settings, only : dp
      implicit none
!     -----------------------------------------------------------------
!     Transmission-free negative imaginary absorbing potential.
!     -----------------------------------------------------------------

      real(dp) :: rmrmid, rabs, eminr, x, y
      complex(dp) :: vopt
      real(dp), parameter :: c = 2.62206d0
      vopt = (0.d0,0.d0)
      if (rmrmid .gt. 0.d0) then
         x = c*rmrmid/rabs
         y = 4.d0/(c-x)**2+4.d0/(c+x)**2-8.d0/c**2  
         vopt = cmplx(0.d0,-eminr*y, kind=dp)
      endif
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xoperats(r,etr,cost,btg,mrflag,voptic)
      use settings, only : field, dk, z, nr, nt, TEST, rmid, dp
      implicit none
      real(dp) :: r(nr),tr(nr,nr)
      real(dp) :: btg(nt,nt),cost(nt),cent(nt)
      complex(dp) etr(nr,nr,nt)
      real(dp) :: voptic(nr)
      integer :: mrflag, ii
! 
!     radial dvr kinetic operator
!
      call coutr(r,nr,dk,z,tr)
!
!     angular grid points,transform matrix,weights,centrifugal
!
      call legdvr(nt,cost,btg,cent)
!
      if (TEST .eqv. .FALSE.) then
!        exponentiate kinetic energy operator
         print*,'exponentiating KE operator'
         call expmih(etr,tr,r,cent,voptic)
!         write(6,'(f10.5f10.5/)')REAL(etr(1, 2, 2)),AIMAG(etr(1, 2, 2))
!   
!        work out absorbing boundary limit
         do ii = 1,nr
            if (r(ii).ge.rmid) then
               mrflag = ii
               goto 555
            endif
         enddo
      endif
555   return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine couzero(r,nr,nmax)
      use settings, only : eta1, xmin, nx, xlmin, xlmax, dx, dp
      implicit none
!     -----------------------------------------------------------------
!     approx location of coulomb function, F_0(eta,x), zeros for DVR grid.
!     -----------------------------------------------------------------
      integer :: nr, ifail, nmax, i
      real(dp) :: r(nmax)
      real(dp) :: x, f(1), g(1), gp(1), fp(1), fo
!
!     First find approximate location of zeros by stepping
      nr = 0
      ifail = 0
!     called wth: xmin=8e-3,eta1=-50/3,xlmin=0,xlmax=0
      call coulfg(xmin,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
      if (ifail.ne.0) stop 'couzero | coulfg 1'


      fo = f(1)
!
      do i = 1,nx+1
         x = xmin+i*dx 
         ifail = 0
         call coulfg(x,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
         if (ifail.ne.0) stop 'couzero | coulfg 2'
         if (f(1)*fo.lt.0.d0) then
            nr = nr+ 1   
            if (nr.gt.nmax) stop 'couzero | nr > nmax'
            r(nr) = x - 0.5d0*dx
         endif
         fo = f(1)
      enddo
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rtnewt(r,a,wt)
      use settings, only : field, dk, zacc, nr, xlmin, xlmax, dx, eta1, pi, dp
      implicit none
      real(dp) :: r(nr),a(nr),wt(nr)
      
      integer :: i, ncount, IFAIL, J
      real(dp) :: DABS, dif, xmaxdif, slimit, abs, dz, sqrt
      real(dp) :: f(1), g(1), fp(1), gp(1), x1, x2, xz

!     Newton-Raphson method to find roots more accurately. 
      do i = 1,nr
         x1 = r(i) - 0.5d0*dx 
         x2 = r(i) + 0.5d0*dx 
         xz = r(i)
         j=0
123      j=j+1
         ifail = 0
         call coulfg(xz,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
         if (ifail.ne.0) stop 'couzero | coulfg 3'
         dz = f(1)/fp(1)
         xz = xz - dz
         if ((x1-xz)*(xz-x2).lt.0.d0) then
           write(6, *) 'jumped out of bracket'
           stop
         endif
         if (abs(dz).gt.zacc) go to 123
!         print*, 'newtr converged after',j,'iterations'
         call coulfg(xz,eta1,xlmin,xlmax,f,g,fp,gp,1,0,ifail)
         if (ifail.ne.0) stop 'couzero | coulfg 3'
         a(i) = fp(1)*sqrt(dk)
         wt(i) = pi/(dk*(fp(1)**2))
         r(i) = xz/dk
         !write(6,'(f8.2/)')xz
      enddo
!     test if seperation of grid points is ok
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
!         pause 'need to inc. dk paramter!!!!!!!!'
      endif
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE COULFG(XX,ETA1,XLMIN,XLMAX, FC,GC,FCP,GCP,MODE1,KFN,IFAIL)
      implicit none
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      !
!  REVISED COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           !
!                                                                      !
!  A. R. BARNETT           MANCHESTER  MARCH   1981                    !
!                                                                      !
!  ORIGINAL PROGRAM 'RCWFN'      IN    CPC  8 (1974) 377-395           !
!                 + 'RCWFF'      IN    CPC 11 (1976) 141-142           !
!  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314           !
!  THIS VERSION WRITTEN UP       IN    CPC XX (1982) YYY-ZZZ           !
!                                                                      !
!  COULFG RETURNS F,G,F',G', FOR REAL XX.GT.0,REAL ETA1 (INCLUDING 0), !
!   AND REAL LAMBDA(XLMIN) .GT. -1 FOR INTEGER-SPACED LAMBDA VALUES    !
!   THUS GIVING POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER   !
!   EQUATION,TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF     !
!   THE DIRAC EQUATION ,ALSO SPHERICAL + CYLINDRICAL BESSEL EQUATIONS  !
!                                                                      !
!  FOR A RANGE OF LAMBDA VALUES (XLMAX - XLMIN) MUST BE AN INTEGER,    !
!  STARTING ARRAY ELEMENT IS M1 = MAX0(  INT(XLMIN+ACCUR),0) + 1       !
!      SEE TEXT FOR MODIFICATIONS FOR INTEGER L-VALUES                 !
!                                                                      !
!  IF 'MODE' = 1  GET F,G,F',G'   FOR INTEGER-SPACED LAMBDA VALUES     !
!            = 2      F,G      UNUSED ARRAYS MUST BE DIMENSIONED IN    !
!            = 3      F               CALL TO AT LEAST LENGTH (1)      !
!  IF 'KFN'  = 0 REAL        COULOMB FUNCTIONS ARE RETURNED            !
!            = 1 SPHERICAL   BESSEL                                    !
!            = 2 CYLINDRICAL BESSEL                                    !
!  THE USE OF 'MODE' AND 'KFN' IS INDEPENDENT                          !
!                                                                      !
!  PRECISION&  RESULTS TO WITHIN 2-3 DECIMALS OF 'MACHINE ACCURACY'    !
!   IN OSCILLATING REGION X .GE. ETA1 + SQRT(ETA1**2 + XLM(XLM+1))     !
!   COULFG IS CODED FOR REAL*8 ON IBM OR EQUIVALENT  ACCUR = 10**-16   !
!   USE AUTODBL + EXTENDED PRECISION ON HX COMPILER  ACCUR = 10**-33   !
!   FOR MANTISSAS OF 56 + 112 BITS. FOR SINGLE PRECISION CDC (48 BITS) !
!   REASSIGN DSQRT=SQRT ETC.  SEE TEXT FOR COMPLEX ARITHMETIC VERSION  !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
      double precision :: XX, ETA1, XLMIN, XLMAX
      double precision :: FC(1),GC(1),FCP(1),GCP(1)
      integer :: MODE1, KFN, IFAIL
      LOGICAL      ETANE0,XLTURN
      double precision :: PACCQ
      integer :: NFP, NPQ, IEXP, M1
      COMMON       /STEED/ PACCQ,NFP,NPQ,IEXP,M1
!***  COMMON BLOCK IS FOR INFORMATION ONLY.  NOT REQUIRED IN CODE
!***  COULFG HAS CALLS TO&  SQRT, ABS, AMOD ,  INT, SIGN, FLOAT, MIN1 
      double precision :: ZERO, ONE, TWO, TEN2, ABORT
      DATA ZERO,ONE,TWO,TEN2,ABORT /0.0d0, 1.0d0, 2.0d0, 1.0d2, 2.0d4/
      double precision :: HALF, TM30
      DATA HALF,TM30 / 0.5d0, 1.0d-30 / 
!      DATA RT2EPI /0.79788 45608 02865/
! *** THIS CONSTANT IS  DSQRT(TWO/PI)&  USE Q0 FOR IBM REAL*16& D0 FOR
! ***  REAL*8 + CDC DOUBLE P&  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
!
!                       ACCUR = R1MACH(4)
      double precision :: pi, accur, RT2EPI, ETA, ACC, ACC4, ACCH, X, XLM, E2MM1, DELL, MOD, ABS, XLL, FLOAT, MAX0, XI
      double precision :: FCL, PK, PX, EK, F, PK1, D, FPL, XL, RL, EL, SL, A, B, C, DP, DQ, DI, DR, BI, BR, AI, AR, Q, WI, GCL1
      double precision :: GPL, GCL, W,  FCM, BETA, ALPHA, P, GAM, TA, FJWKB, FCL1, DF, TK, GJWKB
      integer :: MODE, LXTRA, INT, L1, LP, L, MAXL, SIGN
      pi=dacos(-1.d0)
      accur = 1.d-20
      RT2EPI=dsqrt(2.d0/pi)
! ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
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
! ***    TEST RANGE OF XX, EXIT IF.LE. SQRT(ACCUR) OR IF NEGATIVE
!
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
! ***       LXTRA IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
! ***       XLL  IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
! ***         DETERMINE STARTING ARRAY ELEMENT (M1) FROM XLMIN
      M1  = MAX0(  INT(XLMIN + ACC),0) + 1
      L1  = M1 + LXTRA
!
! ***    EVALUATE CF1  =  F   =  FPRIME(XL,ETA,X)/F(XL,ETA,X)
!
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
    2 EK  = ETA / PK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 =  PK + ONE
! ***   TEST ENSURES B1 .NE. ZERO FOR NEGATIVE ETA; FIXUP IS EXACT.
             IF( ABS(ETA*X + PK*PK1) .GT. ACC)  GO TO 3
             FCL  = (ONE + EK*EK)/(ONE + (ETA/PK1)**2)
             PK   =  TWO + PK 
      GO TO 2
    3 D   =  ONE/((PK + PK1)*(XI + EK/PK1))
      DF  = -FCL*(ONE + EK*EK)*D
            IF(FCL .NE. ONE )  FCL = -ONE
            IF(D   .LT. ZERO)  FCL = -FCL
      F   =  F  + DF
!
! ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
!
!      write(6,'(f15.4,f15.5,f5.4,f5.2/)')PK, F, D, FCL
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
!
! *** DOWNWARD RECURRENCE TO LAMBDA = XLM. ARRAY GC,IF PRESENT,STORES RL
!
!      write(6,'(f15.4,f15.5,f5.4,f5.2/)')DF, F, NFP, EK
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
         XL = XL - ONE 
    6 end do
      IF(FCL .EQ. ZERO) FCL = ACC
      F  = FPL/FCL
! ***    NOW WE HAVE REACHED LAMBDA = XLMIN = XLM 
! ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
! ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM 
!
    7 IF( XLTURN ) write(6, '(A6)')'XLTURN'
      IF( XLTURN ) CALL JWKB(X,ETA,MAX(XLM,ZERO),FJWKB,GJWKB,IEXP)
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
!      write(6,'(f15.5,f15.5,f15.5,f15.5/)')DR, DI, X, ETA
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
!
! *** SOLVE FOR FCM = F AT LAMBDA = XLM,THEN FIND NORM FACTOR W=W/FCM 
!
      GAM = (F - P)/Q
            IF(Q .LE. ACC4* ABS(P))             GO TO 130
      W   = ONE/ dSQRT((F - P)*GAM + Q)
            GO TO 10
! *** ARRIVE HERE IF G(XLM) .GT. 10**6 OR IEXP .GT. 250 + XLTURN = .TRUE.
    9 W   = FJWKB
      GAM = GJWKB*W 
      P   = F
      Q   = ONE
!
! *** NORMALISE FOR SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
!
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
! *** UPWARD RECURRENCE FROM GC(M1),GCP(M1)  STORED VALUE IS RL
! *** RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
! ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
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
         FC(L+1)     = W* FC(L+1)
  12  end do
      RETURN
 1000 FORMAT(/' CF1 ACCURACY LOSS& D,DF,ACCH,K,ETA/K,ETA,X = ',1P,7E9.2/)
!
! ***    ERROR MESSAGES
!
  100 IFAIL = -1
      WRITE(6,2000) XX,ACCH
 2000 FORMAT(' FOR XX = ',1P,E12.3,' TRY SMALL-X  SOLUTIONS',' OR X NEGATIVE',' SQUARE ROOT ACCURACY PARAMETER =  ',E12.3/)
      RETURN
  105 IFAIL = -2
      WRITE (6,2005) XLMAX,XLMIN,XLM
 2005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES&XLMAX,XLMIN,XLM = ',1P,3E15.6/)
      RETURN
  110 IFAIL =  1
      WRITE (6,2010) ABORT,F ,DF,PK,PX,ACC
 2010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',F10.0,' ITERATIONS',' F,DF,PK,PX,ACCUR =  ',1P,5E12.3//)
      RETURN
  120 IFAIL =  2
      WRITE (6,2020) ABORT,P,Q,DP,DQ,ACC
 2020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',F7.0,' ITERATIONS',' P,Q,DP,DQ,ACCUR =  ',1P,4E17.7,E12.3//)
      RETURN
  130 IFAIL =  3
      WRITE (6,2030) P,Q,ACC,DELL,LXTRA,M1
 2030 FORMAT(' FINAL Q.LE. ABS(P)*ACC*10**4 , P,Q,ACC = ',1P,3E12.3,4X,' DELL,LXTRA,M1 = ',E12.3,2I5 /)
      RETURN
 2040 FORMAT(' XLMAX - XLMIN = DELL NOT AN INTEGER ',1P,3E20.10/)
      END 
      
      
      SUBROUTINE JWKB(XX,ETA1,XL,FJWKB,GJWKB,IEXP)
! *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS    FOR XL.GE. 0
! *** AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
! *** CALLS AMAX1,SQRT,ALOG,EXP,ATAN2,FLOAT,INT        BARNETT FEB 1981

      implicit none
      double precision :: XX,ETA1,XL,FJWKB,GJWKB
      double precision :: X, ETA, GH2, SQRT, HLL, XLL1, GH, RL2, SL, PHI10, EXP, FLOAT, PHI, HL
      integer :: IEXP
      double precision :: ZERO, HALF, ONE, SIX, TEN
      double precision :: DZERO, RL35, ALOGE
      DATA   ZERO,HALF,ONE,SIX,TEN/ 0.0E0, 0.5E0, 1.0E0, 6.0E0, 10.0E0 /
      DATA  DZERO, RL35, ALOGE  /0.0E0, 35.0E0, 0.43429 /
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coutr(r,n,dk,z,t)
      use settings, only : dp
      implicit none
!     -----------------------------------------------------------------
!     Coulomb dvr kinetic energy
!     -----------------------------------------------------------------
      integer :: n
      real(dp) :: dk, z
      real(dp) :: r(n),t(n,n)
      real(dp) :: te, ri, rj
      integer :: i, j
      te = dk**2
      do i = 1,n
         ri = r(i)
         do j = 1,n
            rj = r(j) 
            t(i,j) = 1.d0/(ri-rj)**2 
         enddo 
         t(i,i) =(te + 2.d0*z/ri)/6.d0  ! why is this not negative?? (cf eq. 9 in paper!)
      enddo
!     add diag coulomb potential and for case of XE, bardsley pseudo-potential
!      ad = 4.044d0
!      aq = 14.235d0
!      d = 1.d0
      do i = 1,n
         t(i,i) = t(i,i)-1.d0/r(i) 
!     +                            -ad*0.5d0/((r(i)**2+d**2)**2)
!     +                           -aq*0.5d0/((r(i)**2+d**2)**3)
      enddo
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine symevp (a,n,d,ierr)
      use settings, only : dp
      implicit none
!     ----------------------------------------------------------------- 
!     Uses LAPACK DSYEV to diagonalise a real symmetric matrix.
!     ----------------------------------------------------------------- 
      
      integer :: n, ierr
      real(dp) :: a(n,n),d(n)
      real(dp) :: work(34*n)
      integer :: lwork, i, j
      real(dp) :: dsqrt, s
!
      lwork = 34*n 
      call dsyev ('V','U',n,a,n,d,work,lwork,ierr) 
!
!     normalise eigenvectors
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
     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine legdvrw(n,x,btg,twt,cent)
      use settings, only : dp
      implicit none
!     -----------------------------------------------------------------
!     Legendre polynomial discrete variable representation.
!     -----------------------------------------------------------------
      integer :: n
      real(dp) :: x(n),d(n),btg(n,n),work(2*(n-1)),twt(n),cent(n)
      integer :: i, ierr, j
      real(dp) :: ovnorm
!     compute jacobi matrix
      do i = 1,n
         x(i) = 0.d0
         d(i) = (dble(i)/(dsqrt(4.d0*(i**2)-1.d0)))
      enddo
!     compute zeros (eigenvalues)
      ierr = 0
      call dstev('V',n,x,d,btg,n,work,ierr)
      if (ierr.ne.0) stop 'lagrt | dgtsv'
!     work out legendre weights
      ovnorm = dsqrt(2.d0)
      do i = 1,n
         if (btg(1,i).lt.0.d0) then
          do j = 1,n
             btg(j,i) = btg(j,i)*(-1.d0)
          enddo
         endif
         twt(i) = (btg(1,i)*ovnorm)**2
      enddo
!     work out centrifugal term
      j = 0
      do i = 1,n
         cent(i) = 0.5d0*j*(j+1)
         j = j+1
      enddo
      return
      end
      
      
      subroutine legdvr(n,x,btg,cent)
      use settings, only : dp
      implicit none
!     -----------------------------------------------------------------
!     Legendre polynomial discrete variable representation.
!     -----------------------------------------------------------------
      integer :: n
      real(dp) :: x(n),d(n),btg(n,n),work(2*(n-1)),cent(n)
      integer :: i, ierr, j
      real(dp) :: ovnorm
!     compute jacobi matrix
      do i = 1,n
         x(i) = 0.d0
         d(i) = (dble(i)/(dsqrt(4.d0*(i**2)-1.d0)))
      enddo
!     compute zeros (eigenvalues)
      ierr = 0
      call dstev('V',n,x,d,btg,n,work,ierr)
      if (ierr.ne.0) stop 'lagrt | dgtsv'
!     work out legendre weights
      ovnorm = dsqrt(2.d0)
      do i = 1,n
         if (btg(1,i).lt.0.d0) then
          do j = 1,n
             btg(j,i) = btg(j,i)*(-1.d0)
          enddo
         endif
      enddo
!     work out centrifugal term
      j = 0
      do i = 1,n
         cent(i) = 0.5d0*j*(j+1)
         j = j+1
      enddo
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potop (ev2,r,dist,zx)
      use settings, only : field, timestep, nr, nt, dp
      use surface, only : potsurf
      implicit none
!     -------------------------------------------------------------------------
!     compute the time-dependent potential energy operator
!     -------------------------------------------------------------------------
      real(dp) :: dist
      complex(dp) :: ev2(nr,nt)
      real(dp) :: r(nr)
      real(dp) :: zx(nr,nt)
      integer :: i, j
      real(dp) :: dt2, zz, zd, v, ar
!     exp(-iVdt/2)
      dt2 = 0.5d0*timestep
      do i = 1,nr
         do j = 1,nt
            zz = zx(i,j)
            zd = zz + dist
            if (zd.lt.0.d0) then
               v = potsurf(r(i),zz,dist)
            else
               v = potsurf(r(i),zz,dist)+field*zd
            endif
            ar = -v*dt2
            ev2(i,j) = cmplx(dcos(ar), dsin(ar), kind=dp)
         enddo
      enddo
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine split (psi,etr,ev2,btg)
      use settings, only : nr, nt, dp
      implicit none
!     -----------------------------------------------------------------
!     Evolves the wavepacket through a time step dt
!     using the symmetric split operator method.
!     -----------------------------------------------------------------
!     common block
!     input arrays
      complex(dp) :: psi(nr*nt)
      complex(dp) :: etr(nr*nr*nt)
      complex(dp) :: ev2(nr*nt)
      real(dp) :: btg(nt*nt)
!     local array
      complex(dp) :: phi(nr*nt)
      complex(dp) :: c0, c1
      real(dp) :: r0, r1
      integer :: n, ij, k, ke, kp, nrs, ns
      
      
!     multiplication by exp(-iVdt/2)
      n = nr*nt
      do ij = 1,n
         psi(ij) = ev2(ij)*psi(ij)
      enddo
!     angular grid to basis transformation
      r0 = 0.d0
      r1 = 1.d0
      nrs = 2*nr
      call dgemm ('n','t',nrs,nt,nt,r1,psi,nrs,btg,nt,r0,phi,nrs)
!
!     multiplication by exp(-iTdt)
!
      c0 = (0.d0,0.d0)
      c1 = (1.d0,0.d0)
      ns = 1
!$OMP PARALLEL DO
!$OMP& SHARED(nt,nr,c0,c1,etr,phi,psi)
!$OMP& PRIVATE(k,ke,kp)
      do k = 1,nt
         ke = (k-1)*nr*nr + 1
         kp = (k-1)*nr + 1
         call zgemv('n',nr,nr,c1,etr(ke),nr,phi(kp),1,c0,psi(kp),1)
!         call zgemm
!     +   ('n','n',nr,ns,nr,c1,etr(ke),nr,phi(kp),nr,c0,psi(kp),nr)
      enddo
!$OMP END PARALLEL DO
!
!     angular basis to grid transformation
!
      nrs = 2*nr
      call dgemm ('n','n',nrs,nt,nt,r1,psi,nrs,btg,nt,r0,phi,nrs)
!
!     multiplication by exp(-iVdt/2)
!
      do ij = 1,n
         psi(ij) = ev2(ij)*phi(ij)
      enddo
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fion (psi,dedz,dft,vel,zx,ionstep,dedzo)
      use settings, only : mass, v0, nr, nt, dp
      implicit none
      real(dp) :: dedz, dft, vel, dedzo, vel0
      integer :: ionstep
      integer :: n
      complex(dp) :: psi(nr*nt)
      real(dp) :: zx(nr*nt)
      integer :: i
      real(dp) :: dsin, dcos, arg, v1, v5
      
      vel0 = vel
      if (ionstep.ne.0) then
        v1 = vel0 + 0.5d0*dedzo*dft/mass
      else
        v1 = vel0
      endif
      v5 = v1 + 0.5d0*dedz*dft/mass
      vel = v5
!
      n = nr*nt
      do i = 1,n
         arg = (-vel+vel0)*zx(i)
         psi(i) = psi(i)*cmplx(dcos(arg), dsin(arg), kind=dp)
      enddo
!
      dedzo = dedz
      return
      end
      
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine denergy2(psi,r,dist,zx,dedz)
      use settings, only : field, nr, nt, dp
      use surface, only : dpotsurf
      implicit none
      real(dp) :: dist, dedz
      complex(dp) :: psi(nr, nt)
      real(dp) :: r(nr)
      real(dp) :: zx(nr,nt)
      
      real(dp) :: vint, ev, zd, zz
      integer :: i, j
      dedz = 0.d0

      vint = 0.d0
      do i = 1,nt
         do j = 1,nr
            zz = zx(j,i)
            zd = zz + dist
            if (zd.lt.0.d0) then
               ev = dpotsurf(r(j),zz,dist)
            else
               ev = dpotsurf(r(j),zz,dist) + field
            endif
            vint =  vint + dble(conjg(psi(j,i))*ev*psi(j,i))
        enddo
      enddo
      dedz =  -(-field +0.25d0/(dist**2) + vint) 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function en(xn,xm,xk)
      use settings, only : field, dp
      implicit none
      real(dp) :: xn, xm, xk
      real(dp) :: a, b, c, d, e, g, h, p
      a=-1.d0/(2.d0*(xn**2))
      b=3.d0*field*xn*xk/2.d0
      c=-(field**2)*(xn**4)*(17.d0*(xn**2)-3.d0*(xk**2)-9.d0*(xm**2)+19.d0)/16.d0
      d=(3.d0/32.d0)*(xn**7)*xk*(23.d0*(xn**2)-(xk**2)+11.d0*(xm**2)+39.d0)*(field**3)
      e=-(xn**10)*(field**4)*(5487.d0*(xn**4)+35182.d0*(xn**2)-1134.d0*(xm**2)*(xk**2)+1806.d0*(xn**2)*(xk**2))/1024.d0
      g=-(xn**10)*(field**4)* & 
         (-3402.d0*(xn**2)*(xm**2)+147.d0*(xk**4)-549.d0*(xm**4)+5754.d0*(xk**2)-8622.d0*(xm**2)+16211.d0)/1024.d0
      h=3.d0*(xn**13)*xk*(field**5)*(10563.d0*(xn**4)+90708.d0*(xn**2)+220.d0*(xm**2)*(xk**2)+98.d0*(xn**2)*(xk**2))/1024.d0
      p=3.d0*(xn**13)*xk*(field**5)* &
        (772.d0*(xn**2)*(xm**2)-21.d0*(xk**4)+725.d0*(xm**4)+780.d0*(xk**2)+830.d0*(xm**2)+59293.d0)/1024.d0
      en=a+b+c+d+e+g+h+p
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function pseudopot(l,r)
      use settings, only : dp
      implicit none
!     -------------------------------------------------------------------------- 
!     Bardsley Pseudopotential for Xenon
!     ---------------------------------------------------------------------- 
      integer :: l
      real(dp) :: r, a, b
      if (l.eq.0) then
!       a=9.102d0
!       b=0.511d0
           a = -1.1798376d0
           b = 0.01d0
      endif
      if (l.eq.1) then
!       a=2.567d0
!       b=0.224d0
           a = -0.89700294d0
           b = 0.01d0
      endif
      if (l.eq.2) then
!       a=-0.468d0
!       b=0.190d0
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init(pbtg,psi,hwav,x,nx,ntcw,a,nr,basis,nt,velacc)
      use settings, only : v0, TEST, ni, n0, k0, nStates, field, dist0, mass, b, dp
      use surface, only : potsurf
      implicit none
!--------------------------------------------------------------------------
!     intitial diagonalisation on regularised Laguerre DVR to find
!     eigenstates to propagate      
!     NB x(nx) are the cwdvr grid points and r(nr) are the laguerre points
!     cost(nt) are the ang points for initial diagonalisation and cost(ntcw) and those for propagation
!--------------------------------------------------------------------------
      
      integer :: nx, ntcw, nr, nt
      real(dp) :: btg(nt,nt),cost(nt),cent(nt)
      real(dp) :: pbtg(ntcw,ntcw),ptg(nt,ntcw)
      real(dp) :: r(nr),hr(nr,nr)
      real(dp) :: e(nr*nt)
      real(dp), dimension(:, :), allocatable :: h
      real(dp) :: vt(nt,nt),ht(nt,nt,nr),phi(nx,nt)
      real(dp) :: p(nx,ntcw)
      real(dp) :: trans(nx),x(nx),a(nx),basis(nx,nr)
      real(dp) :: hwav(nr*nt,nStates)
      complex(dp) psi(nx,ntcw)
      real(dp) :: velacc(nStates)
      real(dp) :: aintegral, sa, denn, diff, enn, ei, zd, zz, r0, r1, en
      integer :: ic, jc, nwav, istat, ndd, ierr, jt, jr, it, ir, i, j, k, n
!
      allocate(h(nr*nt, nr*nt))
      n = nr*nt
!     work out nt legendre ang functions/points for initial diagonalisation
      call legdvr(nt,cost,btg,cent)
!
!     work out laguerre grid points and kinetic energy operator for regularised laguerre dvr
      call lagdvr(nr,b,r,hr)
!
!     work out L_(N)(x) with x at the nx cwdvr grid points
      call lagtrans(x,nx,trans,nr,b)
!
!     work out f_(i)(x) for i=1,...,nr laguerre basis functions at x(j) (j=1,...,nx) cwdvr points
!     BAYE REGULARISED LAGUERRE: f_(i)(x)=(-1)^(i)*x_(i)^(-0.5)*x*L_(N)(x)*e^(-x/2)/(x-x_(i))
      do j = 1,nx
         do i = 1,nr
            basis(j,i) = ((-1.d0)**i)*trans(j)/(dsqrt(r(i)/b)*(x(j)-r(i))/b)
         enddo
      enddo
      

      
!
!      add the diagonal (wrt radial coord) coulomb potential and for XE the bardsley pseudo-potential
!      ad = 4.044d0
!      aq = 14.235d0
!      d = 1.d0
      do i = 1,nr
          hr(i,i) = hr(i,i)-1.d0/r(i)
!     +                      -ad*0.5d0/((r(i)**2+d**2)**2)
!     +                      -aq*0.5d0/((r(i)**2+d**2)**3)
      enddo
      
!     add the diagonal (wrt angular and radial coord) centrifugal potential
      k = 1
      do i = 1,nr
         do j = 1,nt
            h(k,k) = cent(j)/r(i)**2 !+pseudopot(j-1,r(i)) !l-dependent pseudo-potential for xe
            k = k + 1
         enddo
      enddo
      

      
!     angular potential similarity transform (angular DVR => FBR)
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
               ht(j,j,i) = field*zd+potsurf(r(i),zz,dist0)
            endif
         enddo
      enddo
      do i = 1, nr
         call dgemm ('n','t',nt,nt,nt,r1,ht(1,1,i),nt,btg,nt,r0,vt,nt)
         call dgemm ('n','n',nt,nt,nt,r1,btg,nt,vt,nt,r0,ht(1,1,i),nt)
      enddo
      
      
!     Build total Hamiltonian
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
      

      
!     Diagonalise Hamiltonian
!     print*,'diagonalising intial Hamiltonian'
      ierr = 0
      call symevp (h,n,e,ierr)
      if (ierr.ne.0) stop 'hdvr | dgeev 2'
      if (TEST .eqv. .TRUE.) then
         print*,'picked k states'
         print*,'--------------'
         enn = -0.5d0/dble(ni)**2 + 0.25d0/dist0
         do i = 1,n
            ei = e(i)!*4.35974417d-18/1.60218d-19
!            if ((ei.ge.-0.34d0).and.(ei.le.-0.18d0)) then
            diff = dabs((enn-ei)/enn)
            if (diff.le.0.2d0) then
               ndd = int(dist0)
               write(ndd,*) ei,i,diff,enn
               print*,ei,i,diff,enn
            endif
         enddo
         stop
      endif
      
!     store selection of wavefunctions (eigenvectors)
      do istat = 1,nStates
         nwav = n0(istat)
         do j = 1,n 
            hwav(j,istat) = h(j,nwav)
         enddo
!    work out energy shift from infinite distance
         denn = (e(nwav)-field*dist0-0.25d0/dist0)-en(dble(ni),0.d0,dble(k0(istat)))
!    work out corresponding velocity change from infinite distance
         velacc(istat) = -dsqrt(v0**2-denn*2.d0/mass)
      enddo
!      
!     Tranform Eigenvector to CWDVR basis:
!     INTERPOLATE RADIAL POINTS
      
      do k = 1,nt        !loop over angles
         do j = 1,nx     ! loop over interpolation points (cwdvr)
            jc = nt*(j-1)+k
            phi(j,k) = 0.0d0
            do i = 1,nr  !loop over basis set (laguerre)
               ic = nt*(i-1)+k
               phi(j,k) = phi(j,k)+(hwav(ic,1)*basis(j,i)/dsqrt(b))/a(j)
!   coeff for state ij with angle k and r=j on CWDVR is: 
!   (sum of coeff for state ij over the laguerre basis i=1,nr) x CWDVR weight for r=j 
            enddo
         enddo
      enddo
  
!     interpolate angular functions
      do i = 1,nt
         do j = 1,ntcw
              ptg(i,j) = pbtg(i,j)
         enddo
      enddo
!     coeff of phi (in FBR) x [FBR basis functions(evaluated at ntcw angular points)*gaussian weight]
!     ---> angular DVR with ntcw ang points
      r0 = 0.d0
      r1 = 1.d0
      call dgemm('n','n',nx,ntcw,nt,r1,phi,nx,ptg,nt,r0,p,nx)

!
! normalise wavefunction


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
!
      do i = 1,nx
         do j = 1,ntcw
            psi(i,j) = cmplx(p(i,j), 0.d0, kind=dp)
         enddo
      enddo
      return 
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wavfunc(pbtg,psi,hwav,nx,ntcw,a,nr,basis,nt,istat,nStates)
      use settings, only : b, dp
      implicit none
!--------------------------------------------------------------------------
!     NB x(nx) are the cwdvr grid points and r(nr) are the laguerre points
!     cost(nt) are the ang points for initial diagonalisation and cost(ntcw) and those for propagation
!--------------------------------------------------------------------------
      integer :: nx, ntcw, nr, nt, istat, nStates
      real(dp) :: phi(nx,nt),hwav(nr*nt,nStates),basis(nx,nr)
      real(dp) :: a(nx),ptg(nt,ntcw),pbtg(ntcw,ntcw)
      real(dp) :: p(nx,ntcw)
      complex(dp) psi(nx,ntcw)
      real(dp) :: n
      
      integer :: i, j, k, ic, jc
      real(dp) :: aintegral, sa, r0, r1

      n = nr*nt
      do k = 1,nt        !loop over angles
         do j = 1,nx     ! loop over interpolation points (cwdvr)
            jc = nt*(j-1)+k
            phi(j,k) = 0.0d0
            do i = 1,nr  !loop over basis set (laguerre)
               ic = nt*(i-1)+k
               phi(j,k) = phi(j,k)+(hwav(ic,istat)*basis(j,i)/dsqrt(b))/a(j)
!   coeff for state ij with angle k and r=j on CWDVR is: 
!   (sum of coeff for state ij over the laguerre basis i=1,nr) x CWDVR weight for r=j 
            enddo
         enddo
      enddo
!     interpolate angular functions
      do i = 1,nt
         do j = 1,ntcw
              ptg(i,j) = pbtg(i,j)
         enddo
      enddo
!     coeff of phi (in FBR) x [FBR basis functions(evaluated at ntcw angular points)*gaussian weight]
!     ---> angular DVR with ntcw ang points
      r0 = 0.d0
      r1 = 1.d0
      call dgemm('n','n',nx,ntcw,nt,r1,phi,nx,ptg,nt,r0,p,nx)
!
! normalise wavefunction
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
!
      do i = 1,nx
         do j = 1,ntcw
            psi(i,j) = cmplx(p(i,j), 0.d0, kind=dp)
         enddo
      enddo
      return 
      end

      subroutine lagdvr(n,b,r,t)
      use settings, only : dp
      implicit none
!     -------------------------------------------------------------------------
!     compute the zeroes of Ln(x) to obtain n DVR grid points
!     ans also compute the kinetic energy operator matrix for reg-laguerre DVR
!     -------------------------------------------------------------------------
      integer :: n
      real(dp) :: b
      real(dp) :: r(n),d(n),t(n,n),z(n,n),work(2*n-2)!,wt(n)
      integer :: i, j, ierr
      real(dp) :: obs
      
      do i = 1,n
         r(i) = (2*i-1)
         d(i) = -i
      enddo
!     compute zeros (eigenvalues)
      ierr = 0
      call dstev('V',n,r,d,z,n,work,ierr)
      if (ierr.ne.0) stop 'lagrt | dgtsv'
!
!       do i=1,n
!          rtgwt=z(1,i)
!          wt(i)=rtgwt*sqrt(b)/(dexp(-r(i)/(2.d0)))
!       enddo
!     computer kinetic energy operator matrix
      obs = 1.d0/(b*b)
      do i = 1,n
         do j = 1,n
            if (i.ne.j) then
                if (mod(i-j,2).eq.0) then
                   t(i,j) = ((r(i)+r(j))/((dsqrt(r(j)*r(i)))*((r(i)-r(j))**2)))*obs *0.5d0
                else
                   t(i,j)=((r(i)+r(j))/((dsqrt(r(j)*r(i)))*((r(i)-r(j))**2)))*obs *(-0.5d0)

                endif
            else
                t(i,i)=((4.d0+(4.d0*n+2.d0)*r(i)-r(i)**2)/(12.d0*r(i)**2))*obs *0.5d0
            endif
         enddo
      enddo
      do i = 1,n
         r(i) = r(i)*b
      enddo
      write(6,'(A,i9.2)')'initial diagonalisation nr = ',n
      write(6,'(A,f9.2,/)')'initial diagonalisation grid rmax = ',r(n)
      return
      end

      subroutine lagtrans(x,nx,trans,nr,b)
      use settings, only : dp
      implicit none
!     ---------------------------------------------------------------------------
!     computer value of Baye reg-laguerre basis functions at CWDVR grid points
!     ---------------------------------------------------------------------------
      integer :: nx, nr
      real(dp) :: b
      real(dp) :: x(nx),ts(0:nr),trans(nx)
      integer :: i, j
      real(dp) :: xj
      do i = 1,nx
         ts(0) = 1.d0
         ts(1) = 1.d0-x(i)/b
         do j = 1,nr-1
            xj = dble(j)
            ts(j+1) = (2.d0*xj+1.d0-x(i)/b)*ts(j)/(xj+1.d0)-(xj)*ts(j-1)/(xj+1.d0)
         enddo
         trans(i) = (x(i)/b)*ts(nr)*dexp(-x(i)*0.5d0/b)
      enddo
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine flux2(psi,cost,mrflag,voptic,fg2,fl2,vel)
      use settings, only : nr, nt, dp
      implicit none
      real(dp) :: vel, fl2, fg2
      integer :: mrflag
      complex(dp) :: psi(nr,nt)
      real(dp) :: cost(nt)
      real(dp) :: voptic(nr)
      integer :: k, i
      real(dp) :: dp2, dble
      
      fg2 = 0.d0
      fl2 = 0.d0
      do k = 1,nt
         dp2 = 0.d0
         do i = mrflag,nr
            if (voptic(i).ne.0.d0) then
               dp2 = dp2 + 2.d0*dble(conjg(psi(i,k))*voptic(i)*psi(i,k))/vel
            endif
         enddo
         if (cost(k).gt.0.d0) then
            fg2 = fg2+dp2
         else
            fl2 = fl2+dp2
         endif
      enddo
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine anginterp(npt,ptg,pcost)
      use settings, only : nt, dp
      implicit none
!     --------------------------------------------------------------------------
!     work out angular FBR function at npt angular plotting points
!     --------------------------------------------------------------------------
      integer :: npt
      real(dp) :: pbtg(npt,npt),pcost(npt),pcent(npt),ptwt(npt)
      real(dp) :: ptg(nt,npt)
      integer :: i, j
      call legdvrw(npt,pcost,pbtg,ptwt,pcent)
      do i = 1,nt
         do j = 1,npt
              ptg(i,j) = pbtg(i,j)/sqrt(ptwt(j))
         enddo
      enddo
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

