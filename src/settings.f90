module settings
  implicit none
  integer, save :: mm = 3 ! the number of states to run calcs for
  integer, save :: ni     ! quantum number of interest
  integer, allocatable, save :: n0(:)
  integer, allocatable, save :: k0(:)
  logical, save :: ntest, mplot, meanfield ! options
  double precision, save :: mass
  double precision, save :: field, v0, dist0, min_pop, timestep
  integer, save :: nr, nt, nx
  double precision, save :: dk, z, zacc, xmin

  contains
    subroutine loadSettings()
      allocate(n0(mm))
      allocate(k0(mm))
      
      
      ntest = .FALSE.         ! ntest=1 stops calculation after initial diagonalisation, use this option if you want to find the state number first!
                              ! ntest=0 diagonalises and the propagates the wf      
      meanfield = .TRUE.     ! meanfield = 1 for meanfield calc, meanfield=0 for constant velocity calc
      mplot = .FALSE.         ! output wavefunction for plotting? (note: output files can get big! don't just plot out wf willy nilly!)


      ni = 3            ! principal quantum number of interest
      n0 = (/4, 5, 6/)  ! array index of interest (run with ntest = 1 to determine this number)
      k0 = (/-2, 0, 2/) ! k-state of interest
      
      mass = 1836.d0    ! mass of ion core
      
      min_pop = 1e-2    ! stop when population less than 'min_pop'
      timestep = 1.0d0   ! timestep
      
      field = 0.d0      ! electric field
      v0 = -3d-4        ! velocity at inifinite distance or the velocity for constant vel calc
      
      dist0 = 6.d0*ni**2  ! initial distance to start propagation (where the initial diagonalisation is carried out)
      
!     CWDVR parameters (the propagation grid)
      nt = 20      !no. of cwdvr angular points, can be greater than ntb if desired
      dk = 3.0d0   !coulomb wave parameter: inc. dk-> 1.more points, 2.smaller sep further out, 3.more even distribution at larger dist
      z = 50.d0    !inc.  z-> smaller the first grid point  

!     parameters for initial CWDVR grid point search
      zacc = 1.d-8 ! newton-raphson accuracy
      xmin = 8.d-3 ! Lower bound of 1st zero
      nx = 100000  ! Number of grid points to scan for zero, may need to increase this if using high dk parameter    

      
    end subroutine loadSettings
  
end module settings
