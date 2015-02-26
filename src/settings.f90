module settings
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer, save :: nStates     ! the number of states to run calcs for
  integer, save :: ni     ! quantum number of interest
  integer, allocatable, save :: n0(:)
  integer, allocatable, save :: k0(:)
  logical, save :: TEST, MPLOT, MEANFIELD ! options
  real(dp), save :: mass
  real(dp), save :: field, v0, dist0, min_pop, timestep
  integer, save :: nt, nx, npt, nr
  real(dp), save :: dk, z, zacc, xmin
  real(dp), parameter :: pi = dacos(-1.d0)
  real(dp), save :: b
  integer, save :: nrb, ntb
  real(dp), save :: rmid, rabs
  integer, save :: nwstep, nfstep, npstep
  real(dp), save :: pminr, eminr, xmax, eta1, dx, xlmin, xlmax

  contains
    subroutine loadSettings()
      NAMELIST/CONTROL/ nStates, TEST, MEANFIELD, MPLOT
      NAMELIST/STATE/ ni, n0, k0
      NAMELIST/ENVIRONMENT/ field, v0
      NAMELIST/GRID/ nt, dk, z, b, nrb, ntb
            
      OPEN(UNIT=1, FILE='config.info')
      READ(1, NML=CONTROL)
      
      allocate(n0(nStates))
      allocate(k0(nStates))

      READ(1, NML=STATE)
      READ(1, NML=ENVIRONMENT)
      READ(1, NML=GRID)

      CLOSE(1)
        
      print*,'----------------------------------------------------------------------------------'
      write(6, *) 'velocity: ', v0
      write(6, *) 'field: ', field
      print*,'---------------------------------------'
      
      ! parameters not (yet) in the config file
      dist0 = 6.d0*ni**2  ! initial distance to start propagation (where the initial diagonalisation is carried out)
      
      mass = 1836.d0    ! mass of ion core
      
      min_pop = 1e-2    ! stop when population less than 'min_pop'
      timestep = 1.0d0  ! timestep
      

!     parameters for initial CWDVR grid point search
      zacc = 1.d-8 ! newton-raphson accuracy
      xmin = 8.d-3 ! Lower bound of 1st zero
      nx = 100000  ! Number of grid points to scan for zero, may need to increase this if using high dk parameter    
      
      npt = nt     ! no. of angular points in outputting the wavefunction, setting more than nt will give interpolated values 
      
!     lots more definitions of unclear purpose making this all awefully messy

!     absorbing potential parameters
      rmid = dist0 !absorbing boundary position on radial grid, AND where flux detector plane sits, recommended to set it to dist0
      rabs = 30.d0 !width of absorbing boundary, needs to be wide enough to absorb the lowest energy components

!     outputting parameters
      nwstep = 100 ! number of time steps (+1)  between each output
      nfstep = 1   ! number of time steps (+1) between evaluation of meanfield
      npstep = 500 ! if outputting wf, the time steps between succesive outputs

!     absorbing potential parameters
      pminr = 2.d0*pi/rabs
      eminr = 0.5d0*pminr**2

!     CWDVR parameters and parameters for scanning CWDVR points
      xmax = (rmid+rabs)*dk ! rmid + rabs = maximum radial point of grid
      ! /common/wave/
      eta1 = -z/dk
      dx = (xmax-xmin)/nx 
      xlmin = 0
      xlmax = 0

    end subroutine loadSettings
  
end module settings
