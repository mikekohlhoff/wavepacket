module settings
  implicit none
  integer, save :: mm = 3 ! the number of states to run calcs for
  integer, save :: ni     ! quantum number of interest
  integer, allocatable, save :: n0(:)
  integer, allocatable, save :: k0(:)
  logical, save :: ntest, mplot, meanfield ! options
  double precision, save :: mass
  double precision, save :: field, v0

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
      
      field = 0.d0      ! electric field
      v0 = -3d-4        ! velocity at inifinite distance or the velocity for constant vel calc
      
      
    end subroutine loadSettings
  
end module settings
