EXECUTABLE:= ./../bin/cwdvr

intel: FC = ifort
gnu: FC = gfortran

F_SRCS := $(wildcard *.f90)
#F_OBJS := ${F_SRCS:.f90=.o}
BUILD_DIR := ./../build/
F_OBJS = settings.o surface.o helpers.o cwdvr.o
OBJS := $(addprefix $(BUILD_DIR), $(F_OBJS))

intel: LIBRARY_DIR := -L/opt/intel/mkl/lib/intel64
intel: INCLUDE_DIR := -I/opt/intel/mkl/include -I/opt/intel/mkl/include/intel64/lp64

intel: DEBUGFLAGS = -debug all -check all -warn all -debug-parameters all -traceback
intel: FLFLAGS = -ipo -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
# AT: 
# intel : FCFLAGS = -ipo -xHost -fno-alias -debug -stand f08 -132
# try speed up (linux)
# fast = -ipo, -O3, -no-prec-div, -static, and -xHost
intel : FCFLAGS = -ipo -xHost -O3 -no-prec-div -ftz -ip -march=native -parallel -openmp -fno-alias -stand f08 -132 -complex-limited-range
intel : FCFLAGS = -ipo -xHost -fno-alias -debug -stand f08 -132 -O2 -march=corei7-avx -complex-limited-range -ip
gnu: FCFLAGS = -llapack -lblas -ffast-math -O3 -std=f2008 -ffree-line-length-none

.PHONY: all clean distclean
all:	$(EXECUTABLE)
intel:	$(EXECUTABLE)
gnu:	$(EXECUTABLE)

$(EXECUTABLE) : $(OBJS)
	$(FC) $^ $(LIBRARY_DIR) $(INCLUDE_DIR) $(FLFLAGS) $(FCFLAGS) -o $(EXECUTABLE)

$(BUILD_DIR)%.o : %.f90 
	$(FC) -c $< $(FCFLAGS) -o $@


clean:
	@- $(RM) $(EXECUTABLE)
	@- $(RM) $(OBJS)
	rm -f *.mod
	rm -f *__genmod*


distclean: clean

