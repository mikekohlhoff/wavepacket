EXECUTABLE:= ./../bin/cwdvr

intel: FC = ifort
gnu: FC = gfortran

F_SRCS := $(wildcard *.f90)
#F_OBJS := ${F_SRCS:.f90=.o}
BUILD_DIR := ./../build/
F_OBJS = settings.o surface.o helpers.o cwdvr.o
OBJS := $(addprefix $(BUILD_DIR), $(F_OBJS))

intel: LIBRARY_DIR := -L/home/atreju/intel/mkl/lib/intel64
intel: INCLUDE_DIR := -I/home/atreju/intel/mkl/include -I/home/atreju/intel/mkl/include/em64t/lp64

intel: DEBUGFLAGS = -debug all -check all -warn all -debug-parameters all -traceback
intel: FLFLAGS = -ipo -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -Wl,--start-group /home/atreju/intel/mkl/lib/intel64/libmkl_intel_lp64.a /home/atreju/intel/mkl/lib/intel64/libmkl_intel_thread.a /home/atreju/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread
intel : FCFLAGS = -ipo -xHost -fno-alias -debug -stand f08 -132
gnu: FCFLAGS = -llapack -lblas -ffast-math -O3 -std=f2008 -ffree-line-length-none

.PHONY: all clean distclean
all: $(EXECUTABLE)
intel : $(EXECUTABLE)
gnu : $(EXECUTABLE)

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

