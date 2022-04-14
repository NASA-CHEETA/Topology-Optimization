# Note that in order to run PARDISO with Multithreading, 
# the environment variable OMP_NUM_THREADS needs to be set to
# the desired number of cpus, before compiling and before running ccx
# e.g. by entering 
#	>export OMP_NUM_THREADS=2
# in terminal to run PARDISO with 2 cpus
# More info in ccx user manual section 2

OPT = -g -O2
CFLAGS = -Wall $(OPT) -fopenmp -DARCH="Linux" -DARPACK -DPARDISO -DMATRIXSTORAGE -DUSE_MT=1 -DMKL_LP64 -I{MKL_INCLUDE}
FFLAGS = -Wall $(OPT) $(OPTIONS) -I$(MKL_INCLUDE) -fopenmp

ARPACK          = /home/prateek/Documents/GitHub/MDO/Dependencies/ARPACK

#ARPACK
#CFLAGS+= -I /usr/local/ARPACK -DARPACK

#SPOOLES
#CFLAGS+= -I /usr/local/SPOOLES.2.2 -DSPOOLES

#BLAS and LAPACK
#CFLAGS+= -llpack -lblas

CC=gcc
FC=gfortran

OPTIONS = -w -fno-second-underscore -fcray-pointer -x f77-cpp-input
FPPSTOP= -x none

#INTEL MKL
MKL_LIB= //opt/intel/oneapi/mkl/2022.0.2/lib/intel64
MKL_INCLUDE=//opt/intel/oneapi/mkl/2022.0.2/include
CMPLR_PATH= //opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin



.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.15.c

OCCXF = $(SCCXF:.f=.o) 
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)


#LIBS= ../TOP_OPT_DEV/ARPACK/libarpack_UBUNTU.a -lpthread -lm -lc
# OS-specific options

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
        LIBS += $(ARPACK)/libarpack_MAC.a
else
        LIBS += $(ARPACK)/libarpack_INTEL.a
        LIBS += -lpthread -lm -lc
endif

caltop_gm: $(OCCXMAIN) caltop_gm.a $(LDFLAGS) $(LIBS)\
        ./date.pl; $(CC) $(CFLAGS) -c ccx_2.15.c; \
        $(FC) $(OPT)p $(FFLAGS) -o $@ $(FPPSTOP) $(OCCXMAIN) caltop_gm.a $(LDFLAGS) $(LIBS)\
         -Wl,--start-group \
        $(MKL_LIB)/libmkl_gf_lp64.a \
        $(MKL_LIB)/libmkl_gnu_thread.a \
        $(MKL_LIB)/libmkl_core.a \
        $(MKL_LIB)/libmkl_blacs_openmpi_lp64.a \
         -Wl,--end-group \
        -L$(CMPLR_PATH) -lgomp -lpthread -lm -ldl

caltop_gm.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
clean:
	rm *.a *.o caltop_gm

