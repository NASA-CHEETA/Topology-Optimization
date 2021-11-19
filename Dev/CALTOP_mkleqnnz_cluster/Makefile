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
MKL_LIB= /usr/local/intel/parallel_studio_xe_2018/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
MKL_INCLUDE= /usr/local/intel/parallel_studio_xe_2018/compilers_and_libraries_2018.1.163/linux/mkl/include
CMPLR_PATH= /usr/local/intel/parallel_studio_xe_2018/compilers_and_libraries/linux/lib/intel64_lin


.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.15.c

OCCXF = $(SCCXF:.f=.o) 
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)


LIBS= ../TOP_OPT_DEV/ARPACK/libarpack_UBUNTU.a -lpthread -lm -lc

ccx_2.15_mkleqnnz: $(OCCXMAIN) ccx_2.15_mkleqnnz.a $(LDFLAGS) $(LIBS)\
        ./date.pl; $(CC) $(CFLAGS) -c ccx_2.15.c; \
        $(FC) $(OPT)p $(FFLAGS) -o $@ $(FPPSTOP) $(OCCXMAIN) ccx_2.15_mkleqnnz.a $(LDFLAGS) $(LIBS)\
         -Wl,--start-group \
        $(MKL_LIB)/libmkl_gf_lp64.a \
        $(MKL_LIB)/libmkl_gnu_thread.a \
        $(MKL_LIB)/libmkl_core.a \
        $(MKL_LIB)/libmkl_blacs_openmpi_lp64.a \
         -Wl,--end-group \
        -L$(CMPLR_PATH) -lgomp -lpthread -lm -ldl

ccx_2.15_mkleqnnz.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
clean:
	rm *.a *.o ccx_2.15_mkleqnnz

