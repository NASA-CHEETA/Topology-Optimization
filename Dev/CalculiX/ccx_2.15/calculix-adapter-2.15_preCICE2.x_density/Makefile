# See our wiki for getting the CalculiX dependencies:
# https://github.com/precice/calculix-adapter/wiki/Installation-instructions-for-CalculiX
# Set the following variables before building:
# Path to original CalculiX source (e.g. $(HOME)/ccx_2.15/src )
CCX             = /home/prateek/MDO/nasa-cheeta/AEROELASTIC_DEV/Modules/Dev/CalculiX/ccx_2.15/src_density
# Path to SPOOLES main directory (e.g. $(HOME)/SPOOLES.2.2 )
# SPOOLES         = /home/prateek/MDO/nasa-cheeta/FSI/Modules/SPOOLES.2.2
# Path to ARPACK main directory (e.g. $(HOME)/ARPACK )
ARPACK          = /home/prateek/MDO/nasa-cheeta/AEROELASTIC_DEV/Modules/Dependencies/ARPACK
# Path to yaml-cpp prefix (e.g. $(HOME)/yaml-cpp, should contain "include" and "build")
YAML            = /home/prateek/MDO/nasa-cheeta/AEROELASTIC_DEV/Modules/Dependencies/yaml-cpp-yaml-cpp-0.6.2


# Get the CFLAGS and LIBS from pkg-config (preCICE version >= 1.4.0).
# If pkg-config cannot find the libprecice.pc meta-file, you may need to set the
# path where this is stored into PKG_CONFIG_PATH when building the adapter.
PKGCONF_CFLAGS  = $(shell pkg-config --cflags libprecice)
PKGCONF_LIBS    = $(shell pkg-config --libs libprecice)

# Specify where to store the generated .o files
OBJDIR = bin

# Includes and libs
INCLUDES = \
	-I./ \
	-I./adapter \
	-I$(CCX) \
	-I$(ARPACK) \
	$(PKGCONF_CFLAGS) \
	-I$(YAML)/include
#	-I /home/prateek/software/Precice/include/precice


LIBS = \
	-L$(YAML)/build -lyaml-cpp \
		$(PKGCONF_LIBS) \
		-lstdc++ \

# OS-specific options
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIBS += $(ARPACK)/libarpack_MAC.a
else
	LIBS += $(ARPACK)/libarpack_INTEL.a
	LIBS += -lpthread -lm -lc
endif


#LIBS += /home/prateek/MDO/nasa-cheeta/AEROELASTIC_DEV/Modules/Dependencies/ARPACK/libarpack_INTEL.a

# Compilers and flags
#CFLAGS = -g -Wall -std=c++11 -O0 -fopenmp $(INCLUDES) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE
#FFLAGS = -g -Wall -O0 -fopenmp $(INCLUDES)
OPTIONS = -w -fno-second-underscore -fcray-pointer -x f77-cpp-input
FPPSTOP= -x none


OPT = -g -O2
CFLAGS = -Wall $(OPT) -fopenmp $(INCLUDES) -DARCH="Linux" -DARPACK -DPARDISO -DMATRIXSTORAGE -DUSE_MT=1 -DMKL_LP64 -I{MKL_INCLUDE}
FFLAGS = -Wall $(OPT) $(OPTIONS) -I$(MKL_INCLUDE) -fopenmp $(INCLUDES)

# CFLAGS = -Wall -O3 -fopenmp $(INCLUDES) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DUSE_MT=1

# OS-specific options
ifeq ($(UNAME_S),Darwin)
	CC = /usr/local/bin/gcc
else
	CC = mpicc
endif

#FFLAGS = -Wall -O3 -fopenmp $(INCLUDES)
#FC = mpifort
# FC = mpif90
FC = gfortran



#INTEL MKL
MKL_LIB= //opt/intel/oneapi/mkl/2021.2.0/lib/intel64
MKL_INCLUDE= //opt/intel/oneapi/mkl/2021.2.0/include
CMPLR_PATH= //opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin

# Include a list of all the source files
include $(CCX)/Makefile.inc
SCCXMAIN = ccx_2.15.c

# Append additional sources
SCCXC += nonlingeo_precice.c CCXHelpers.c PreciceInterface.c
SCCXF += getflux.f getkdeltatemp.f



# Source files in this folder and in the adapter directory
$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@
$(OBJDIR)/%.o : adapter/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : adapter/%.cpp
	g++ -std=c++11 -I$(YAML)/include -c $< -o $@ $(LIBS)
	#$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@ $(LIBS)

# Source files in the $(CCX) folder
$(OBJDIR)/%.o : $(CCX)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
$(OBJDIR)/%.o : $(CCX)/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Generate list of object files from the source files, prepend $(OBJDIR)
OCCXF = $(SCCXF:%.f=$(OBJDIR)/%.o)
OCCXC = $(SCCXC:%.c=$(OBJDIR)/%.o)
OCCXMAIN = $(SCCXMAIN:%.c=$(OBJDIR)/%.o)
OCCXC += $(OBJDIR)/ConfigReader.o

$(OBJDIR)/ccx_2.15_mkle: $(OBJDIR) $(OCCXMAIN) $(OBJDIR)/ccx_2.15_mkle.a
	$(FC) $(OPT)p $(FFLAGS) -o $@ $(FPPSTOP) $(OCCXMAIN) $(OBJDIR)/ccx_2.15_mkle.a $(LDFLAGS) $(LIBS)\
         -Wl,--start-group \
        $(MKL_LIB)/libmkl_gf_lp64.a \
        $(MKL_LIB)/libmkl_gnu_thread.a \
        $(MKL_LIB)/libmkl_core.a \
        $(MKL_LIB)/libmkl_blacs_openmpi_lp64.a \
         -Wl,--end-group \
        -L$(CMPLR_PATH) -lgomp -lpthread -lm -ldl

#$(OBJDIR)/ccx_2.15_mkle: $(OBJDIR) $(OCCXMAIN) $(OBJDIR)/ccx_2.15_mkle.a $(LDFLAGS) $(LIBS)\
				./date.pl; $(CC) $(CFLAGS) -c $(CC)/ccx_2.15.c; \
				$(FC) $(OPT)p $(FFLAGS) -o $@ $(FPPSTOP) $(OCCXMAIN) $(OBJDIR)/ccx_2.15_mkle.a $(LDFLAGS) $(LIBS)\
					-Wl,--start-group \
				 $(MKL_LIB)/libmkl_gf_lp64.a \
				 $(MKL_LIB)/libmkl_gnu_thread.a \
				     $(MKL_LIB)/libmkl_core.a \
				        $(MKL_LIB)/libmkl_blacs_openmpi_lp64.a \
				         -Wl,--end-group \
				        -L$(CMPLR_PATH) -lgomp -lpthread -lm -ldl

#$(OBJDIR)/ccx_preCICE: $(OBJDIR) $(OCCXMAIN) $(OBJDIR)/ccx_2.15.a $(SPOOLES)/MT/src/spoolesMT.a
#	$(FC) -fopenmp -Wall -O3 -o $@ $(OCCXMAIN) $(OBJDIR)/ccx_2.15.a $(LIBS)


$(OBJDIR)/ccx_2.15_mkle.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/ccx_2.15_mkle.a $(OBJDIR)/ccx_2.15_mkle
