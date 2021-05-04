#=======================================================================
# Makefile
#=======================================================================


# CMP = gcc #intel, gcc

ifeq ($(CMP),gcc)
FC = mpif90
FFLAGS =
NETCDFloc = /usr/local
NETCDFlib = -I${NETCDFloc}/include -L${NETCDFloc}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl -lm
else
FC = mpiifort
# FFLAGS = -fpp -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr -I$(MPI_ROOT)/libFFLAGS = -fpp -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr -I$(MPI_ROOT)/lib
NETCDFloc =
NETCDFlib = -lnetcdf -lnetcdff
endif

MODDIR = ./mod
SRCDIR = ./src

### List of all files for the main code
SRC = $(SRCDIR)/interpolate_mod.f90 $(SRCDIR)/particleTracer.f90 $(SRCDIR)/initialize.f90 $(SRCDIR)/load_velocity.f90 $(SRCDIR)/save.f90
OBJ = $(SRC:%.f90=%.o)

###### OPTIONS settins ########
OPT = -I$(SRCDIR) $(NETCDFlib) $(FFLAGS)
LINKOPT = $(NETCDFlib) $(FFLAGS)


# -------------------------------------------------

all: particleTracer

particleTracer : $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJ) $(NETCDFlib)

$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) -c $<
	mv $(@F) ${SRCDIR}
	# mv *.mod ${SRCDIR}

%.o : %.f90
	$(FC) $(FFLAGS) -c $<


.PHONY: clean

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f *.o *.mod particleTracer
