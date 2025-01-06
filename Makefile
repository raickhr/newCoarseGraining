CC=mpiicc
FC=mpiifort
CFLAGS=`nc-config --cflags`
FFLAGS=-cpp `nf-config --fflags`
CLIBS=`nc-config --libs`
FLIBS=`nf-config --flibs`
DEBUG= 
#=-g -O0 -traceback -check all -check bounds -check-uninit 

PETSC_DIR=/home/shikhar.rai/myLibraries/petsc
# PETSC_DIR=/opt/anaconda3/envs/fortran/
PETSC_ARCH=

PETSC_INCLUDE=-I${PETSC_DIR}/include
PETSC_LIBS=-L${PETSC_DIR}/lib -lpetsc

NC_INCLUDE=$(nf-config --fflags)
NC_LIBS=`nf-config --flibs`

#-ftrapuv

SRC = $(wildcard *.F90)
OBJ = $(SRC:.F90=.o)

%.o :%.F90
	$(FC) -c $(FFLAGS) $(DEBUG) $< 

main.o: fields.o configurationMod.o gridModule.o mpiwrapper.o filterparallel.o input_data_info.o read_write.o multiGridHelmHoltz.o 
	$(FC) -c $(FFLAGS) $(DEBUG) main.F90 ${PETSC_INCLUDE}

configurationMod.o: kinds.o mpiwrapper.o
	$(FC) -c $(FFLAGS) $(DEBUG) configurationMod.F90

gridModule.o : kinds.o constants.o mpiwrapper.o configurationMod.o
	$(FC) -c $(FFLAGS) $(DEBUG) gridModule.F90

netcdf_io.o : ncdf_wrapper.o gridModule.o
	$(FC) -c $(FFLAGS) $(DEBUG) netcdf_io.F90

interpolation.o : kinds.o
	$(FC) -c $(FFLAGS) $(DEBUG) interpolation.F90

operators.o : kinds.o
	$(FC) -c $(FFLAGS) $(DEBUG) operators.F90

constants.o: kinds.o
	$(FC) -c $(FFLAGS) $(DEBUG) constants.F90

coarsening.o : kinds.o
	$(FC) -c $(FFLAGS) $(DEBUG) coarsening.F90

helmHoltzDecomp.o :	kinds.o mpiwrapper.o operators.o
	$(FC) -c $(FFLAGS) ${PETSC_INCLUDE} $(DEBUG) helmHoltzDecomp.F90 

multiGridHelmHoltz.o :   kinds.o coarsening.o interpolation.o mpiwrapper.o operators.o helmHoltzDecomp.o
	$(FC) -c $(FFLAGS) ${PETSC_INCLUDE} $(DEBUG)  multiGridHelmHoltz.F90 ${PETSC_INCLUDE} ${NC_INCLUDE}

read_write.o: kinds.o ncdf_wrapper.o netcdf_io.o fields.o input_data_info.o gridModule.o	
	$(FC) -c $(FFLAGS) $(DEBUG) read_write.F90

fields.o: kinds.o gridModule.o
	$(FC) -c $(FFLAGS) $(DEBUG) fields.F90	

mpiwrapper.o: mpiwrapper.F90
	$(FC) -c $(FFLAGS) $(DEBUG) mpiwrapper.F90	

filterparallel.o : netcdf_io.o kinds.o mpiwrapper.o fields.o configurationMod.o 
	$(FC) -c $(FFLAGS) $(DEBUG) filterparallel.F90	


all: $(OBJ)
	$(FC) $^ -o main.exe  ${PETSC_LIBS} ${NC_LIBS}
	mv main.exe tests/

clean:
	rm -f *.mod *.o



