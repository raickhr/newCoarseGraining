CC=mpiicc
FC=mpiifort
CFLAGS=`nc-config --cflags`
FFLAGS=`nf-config --fflags`
CLIBS=`nc-config --libs`
FLIBS=`nf-config --flibs`
DEBUG=-g -O0 -traceback -check all -check bounds -check-uninit 
#-ftrapuv

SRC = $(wildcard *.F90)
OBJ = $(SRC:.F90=.o)

%.o :%.F90
	$(FC) -c $(FFLAGS) $(DEBUG) $< 

main.o: fields.o configurationMod.o gridModule.o mpiwrapper.o filterparallel.o input_data_info.o read_write.o
	$(FC) -c $(FFLAGS) $(DEBUG) main.F90

configurationMod.o: kinds.o
	$(FC) -c $(FFLAGS) $(DEBUG) configurationMod.F90

gridModule.o : kinds.o constants.o mpiwrapper.o
	$(FC) -c $(FFLAGS) $(DEBUG) gridModule.F90

netcdf_io.o : ncdf_wrapper.o gridModule.o
	$(FC) -c $(FFLAGS) $(DEBUG) netcdf_io.F90

constants.o: kinds.o
	$(FC) -c $(FFLAGS) $(DEBUG) constants.F90

read_write.o: kinds.o ncdf_wrapper.o netcdf_io.o fields.o input_data_info.o gridModule.o	
	$(FC) -c $(FFLAGS) $(DEBUG) read_write.F90

fields.o: kinds.o gridModule.o
	$(FC) -c $(FFLAGS) $(DEBUG) fields.F90	

mpiwrapper.o: mpiwrapper.F90
	$(FC) -c $(FFLAGS) $(DEBUG) mpiwrapper.F90	

all: $(OBJ)
	$(FC) $^ -o main.exe $(FFLAGS) $(FLIBS)
	mv main.exe tests/

clean:
	rm -f *.mod *.o



