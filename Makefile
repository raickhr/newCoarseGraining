CC=mpiicc
FC=mpiifort
CFLAGS=`nc-config --cflags`
FFLAGS=`nf-config --fflags`
CLIBS=`nc-config --libs`
FLIBS=`nf-config --flibs`

SRC = $(wildcard *.F90)
OBJ = $(SRC:.F90=.o)

%.o :%.F90
	$(FC) -c $(FFLAGS) $< 

main.o: fields.o configurationMod.o gridModule.o mpiwrapper.o filterparallel.o
	$(FC) -c $(FFLAGS) main.F90

gridModule.o : kinds.o constants.o netcdf_io.o
	$(FC) -c $(FFLAGS) gridModule.F90

netcdf_io.o : ncdf_wrapper.o
	$(FC) -c $(FFLAGS) netcdf_io.F90

constants.o: kinds.o
	$(FC) -c $(FFLAGS) constants.F90	

fields.o: kinds.o gridModule.o
	$(FC) -c $(FFLAGS) fields.F90	

all: $(OBJ)
	$(FC) $^ -o main.exe $(FFLAGS) $(FLIBS)

clean:
	rm -f *.mod *.o