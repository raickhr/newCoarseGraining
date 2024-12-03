set -ex
FC=mpiifort
INCLUDE=-I/home/shikhar.rai/myLibraries/petsc/include
LIBS="-L/home/shikhar.rai/myLibraries/petsc/lib -lpetsc"

PETSC_DIR=/home/shikhar.rai/myLibraries/petsc
PETSC_ARCH=

# INCLUDE+=$(nf-config --fflags)
# LIBS+=$(nf-config --flibs)

FLAGS=-qopenmp

# ${FC} ${FLAGS} -c kinds.F90 ${INCLUDE}
# ${FC} ${FLAGS} -c forTestReadWrite.F90 ${INCLUDE}
${FC} ${FLAGS} -c backUPsolvePoission.F90 ${INCLUDE}

${FC} ${FLAGS} *.o -o test.exe ${LIBS}

rm -rf *.o *.mod
