set -ex
FC=mpiifort

PETSC_DIR=/home/shikhar.rai/myLibraries/petsc
PETSC_ARCH=

PETSC_INCLUDE=-I${PETSC_DIR}/include
PETSC_LIBS="-L${PETSC_DIR}/lib -lpetsc"

NC_INCLUDE=$(nf-config --fflags)
NC_LIBS=$(nf-config --flibs)

FLAGS= #-qopenmp

${FC} ${FLAGS} -c kinds.F90 
${FC} ${FLAGS} -c operator.F90 
${FC} ${FLAGS} -c coarsening.F90 
${FC} ${FLAGS} -c interpolation.F90 
${FC} ${FLAGS} -c forTestReadWrite.F90 ${NC_INCLUDE}
${FC} ${FLAGS} -c forTestMain.F90 ${NC_INCLUDE}
#${FC} ${FLAGS} -c solvePoission.F90 ${PETSC_INCLUDE} ${NC_INCLUDE}

${FC} ${FLAGS} *.o -o test.exe ${PETSC_LIBS} ${NC_LIBS}

# ${FC} ${FLAGS} -c gather_vector.F90 ${PETSC_INCLUDE} #${NC_INCLUDE}
# ${FC} ${FLAGS} *.o -o test.exe ${PETSC_LIBS} #${NC_LIBS}

rm -rf *.o *.mod
