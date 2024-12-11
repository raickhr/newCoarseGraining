set -ex
FC=mpif90

PETSD_DIR="/home/shikhar.rai/myLibraries/petsc"
PETSC_DIR="/opt/anaconda3/envs/fortran/"
PETSC_ARCH=

PETSC_INCLUDE=-I${PETSC_DIR}/include
PETSC_LIBS="-L${PETSC_DIR}/lib -lpetsc"

NC_INCLUDE=$(nf-config --fflags)
NC_LIBS=$(nf-config --flibs)

FLAGS= #-qopenmp

${FC} ${FLAGS} -c kinds.F90 
${FC} ${FLAGS} -c operator.F90 
${FC} ${FLAGS} -c interpolation.F90 
${FC} ${FLAGS} -c coarsening.F90 
${FC} ${FLAGS} -c mpiwrapper.F90
${FC} ${FLAGS} -c forTestReadWrite.F90 ${NC_INCLUDE}
${FC} ${FLAGS} -c helmHoltzDecomp.F90 ${PETSC_INCLUDE} ${NC_INCLUDE}
${FC} ${FLAGS} -c multiGridHelmHoltz.F90 ${PETSC_INCLUDE} ${NC_INCLUDE}
${FC} ${FLAGS} -c forTestMain.F90 ${PETSC_INCLUDE} ${NC_INCLUDE}

${FC} ${FLAGS} *.o -o test.exe ${PETSC_LIBS} ${NC_LIBS}

################ scatter test ##################3
# ${FC} ${FLAGS} -c scatter_test.F90 ${PETSC_INCLUDE} #${NC_INCLUDE}
# ${FC} ${FLAGS} *.o -o test.exe ${PETSC_LIBS} #${NC_LIBS}

################ gather test ##################3
# ${FC} ${FLAGS} -c gather_vector.F90 ${PETSC_INCLUDE} #${NC_INCLUDE}
# ${FC} ${FLAGS} *.o -o test.exe ${PETSC_LIBS} #${NC_LIBS}

rm -rf *.o *.mod