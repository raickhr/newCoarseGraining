program gather_vector
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscsys.h"
   use petscsys
   use petscvec
   implicit none
 
   PetscInt    :: i, N, M, low, high, onlylocal
   PetscMPIInt :: size, rank
   Vec         :: x, y
   VecScatter  :: vscat
   PetscErrorCode :: ierr
   PetscScalar :: val

   N = 10
   M = PETSC_DECIDE
   onlylocal = 2
 
   ! Initialize PETSc
   call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
   call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
   !call PetscOptionsGetInt(PETSC_NULL_CHARACTER, PETSC_NULL_CHARACTER, "-n", M, ierr)
 
   ! Trigger special case in VecScatterCreateToAll to deal with the one-to-all pattern
   !call PetscOptionsGetInt(PETSC_NULL_CHARACTER, PETSC_NULL_CHARACTER, "-onlylocal", onlylocal, ierr)
   ! if (onlylocal >= 0 .and. onlylocal < size) then
   !    if (rank == onlylocal) then
   !       M = N
   !    else
   !       M = 0
   !    end if
   ! end if
 
   ! Create vector x
   call VecCreate(PETSC_COMM_WORLD, x, ierr)
   call VecSetFromOptions(x, ierr)
   call VecSetSizes(x, PETSC_DECIDE, N, ierr)
   call VecGetOwnershipRange(x, low, high, ierr)
   call PetscObjectSetName(x, "x", ierr)

   ! if (rank == 0) then
   !    call VecCreate(PETSC_COMM_WORLD, y, ierr)
   !    call VecSetFromOptions(y, ierr)
   !    call VecSetSizes(y, N, N, ierr)
   !    call PetscObjectSetName(y, "y", ierr)
   ! endif
 
   ! -------------------------------------
   ! VecScatterCreateToZero
   ! -------------------------------------
   ! MPI vec x = [0, 1, 2, ..., N-1]
   !val = 100
   do i = low, high-1
      val = i
      call VecSetValue(x, i, val, INSERT_VALUES, ierr)
   end do
   
   !call VecSet(x, val, ierr)

   call VecAssemblyBegin(x, ierr)
   call VecAssemblyEnd(x, ierr)

   call PetscPrintf(PETSC_COMM_WORLD, "value of x in each processor after assigning\n", ierr)
   call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)
 
   call PetscPrintf(PETSC_COMM_WORLD, "Testing VecScatterCreateToZero\n", ierr)
   call VecScatterCreateToZero(x, vscat, y, ierr)
   call PetscObjectSetName(y, "y", ierr)
 
   ! Test PetscSFBcastAndOp with op = MPI_REPLACE, which does y = x on rank 0
   call VecScatterBegin(vscat, x, y, INSERT_VALUES, SCATTER_FORWARD, ierr)
   call VecScatterEnd(vscat, x, y, INSERT_VALUES, SCATTER_FORWARD, ierr)
   if (rank == 0) call VecView(y, PETSC_VIEWER_STDOUT_SELF, ierr)
 
   ! Test PetscSFBcastAndOp with op = MPI_SUM, which does y += x
   call VecScatterBegin(vscat, x, y, ADD_VALUES, SCATTER_FORWARD, ierr)
   call VecScatterEnd(vscat, x, y, ADD_VALUES, SCATTER_FORWARD, ierr)
   if (rank == 0) call VecView(y, PETSC_VIEWER_STDOUT_SELF, ierr)
 
   ! Test PetscSFReduce with op = MPI_REPLACE, which does x = y
   call VecScatterBegin(vscat, y, x, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call VecScatterEnd(vscat, y, x, INSERT_VALUES, SCATTER_REVERSE, ierr)
   call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)
 
   ! Test PetscSFReduce with op = MPI_SUM, which does x += y on x's local part on rank 0
   call VecScatterBegin(vscat, y, x, ADD_VALUES, SCATTER_REVERSE_LOCAL, ierr)
   call VecScatterEnd(vscat, y, x, ADD_VALUES, SCATTER_REVERSE_LOCAL, ierr)
   call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)
 
   call VecDestroy(x, ierr)
   call VecDestroy(y, ierr)
   call VecScatterDestroy(vscat, ierr)
 
   ! -------------------------------------
   ! VecScatterCreateToAll
   ! -------------------------------------
   ! do i = low, high - 1
   !    call VecSetValue(x, i, i, INSERT_VALUES, ierr)
   ! end do
   ! call VecAssemblyBegin(x, ierr)
   ! call VecAssemblyEnd(x, ierr)
 
   ! call PetscPrintf(PETSC_COMM_WORLD, "Testing VecScatterCreateToAll\n", ierr)
 
   ! call VecScatterCreateToAll(x, vscat, y, ierr)
   ! call PetscObjectSetName(y, "y", ierr)
 
   ! ! Test PetscSFBcastAndOp with op = MPI_REPLACE, which does y = x on all ranks
   ! call VecScatterBegin(vscat, x, y, INSERT_VALUES, SCATTER_FORWARD, ierr)
   ! call VecScatterEnd(vscat, x, y, INSERT_VALUES, SCATTER_FORWARD, ierr)
   ! if (rank == 0) call VecView(y, PETSC_VIEWER_STDOUT_SELF, ierr)
 
   ! ! Test PetscSFBcastAndOp with op = MPI_SUM, which does y += x
   ! call VecScatterBegin(vscat, x, y, ADD_VALUES, SCATTER_FORWARD, ierr)
   ! call VecScatterEnd(vscat, x, y, ADD_VALUES, SCATTER_FORWARD, ierr)
   ! if (rank == 0) call VecView(y, PETSC_VIEWER_STDOUT_SELF, ierr)
 
   ! ! Test PetscSFReduce with op = MPI_REPLACE, which does x = y
   ! call VecScatterBegin(vscat, y, x, INSERT_VALUES, SCATTER_REVERSE, ierr)
   ! call VecScatterEnd(vscat, y, x, INSERT_VALUES, SCATTER_REVERSE, ierr)
   ! call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)
 
   ! ! Test PetscSFReduce with op = MPI_SUM, which does x += size*y
   ! call VecScatterBegin(vscat, y, x, ADD_VALUES, SCATTER_REVERSE, ierr)
   ! call VecScatterEnd(vscat, y, x, ADD_VALUES, SCATTER_REVERSE, ierr)
   ! call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)
 
   call VecDestroy(x, ierr)
   call VecDestroy(y, ierr)
   call VecScatterDestroy(vscat, ierr)
 
   call PetscFinalize(ierr)
 end program gather_vector
 