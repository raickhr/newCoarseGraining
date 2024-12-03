program solvePoission
#include "petsc/finclude/petscdmda.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscsys.h"
  use petscsys
  use petscmpi
  use petscksp
  use petscdmda  
  implicit none

  PetscErrorCode :: ierr
  PetscMPIInt :: rank, size
  Vec :: u, b
  Mat :: A
  KSP :: ksp
  PetscInt :: Istart, Iend, Ii, i, j, mx, my
  PetscScalar :: hx, hy, h2x, h2y, rhs
  PetscScalar :: norm

  ! Initialize PETSc
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

  ! Grid dimensions (example: large grid)
  mx = 2600
  my = 600
  hx = 1.0d0 / (mx - 1)
  hy = 1.0d0 / (my - 1)
  h2x = 1.0d0 / (hx * hx)
  h2y = 1.0d0 / (hy * hy)

  ! Create distributed matrix A
  call MatCreate(PETSC_COMM_WORLD, A, ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, mx * my, mx * my, ierr)
  call MatSetFromOptions(A, ierr)
  call MatSetUp(A, ierr)

  ! Create distributed vectors b and u
  call VecCreate(PETSC_COMM_WORLD, b, ierr)
  call VecSetSizes(b, PETSC_DECIDE, mx * my, ierr)
  call VecSetFromOptions(b, ierr)
  call VecDuplicate(b, u, ierr)

  ! Get local ownership range of rows
  call MatGetOwnershipRange(A, Istart, Iend, ierr)

  ! Assemble the matrix A and vector b
  do Ii = Istart, Iend - 1
     i = mod(Ii, mx)
     j = Ii / mx
     rhs = 1.0d0  ! Example RHS value (modify as needed)
     call VecSetValue(b, Ii, rhs * hx * hy, INSERT_VALUES, ierr)

     ! Interior grid points: central difference stencil
     if (i > 0) call MatSetValue(A, Ii, Ii - 1, -h2x, INSERT_VALUES, ierr)
     if (i < mx - 1) call MatSetValue(A, Ii, Ii + 1, -h2x, INSERT_VALUES, ierr)
     if (j > 0) call MatSetValue(A, Ii, Ii - mx, -h2y, INSERT_VALUES, ierr)
     if (j < my - 1) call MatSetValue(A, Ii, Ii + mx, -h2y, INSERT_VALUES, ierr)
     call MatSetValue(A, Ii, Ii, 2.0d0 * (h2x + h2y), INSERT_VALUES, ierr)
  end do

  ! Finalize assembly
  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
  call VecAssemblyBegin(b, ierr)
  call VecAssemblyEnd(b, ierr)

  ! Create the KSP linear solver
  call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
  call KSPSetOperators(ksp, A, A, ierr)
  call KSPSetFromOptions(ksp, ierr)

  ! Solve the linear system A * u = b
  call KSPSolve(ksp, b, u, ierr)

  ! Output the solution norm as a summary
  
  call VecNorm(u, NORM_2, norm, ierr)
  if (rank == 0) then
     write(*,*) 'Solution vector 2-norm:', norm
  end if

  ! Clean up
  call KSPDestroy(ksp, ierr)
  call VecDestroy(b, ierr)
  call VecDestroy(u, ierr)
  call MatDestroy(A, ierr)
  call PetscFinalize(ierr)
end program solvePoission

