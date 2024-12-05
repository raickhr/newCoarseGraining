program scatter_vector
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscis.h"
    use petscvec
    implicit none
  
    Vec :: vec_global, vec_local
    IS :: is_localVec, is_globalVec
    VecScatter :: scatter
    PetscErrorCode :: ierr
    PetscMPIInt :: rank, size, low, high
    PetscInt :: n_global, n_local, start, end, i
    PetscScalar :: val

    real, allocatable :: arr(:)
  
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  
    ! Get the rank and size of the communicator
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
  
    ! Ensure the number of processors is exactly 5
    ! if (size /= 5) then
    !   if (rank == 0) write(*, *) 'Error: This program must be run with exactly 5 processors!'
    !   call PetscFinalize(ierr)
    !   stop
    ! end if
  
    ! Define the global size of the vector
    n_global = 375 * 75  ! Example: Total size of the vector
    !n_local = n_global / size  ! Local size on each processor
  
    ! Create a global vector
    
    call VecCreate(PETSC_COMM_WORLD, vec_global, ierr)
    ! if (rank == 0) then 
    !     call VecSetSizes(vec_global, n_global, n_global, ierr)
    ! else
    !     call VecSetSizes(vec_global, 0, n_global, ierr)
    ! end if
    
    ! Initialize the global vector on rank 0
    if (rank == 0) then
        allocate(arr(18))
        call random_number(arr)
    endif


    call VecSetSizes(vec_global, PETSC_DECIDE, n_global, ierr)
    call VecSetFromOptions(vec_global, ierr)

    if (rank == 0) then
        
        do i = 1, n_global
            val = i  ! Assign values 1, 2, ..., N
            call VecSetValue(vec_global, i-1, val, INSERT_VALUES, ierr)  ! PETSc uses 0-based indexing
        end do
    ! else
    !     call VecSetSizes(vec_global, 0, n_global, ierr)
    !     call VecSetFromOptions(vec_global, ierr)
    end if
  
    ! Finalize the global vector assembly
    call VecAssemblyBegin(vec_global, ierr)
    call VecAssemblyEnd(vec_global, ierr)

    if (rank == 0) print *, 'Assembly done ! vec created '
  
    ! Create a local vector for each processor
    call VecCreate(PETSC_COMM_WORLD, vec_local, ierr)
    call VecSetSizes(vec_local, PETSC_DECIDE, n_global, ierr)
    call VecSetFromOptions(vec_local, ierr)

    call VecGetOwnershipRange(vec_local, low, high, ierr) ! /* low, high are global indices */
    call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, is_localVec, ierr)

    call VecGetOwnershipRange(vec_global, low, high, ierr) ! /* low, high are global indices */
    call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, is_globalVec, ierr)

    ! if ( rank == 0) then
    !     call ISCreateStride(PETSC_COMM_SELF, n_global, 0, 1, is_globalVec, ierr)
    ! else
    !     is_globalVec = PETSC_NULL_IS
    ! endif
    !call ISDuplicate(is_localVec, is_globalVec, ierr)
    
    print *, rank, n_global, low, high, is_localVec

    ! Create a scatter context
    call VecScatterCreate(vec_global, is_globalVec, vec_local, is_localVec, scatter, ierr)
  
    ! Scatter the global vector to local vectors on each processor
    call VecScatterBegin(scatter, vec_global, vec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecScatterEnd(scatter, vec_global, vec_local, INSERT_VALUES, SCATTER_FORWARD, ierr)

    if (rank == 0) print *, 'Scatter Complete!'
  
    ! Each processor prints its portion of the vector
    !if (rank == 0 ) then
        !write(*, '(A,I2)') 'Processor rank: ', rank
        !call VecView(vec_global, PETSC_VIEWER_STDOUT_SELF, ierr)
        call VecView(vec_local, PETSC_VIEWER_STDOUT_SELF, ierr)
    !endif
    if (rank == 1) print *, arr
  
    ! Clean up
    !call VecScatterDestroy(scatter, ierr)
    call VecDestroy(vec_global, ierr)
    !call VecDestroy(vec_local, ierr)
    call PetscFinalize(ierr)
  
  end program scatter_vector
  