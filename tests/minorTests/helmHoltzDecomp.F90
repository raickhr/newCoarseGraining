module helmHoltzDecomp
    #include "petsc/finclude/petscdmda.h"
    #include "petsc/finclude/petscksp.h"
    #include "petsc/finclude/petscsys.h"
    use kinds
    use operators
    use petscsys
    use petscmpi
    use petscksp
    use petscdmda  
    implicit none

    contains

    ! subroutine decomposeHelmholtz(uvel, vvel, psi, phi, dxBottom, dxTop, dyLeft, dyRight, cellArea)

    !     real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, & ! For right hand side
    !                                                      &  dxBottom, dxTop, dyLeft, dyRight, cellArea
    !     real(kind=real_kind), intent(inout), dimension(:,:) :: psi, phi

    !     real(kind=real_kind), allocatable :: divU(:,:), curlU(:,:)
    !     integer :: shapeArr(2), nx, ny


    !     PetscErrorCode :: ierr
    !     PetscMPIInt :: rank, size
    !     Vec :: x_globalOnZero, x_local, y_globalOnZero, y_local
    !     VecScatter :: scatter
    !     Mat :: A
    !     KSP :: ksp
    !     PetscInt :: Istart, Iend, Ii, i, j, mx, my, localSize, globalSize
    !     PetscScalar :: leftCoeff, rightCoeff, topCoeff, bottomCoeff, centerCoeff
    !     PetscScalar :: norm
    !     real(kind=8), pointer :: collected_Pointer(:)
    !     real(kind=real_kind), allocatable :: collected_array(:)

    !     ! Initialize PETSc
    !     call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    !     call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    !     call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

    !     shapeArr = shape(uvel)
    !     nx = shapeArr(1)
    !     ny = shapeArr(2)




    ! end subroutine

    subroutine decomposeHelmholtz(uvel, vvel, psi, phi, dxBottom, dxTop, dyLeft, dyRight, cellArea)

        real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, & ! For right hand side
                                                         &  dxBottom, dxTop, dyLeft, dyRight, cellArea
        real(kind=real_kind), intent(inout), dimension(:,:) :: psi, phi

        real(kind=real_kind), allocatable :: divU(:,:), curlU(:,:)
        

        call calcHozDivVertCurl(uvel, vvel, dxBottom, dxTop, dyLeft, dyRight, cellArea, divU, curlU)

        call solvepoission(psi, -culrU, dxBottom, dxTop, dyLeft, dyRight, cellArea)
        call solvepoission(phi, -divU, dxBottom, dxTop, dyLeft, dyRight, cellArea)

        deallocate(divU, curlU)

    end subroutine

    subroutine solvepoission(phi, RHS, dxBottom, dxTop, dyLeft, dyRight, cellArea)
        real(kind=real_kind), intent(in), dimension(:,:) :: phi, RHS, dxBottom, dxTop, dyLeft, dyRight, cellArea
        integer :: shapeArr(2)

        PetscErrorCode :: ierr
        PetscMPIInt :: rank, size
        Vec :: x_globalOnZero, x_local, y_globalOnZero, y_local
        VecScatter :: scatter
        Mat :: A
        KSP :: ksp
        PetscInt :: Istart, Iend, Ii, i, j, mx, my, localSize, globalSize
        PetscScalar :: leftCoeff, rightCoeff, topCoeff, bottomCoeff, centerCoeff
        PetscScalar :: val, norm

        real(kind=8), pointer :: collected_Pointer(:)
        real(kind=real_kind), allocatable :: collected_array(:)

        ! Initialize PETSc
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

        shapeArr = shape(phi)
        mx = shapeArr(1)
        my = shapeArr(2)

        ! create a distributed vector for phi
        ! Create the source vector (fully defined on rank 0)

        ! The poission equation is solved in the matrix form as [ A ][x] = [y]

        if (rank == 0) then
            call VecCreateSeq(PETSC_COMM_SELF, mx * my, x_globalOnZero, ierr)
            call VecCreateSeq(PETSC_COMM_SELF, mx * my, y_globalOnZero, ierr)
            do Ii = 0, mx * my + 1
                i = mod(Ii, mx) + 1
                j = Ii / mx + 1
                val = phi(i,j)
                call VecSetValue(x_globalOnZero, Ii-1, val, INSERT_VALUES, ierr)

                val = RHS(i,j)
                call VecSetValue(y_globalOnZero, Ii-1, val, INSERT_VALUES, ierr)
            end do
        else
            call VecCreateSeq(PETSC_COMM_SELF, 0, x_globalOnZero, ierr)
            call VecCreateSeq(PETSC_COMM_SELF, 0, y_globalOnZero, ierr)
        end if

        call VecAssemblyBegin(x_globalOnZero, ierr)
        call VecAssemblyEnd(x_globalOnZero, ierr)

        call VecAssemblyBegin(y_globalOnZero, ierr)
        call VecAssemblyEnd(y_globalOnZero, ierr)

        ! Create the distributed destination vector
        call VecCreate(PETSC_COMM_WORLD, x_local, ierr)
        call VecSetSizes(x_local, PETSC_DECIDE, mx * my, ierr)
        ! call VecSetFromOptions(b, ierr) ! commenting this out because using default mpi vector

        ! Create scatter context
        call VecScatterCreateToAll(x_globalOnZero, scatter, x_local, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scatter, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scatter, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)

        ! NOW creating a coefficient distributed matrix A
        call MatCreate(PETSC_COMM_WORLD, A, ierr)
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, mx * my, mx * my, ierr)
        ! call MatSetFromOptions(A, ierr) ! commenting this because we are using the default sparse matrix config
        call MatSetUp(A, ierr)



    end subroutine

    

    

    ! ! Grid dimensions (example: large grid)
    ! mx = 2600
    ! my = 600
    ! hx = 1.0d0 / (mx - 1)
    ! hy = 1.0d0 / (my - 1)
    ! h2x = 1.0d0 / (hx * hx)
    ! h2y = 1.0d0 / (hy * hy)

    ! ! Create distributed matrix A
    ! call MatCreate(PETSC_COMM_WORLD, A, ierr)
    ! call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, mx * my, mx * my, ierr)
    ! call MatSetFromOptions(A, ierr)
    ! call MatSetUp(A, ierr)

    ! ! Create distributed vectors b and u
    ! call VecCreate(PETSC_COMM_WORLD, b, ierr)
    ! call VecSetSizes(b, PETSC_DECIDE, mx * my, ierr)
    ! call VecSetFromOptions(b, ierr)
    ! call VecDuplicate(b, u, ierr)

    ! ! Get local ownership range of rows
    ! call MatGetOwnershipRange(A, Istart, Iend, ierr)

    ! ! Assemble the matrix A and vector b
    ! do Ii = Istart, Iend - 1
    ! i = mod(Ii, mx)
    ! j = Ii / mx
    ! rhs = 1.0d0  ! Example RHS value (modify as needed)
    ! call VecSetValue(b, Ii, rhs * hx * hy, INSERT_VALUES, ierr)

    ! ! Interior grid points: central difference stencil
    ! if (i > 0) call MatSetValue(A, Ii, Ii - 1, -h2x, INSERT_VALUES, ierr)
    ! if (i < mx - 1) call MatSetValue(A, Ii, Ii + 1, -h2x, INSERT_VALUES, ierr)
    ! if (j > 0) call MatSetValue(A, Ii, Ii - mx, -h2y, INSERT_VALUES, ierr)
    ! if (j < my - 1) call MatSetValue(A, Ii, Ii + mx, -h2y, INSERT_VALUES, ierr)
    ! call MatSetValue(A, Ii, Ii, 2.0d0 * (h2x + h2y), INSERT_VALUES, ierr)
    ! end do

    ! ! Finalize assembly
    ! call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    ! call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    ! call VecAssemblyBegin(b, ierr)
    ! call VecAssemblyEnd(b, ierr)

    ! ! Create the KSP linear solver
    ! call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    ! call KSPSetOperators(ksp, A, A, ierr)
    ! call KSPSetFromOptions(ksp, ierr)

    ! ! Solve the linear system A * u = b
    ! call KSPSolve(ksp, b, u, ierr)

    ! ! Output the solution norm as a summary

    ! call VecNorm(u, NORM_2, norm, ierr)
    ! if (rank == 0) then
    ! write(*,*) 'Solution vector 2-norm:', norm
    ! end if

    ! ! Create a local vector for gathering the solution on rank 0
    ! ! if (rank == 0) then
    ! !    call VecCreateSeq(PETSC_COMM_SELF, mx * my, u_local, ierr)
    ! ! else
    ! !    u_local = PETSC_NULL_VEC
    ! ! end if

    ! ! Create a scatter context to gather the global vector to rank 0
    ! call VecScatterCreateToZero(u, scatter, u_local, ierr)
    ! call VecScatterBegin(scatter, u, u_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
    ! call VecScatterEnd(scatter, u, u_local, INSERT_VALUES, SCATTER_FORWARD, ierr)

    ! call VecGetSize(u_local, globalSize, ierr)

    ! if (rank == 0) then
    ! print *, globalSize, mx, my, mx*my
    ! ! Get the array data from u_local
    ! allocate(u_array(globalSize))
    ! call VecGetArrayReadF90(u_local, petsc_array, ierr)
    ! u_array = petsc_array
    ! call write2dVar('poissionSol.nc', 'test', reshape(u_array, (/mx, my/)))
    ! deallocate(u_array)
    ! endif



    ! ! Clean up
    ! call KSPDestroy(ksp, ierr)
    ! call VecDestroy(b, ierr)
    ! call VecDestroy(u, ierr)
    ! call VecDestroy(u_local, ierr)
    ! call MatDestroy(A, ierr)
    ! call VecScatterDestroy(scatter, ierr)
    ! call PetscFinalize(ierr)

end module helmHoltzDecomp

