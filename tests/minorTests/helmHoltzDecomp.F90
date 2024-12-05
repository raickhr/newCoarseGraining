module helmHoltzDecomp
#include "petsc/finclude/petscdmda.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscsys.h"
    use kinds
    use mpiwrapper
    use operators
    use petscsys
    use petscmpi
    use petscksp
    use petscdmda  
    implicit none

    contains

    subroutine decomposeHelmholtz(uvel, vvel, psi, phi, dxBottomE, dxTopE, dyLeftE, dyRightE, &
                                & dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
        real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, & ! For right hand side
                                                         &  dxBottomE, dxTopE, dyLeftE, dyRightE, & !at cell edges
                                                         &  dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea  ! stencil joining nodes

        real(kind=real_kind), intent(inout), dimension(:,:) :: psi, phi
        real(kind=real_kind), allocatable :: divU(:,:), curlU(:,:)

        PetscErrorCode :: ierr
        PetscMPIInt :: rank, size

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

        if (rank == 0) call calcHozDivVertCurl(uvel, vvel, dxBottomE, dxTopE, dyLeftE, dyRightE, cellArea, divU, curlU)

        if (rank == 0) print *, 'solving poission equation for toroidal part'
        call solvepoission(psi, -curlU, dxBottomE, dxTopE, dyLeftE, dyRightE, &
                          &dxLeftN, dxRightN, dyBottomN, dyTopN,cellArea)
        if (rank == 0) print *, 'solving poission equation for toroidal part complete'

        if (rank == 0) print *, 'solving poission equation for poloidal part'
        call solvepoission(phi, -divU, dxBottomE, dxTopE, dyLeftE, dyRightE,&
                        &  dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
        if (rank == 0) print *, 'solving poission equation for poloidal part complete'
        if (rank == 0) deallocate(divU, curlU)
    end subroutine

    subroutine solvepoission(phi, RHS, dxBottomE, dxTopE, dyLeftE, dyRightE, &
                             dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
        real(kind=real_kind), intent(in), dimension(:,:) :: RHS, &  ! RHS of the poission equation
                                                            dxBottomE, dxTopE, dyLeftE, dyRightE, & !at cell edges
                                                            dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea  ! stencil joining nodes
        real(kind=real_kind), intent(inout), dimension(:,:) :: phi
        integer :: shapeArr(2)

        PetscErrorCode :: ierr
        PetscMPIInt :: rank, size
        Vec :: x_globalOnZero, x_local, y_globalOnZero, y_local, sol_globalOnZero
        VecScatter :: scatterx, scattery, gather
        IS :: xis_localVec, xis_globalVec, yis_localVec, yis_globalVec
        Mat :: A
        KSP :: ksp
        PetscInt :: Istart, Iend, Ii, i, j, mx, my, localSize, globalSize, low, high
        PetscScalar :: leftCoeff, rightCoeff, topCoeff, bottomCoeff, centerCoeff
        PetscScalar :: val, norm, rel_tol, abs_tol, div_tol
        PetscInt :: max_iter

        real(kind=8), pointer :: collected_xPointer(:)
        real(kind=real_kind), allocatable :: collected_xarray(:)

        ! Initialize PETSc
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

        shapeArr = shape(cellArea)
        mx = shapeArr(1)
        my = shapeArr(2)

        ! create a distributed vector for phi
        ! Create the source vector (fully defined on rank 0)

        ! The poission equation is solved in the matrix form as [ A ][x] = [y]

        call VecCreate(PETSC_COMM_WORLD, x_globalOnZero,  ierr)
        call VecSetFromOptions(x_globalOnZero, ierr)
        call VecCreate(PETSC_COMM_WORLD, y_globalOnZero, ierr)
        call VecSetFromOptions(y_globalOnZero, ierr)

        call VecSetSizes(x_globalOnZero, PETSC_DECIDE, mx * my, ierr)
        call VecSetSizes(y_globalOnZero, PETSC_DECIDE, mx * my, ierr)

        if (rank == 0) then
            print *, taskid
            print *, 'creating vectors of size', mx,'x', my
            
            do Ii = 0, mx * my - 1
                i = mod(Ii, mx)
                j = Ii / mx
                val = phi(i+1,j+1)
                ! print *,'i=',i,'j=', j,'Ii=', Ii,'val =', val
                call VecSetValue(x_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                
                ! if ((i == 0) .and. (j > 0) .and. (j < my - 1)) then  
                !     i = i+1
                ! else if ((i == mx - 1) .and. (j > 0) .and. (j < my - 1)) then
                !     i = i-1
                ! else if ((j == 0) .and. (i > 0) .and. (i < mx - 1)) then
                !     j = j+1
                ! else if ((j == my - 1) .and. (i > 0) .and. (i < mx - 1)) then
                !     j = j-1
                ! end if
                val = RHS(i+1,j+1) * cellArea(i+1,j+1)

                call VecSetValue(y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
            end do

            print *,'i=',i,'j=', j,'Ii=', Ii,'val =', phi(i+1,j+1), RHS(i+1,j+1) * cellArea(i+1,j+1)

        end if

        call VecAssemblyBegin(x_globalOnZero, ierr)
        call VecAssemblyEnd(x_globalOnZero, ierr)

        call VecAssemblyBegin(y_globalOnZero, ierr)
        call VecAssemblyEnd(y_globalOnZero, ierr)


        

        if (rank == 0) print *, 'End vector create'

        

        ! Create the distributed destination vector
        if (rank == 0) print *, 'creating vector for distributing'
        call VecCreate(PETSC_COMM_WORLD, x_local, ierr)
        call VecCreate(PETSC_COMM_WORLD, y_local, ierr)
        
        call VecSetSizes(x_local, PETSC_DECIDE, mx * my, ierr)
        call VecSetSizes(y_local, PETSC_DECIDE, mx * my, ierr)

        call VecSetFromOptions(y_local, ierr)
        call VecSetFromOptions(x_local, ierr)
        if (rank == 0) print *, 'creating vector for distributing COMPLETE'
        
        
        call VecGetOwnershipRange(x_local, low, high, ierr) ! /* low, high are global indices */
        call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, xis_localVec, ierr)

        call VecGetOwnershipRange(y_local, low, high, ierr) ! /* low, high are global indices */
        call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, yis_localVec, ierr)
        
        ! call VecGetOwnershipRange(vec_global, low, high, ierr) ! /* low, high are global indices */
        ! call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, xis_globalVec, ierr)


        ! Create scatter context
        if (rank == 0) print *, 'scattering vectors'
        call VecScatterCreate(x_globalOnZero, xis_localVec, x_local, xis_localVec, scatterx, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scatterx, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scatterx, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'scattering X vector COMPLETE'

        ! Create scatter context
        call VecScatterCreate(y_globalOnZero, yis_localVec, y_local, yis_localVec, scattery, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'scattering vectors COMPLETE'

        ! NOW creating a coefficient distributed matrix A
        call MatCreate(PETSC_COMM_WORLD, A, ierr)
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, mx * my, mx * my, ierr)
        call MatSetFromOptions(A, ierr) ! commenting this because we are using the default sparse matrix config
        call MatSetUp(A, ierr)

        ! Get local ownership range of rows
        call MatGetOwnershipRange(A, Istart, Iend, ierr)

        ! Assemble the matrix A and vector b
        do Ii = Istart, Iend - 1
            i = mod(Ii, mx)
            j = Ii / mx

            ! Interior grid points: central difference stencil
            if ((i > 0) .and. (i < mx - 1) .and. (j > 0) .and. (j < my - 1)) then
                rightCoeff = dyRightE(i,j)/dxRightN(i,j)
                leftCoeff = dyLeftE(i,j)/dxLeftN(i,j)
                topCoeff = dxTopE(i,j)/dyTopN(i,j)
                bottomCoeff = dxBottomE(i,j)/dyBottomN(i,j)
                centerCoeff = -dyLeftE(i,j)/dxLeftN(i,j) - dyRightE(i,j)/dxRightN(i,j) - dxBottomE(i,j)/dyBottomN(i,j) - dxTopE(i,j)/dyBottomN(i,j)

                call MatSetValue(A, Ii, Ii - 1, leftCoeff, ADD_VALUES, ierr)
                call MatSetValue(A, Ii, Ii + 1, rightCoeff, ADD_VALUES, ierr)
                call MatSetValue(A, Ii, Ii - mx, bottomCoeff, ADD_VALUES, ierr) 
                call MatSetValue(A, Ii, Ii + mx, topCoeff, ADD_VALUES, ierr)
                call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)

            else if ((i == 0)) then
                centerCoeff = 1
                rightCoeff = -1
                call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                call MatSetValue(A, Ii, Ii + 1, rightCoeff, ADD_VALUES, ierr)
            else if ((i == mx-1)) then
                centerCoeff = 1
                leftCoeff = -1
                call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                call MatSetValue(A, Ii, Ii - 1, leftCoeff, ADD_VALUES, ierr)
            else if ((j == my-1)) then
                centerCoeff = 1
                bottomCoeff = -1
                call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                call MatSetValue(A, Ii, Ii -mx, rightCoeff, ADD_VALUES, ierr)
            else if ((j == 0)) then
                centerCoeff = 1
                topCoeff = -1
                call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                call MatSetValue(A, Ii, Ii +mx, topCoeff, ADD_VALUES, ierr)
            endif

            ! else if ((i == 0) .and. (j > 0) .and. (j < my - 1)) then
            !     ! left boundary on physical lat, lon, grid
            !     i = i+1
            !     rightCoeff = dyRightE(i,j)/dxRightN(i,j)
            !     leftCoeff = dyLeftE(i,j)/dxLeftN(i,j)
            !     topCoeff = dxTopE(i,j)/dyTopN(i,j)
            !     bottomCoeff = dxBottomE(i,j)/dyBottomN(i,j)
            !     centerCoeff = -dyLeftE(i,j)/dxLeftN(i,j) - dyRightE(i,j)/dxRightN(i,j) - dxBottomE(i,j)/dyBottomN(i,j) - dxTopE(i,j)/dyBottomN(i,j)
                
            !     call MatSetValue(A, Ii, Ii +1 - 1, leftCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii +1 + 1, rightCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii +1 - mx, bottomCoeff, ADD_VALUES, ierr) 
            !     call MatSetValue(A, Ii, Ii +1 + mx, topCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii +1, centerCoeff, ADD_VALUES, ierr)


            ! else if ((i == mx - 1) .and. (j > 0) .and. (j < my - 1)) then
            !     ! right boundary on physical lat, lon, grid
            !     i = i-1
            !     rightCoeff = dyRightE(i,j)/dxRightN(i,j)
            !     leftCoeff = dyLeftE(i,j)/dxLeftN(i,j)
            !     topCoeff = dxTopE(i,j)/dyTopN(i,j)
            !     bottomCoeff = dxBottomE(i,j)/dyBottomN(i,j)
            !     centerCoeff = -dyLeftE(i,j)/dxLeftN(i,j) - dyRightE(i,j)/dxRightN(i,j) - dxBottomE(i,j)/dyBottomN(i,j) - dxTopE(i,j)/dyBottomN(i,j)

            !     call MatSetValue(A, Ii, Ii -1 - 1, leftCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii -1 + 1, rightCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii -1 - mx, bottomCoeff, ADD_VALUES, ierr) 
            !     call MatSetValue(A, Ii, Ii -1 + mx, topCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii -1, centerCoeff, ADD_VALUES, ierr)



            ! else if ((j == 0) .and. (i > 0) .and. (i < mx - 1)) then
            !     ! bottom boundary physical lat, lon, grid
            !     j = j+1
            !     rightCoeff = dyRightE(i,j)/dxRightN(i,j)
            !     leftCoeff = dyLeftE(i,j)/dxLeftN(i,j)
            !     topCoeff = dxTopE(i,j)/dyTopN(i,j)
            !     bottomCoeff = dxBottomE(i,j)/dyBottomN(i,j)
            !     centerCoeff = -dyLeftE(i,j)/dxLeftN(i,j) - dyRightE(i,j)/dxRightN(i,j) - dxBottomE(i,j)/dyBottomN(i,j) - dxTopE(i,j)/dyBottomN(i,j)

            !     call MatSetValue(A, Ii, Ii +mx - 1, leftCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii +mx + 1, rightCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii +mx - mx, bottomCoeff, ADD_VALUES, ierr) 
            !     call MatSetValue(A, Ii, Ii +mx + mx, topCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii +mx, centerCoeff, ADD_VALUES, ierr)


            ! else if ((j == my - 1) .and. (i > 0) .and. (i < mx - 1)) then
            !     ! top boundary last column
            !     j = j-1
            !     rightCoeff = dyRightE(i,j)/dxRightN(i,j)
            !     leftCoeff = dyLeftE(i,j)/dxLeftN(i,j)
            !     topCoeff = dxTopE(i,j)/dyTopN(i,j)
            !     bottomCoeff = dxBottomE(i,j)/dyBottomN(i,j)
            !     centerCoeff = -dyLeftE(i,j)/dxLeftN(i,j) - dyRightE(i,j)/dxRightN(i,j) - dxBottomE(i,j)/dyBottomN(i,j) - dxTopE(i,j)/dyBottomN(i,j)

            !     call MatSetValue(A, Ii, Ii -mx - 1, leftCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii -mx + 1, rightCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii -mx - mx, bottomCoeff, ADD_VALUES, ierr) 
            !     call MatSetValue(A, Ii, Ii -mx + mx, topCoeff, ADD_VALUES, ierr)
            !     call MatSetValue(A, Ii, Ii -mx, centerCoeff, ADD_VALUES, ierr)

            ! else if (i == 0 .and. j == 0 ) then
            !     centerCoeff = 1
            !     call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
            !     centerCoeff = -1
            !     call MatSetValue(A, Ii, Ii+1,centerCoeff, ADD_VALUES, ierr)
            ! else if (i == mx -1 .and. j == 0 ) then
            !     centerCoeff = -1
            !     call MatSetValue(A, Ii, Ii,centerCoeff, ADD_VALUES, ierr)
            !     centerCoeff = 1
            !     call MatSetValue(A, Ii, Ii-1,centerCoeff, ADD_VALUES, ierr)
            
            ! else if (i == 0 .and. j == my - 1 ) then
            !     centerCoeff = 1
            !     call MatSetValue(A, Ii, Ii , centerCoeff, ADD_VALUES, ierr)
            !     centerCoeff = -1
            !     call MatSetValue(A, Ii, Ii+1, centerCoeff, ADD_VALUES, ierr)

            ! else if (i == mx -1 .and. j == my - 1 ) then
            !     centerCoeff = -1
            !     call MatSetValue(A, Ii, Ii , centerCoeff, ADD_VALUES, ierr)
            !     centerCoeff = 1
            !     call MatSetValue(A, Ii, Ii-1, centerCoeff, ADD_VALUES, ierr)
            ! end if
        end do

        ! Finalize assembly
        call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
        
        ! Create the KSP linear solver
        call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
        call KSPSetOperators(ksp, A, A, ierr)
        call KSPSetFromOptions(ksp, ierr)

        ! Set solver tolerances: rel_tol, abs_tol, div_tol, and max iterations
        max_iter = 1000
        rel_tol = 1d-1
        abs_tol = 1d2
        div_tol = 1d5
        call KSPSetTolerances(ksp, rel_tol, abs_tol, div_tol,  max_iter, ierr)

        ! Solve the linear system A * u = b
        call KSPSolve(ksp, y_local, x_local, ierr)

        ! Output the solution norm as a summary

        call VecNorm(x_local, NORM_2, norm, ierr)
        if (rank == 0) then
            write(*,*) 'Solution vector 2-norm:', norm
        end if

        ! Create a scatter context to gather the global vector to rank 0
        call VecScatterCreateToZero(x_local, gather, sol_globalOnZero, ierr)
        if (rank == 0) print *, 'gather created'
        call VecScatterBegin(gather, x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'gather began'
        call VecScatterEnd(gather, x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)

        if (rank == 0) print *, 'gather complete'

        if (rank == 0) then
            allocate(collected_xarray(mx*my))
            call VecGetArrayReadF90(sol_globalOnZero, collected_xPointer, ierr)
            collected_xarray = collected_xPointer
            phi = reshape(collected_xarray, (/mx, my/))
            deallocate(collected_xarray)
        endif

        ! Clean up
        call KSPDestroy(ksp, ierr)
        call VecDestroy(x_local, ierr)
        call VecDestroy(x_globalOnZero, ierr)
        call VecDestroy(y_local, ierr)
        call VecDestroy(y_globalOnZero, ierr)
        call VecDestroy(sol_globalOnZero, ierr)
        call MatDestroy(A, ierr)
        call VecScatterDestroy(scatterx, ierr)
        call VecScatterDestroy(scattery, ierr)
        call VecScatterDestroy(gather, ierr)
        call ISDestroy(xis_localVec, ierr)
        call ISDestroy(xis_globalVec, ierr)
        call ISDestroy(yis_localVec, ierr)
        call ISDestroy(yis_globalVec, ierr)
        call PetscFinalize(ierr)
    end subroutine

end module helmHoltzDecomp

