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
                    
                    val = RHS(i+1,j+1) * cellArea(i+1,j+1)
                    if (i == 0 .or. i == mx -1 .or. j == 0 .or. j == my -1) val = 0.0
    
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
                    rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
                    leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
                    topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
                    bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
                    centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)    - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
                            &     -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)
    
                    call MatSetValue(A, Ii, Ii - 1, leftCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, Ii + 1, rightCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, Ii - mx, bottomCoeff, ADD_VALUES, ierr) 
                    call MatSetValue(A, Ii, Ii + mx, topCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
    
                else if ((i == 0)) then
                    centerCoeff = -1
                    rightCoeff = 1
                    call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, Ii + 1, rightCoeff, ADD_VALUES, ierr)
                else if ((i == mx-1)) then
                    centerCoeff = -1
                    leftCoeff = 1
                    call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, Ii - 1, leftCoeff, ADD_VALUES, ierr)
                else if ((j == my-1)) then
                    centerCoeff = -1
                    bottomCoeff = 1
                    call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, Ii -mx, rightCoeff, ADD_VALUES, ierr)
                else if ((j == 0)) then
                    centerCoeff = -1
                    topCoeff = 1
                    call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, Ii +mx, topCoeff, ADD_VALUES, ierr)
                endif
            end do
    
            ! Finalize assembly
            call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
            
            ! Create the KSP linear solver
            call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
            call KSPSetOperators(ksp, A, A, ierr)
            call KSPSetFromOptions(ksp, ierr)
    
            ! Set solver tolerances: rel_tol, abs_tol, div_tol, and max iterations
            max_iter = 5000
            rel_tol = 1d-4
            abs_tol = 1d-4
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
    
    
        subroutine decomposeHelmholtz_2(uvel, vvel, psi, phi, dxBottomE, dxTopE, dyLeftE, dyRightE, &
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
    
            if (rank == 0) print *, 'solving poission equation'
            call solvepoissionBig(psi, phi, uvel, vvel, dxBottomE, dxTopE, dyLeftE, dyRightE, &
                                  &dxLeftN, dxRightN, dyBottomN, dyTopN,cellArea)
            if (rank == 0) print *, 'solving poission equation'
    
            if (rank == 0) deallocate(divU, curlU)
        end subroutine
    
        
        subroutine solvepoissionBig(phi, psi, uvel, vvel, dxBottomE, dxTopE, dyLeftE, dyRightE, &
            dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
            real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, &  ! RHS of the poission equation
                                        dxBottomE, dxTopE, dyLeftE, dyRightE, & !at cell edges
                                        dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea  ! stencil joining nodes
            real(kind=real_kind), intent(inout), dimension(:,:) :: phi, psi
            integer :: shapeArr(2)
    
            PetscErrorCode :: ierr
            PetscMPIInt :: rank, size
            Vec :: x_globalOnZero, x_local, y_globalOnZero, y_local, sol_globalOnZero
            VecScatter :: scatterx, scattery, gather
            IS :: xis_localVec, xis_globalVec, yis_localVec, yis_globalVec
            Mat :: A
            KSP :: ksp
            PetscInt :: Istart, Iend, Jstart, Jend, Ii, Jj, i, j, colIndex, &
                        colIndex2, mx, my, localSize, globalSize, low, high
            PetscScalar :: leftCoeff, rightCoeff, topCoeff, bottomCoeff, centerCoeff
            PetscScalar :: val, norm, rel_tol, abs_tol, div_tol
            PetscInt :: max_iter
    
            real(kind=8), pointer :: collected_xPointer(:)
            real(kind=real_kind), allocatable :: collected_xarray(:), divU(:,:), curlU(:,:)
            logical :: inColRange
    
            ! Initialize PETSc
            call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
            call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
            call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
    
            if (rank == 0) call calcHozDivVertCurl(uvel, vvel, dxBottomE, dxTopE, dyLeftE, dyRightE, cellArea, divU, curlU)
    
            shapeArr = shape(cellArea)
            mx = shapeArr(1)
            my = shapeArr(2)
    
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SYSTEM OF EQUATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !          
            !         |    -dx          dy      | |       | |   u    |
            !         |    -dy         -dx      | |  phi  | |   v    |
            !         | laplacian        0      | |  psi  | | -div   | 
            !         |     0        laplacian  | |       | | -curl  |
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
            ! create a distributed vector for phi
            ! Create the source vector (fully defined on rank 0)
    
            ! The poission equation is solved in the matrix form as [ A ][x] = [y]
    
            call VecCreate(PETSC_COMM_WORLD, x_globalOnZero,  ierr)
            call VecSetFromOptions(x_globalOnZero, ierr)
            call VecCreate(PETSC_COMM_WORLD, y_globalOnZero, ierr)
            call VecSetFromOptions(y_globalOnZero, ierr)
    
            call VecSetSizes(x_globalOnZero, PETSC_DECIDE, 2*mx * my, ierr)
            call VecSetSizes(y_globalOnZero, PETSC_DECIDE, 4*mx * my, ierr)
    
            if (rank == 0) then
                print *, taskid
                print *, 'creating vectors of size', mx,'x', my
    
                do Ii = 0, mx * my - 1
                    i = mod(Ii, mx)
                    j = Ii / mx
                    
    
                    !!!!!!!!!!!!!!!!!!!!!!!!!! setting phi and psi  !!!!!!!!!!!!!!!!!!!!!!!
                    val = phi(i+1,j+1)
                    call VecSetValue(x_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                    val = psi(i+1,j+1)
                    call VecSetValue(x_globalOnZero, Ii + mx * my, val, INSERT_VALUES, ierr)
    
                    !!!!!!!!!!!!!!!!!!!!!!!!!! setting RHS !!!!!!!!!!!!!!!!!!!!!!!
                    val = uvel(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                    val = vvel(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii + (mx * my), val, INSERT_VALUES, ierr)
                    
                    val = -divU(i+1, j+1)  * cellArea(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii + 2*(mx * my), val, INSERT_VALUES, ierr)
                    val = -curlU(i+1, j+1) * cellArea(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii + 3*(mx * my), val, INSERT_VALUES, ierr)
    
                end do
    
                print *, 'LHS and RHS vectors assigned values Set'
    
            end if
    
            call VecAssemblyBegin(x_globalOnZero, ierr)
            call VecAssemblyEnd(x_globalOnZero, ierr)
    
            call VecAssemblyBegin(y_globalOnZero, ierr)
            call VecAssemblyEnd(y_globalOnZero, ierr)
    
            ! Create the distributed destination vector
            if (rank == 0) print *, 'creating vector for distributing across processors'
    
            call VecCreate(PETSC_COMM_WORLD, x_local, ierr)
            call VecCreate(PETSC_COMM_WORLD, y_local, ierr)
    
            call VecSetSizes(x_local, PETSC_DECIDE, 2*mx * my, ierr)
            call VecSetSizes(y_local, PETSC_DECIDE, 4*mx * my, ierr)
    
            call VecSetFromOptions(y_local, ierr)
            call VecSetFromOptions(x_local, ierr)
    
            if (rank == 0) print *, 'creating vector for distributing across processors COMPLETE'
    
    
            call VecGetOwnershipRange(x_local, low, high, ierr) ! /* low, high are global indices */
            call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, xis_localVec, ierr)
    
            call VecGetOwnershipRange(y_local, low, high, ierr) ! /* low, high are global indices */
            call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, yis_localVec, ierr)
    
            ! Create scatter context
            if (rank == 0) print *, 'scattering LHS vector  ... '
            call VecScatterCreate(x_globalOnZero, xis_localVec, x_local, xis_localVec, scatterx, ierr)
    
            ! Perform the scatter operation
            call VecScatterBegin(scatterx, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call VecScatterEnd(scatterx, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
            if (rank == 0) print *, 'scattering LHS vector COMPLETE'
    
            if (rank == 0) print *, 'scattering RHS vector  ... '
            ! Create scatter context
            call VecScatterCreate(y_globalOnZero, yis_localVec, y_local, yis_localVec, scattery, ierr)
    
            ! Perform the scatter operation
            call VecScatterBegin(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call VecScatterEnd(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
            if (rank == 0) print *, 'scattering RHS vectors COMPLETE'
    
            ! NOW creating a coefficient distributed matrix A
            if (rank == 0) print *, 'Creating matrix for linear problem  ... '
            call MatCreate(PETSC_COMM_WORLD, A, ierr)
            call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 4 * mx * my, 2 * mx * my, ierr)
            call MatSetFromOptions(A, ierr) 
            call MatSetUp(A, ierr)
    
            if (rank == 0) print *, 'Matrix Memory Allocated', 4 * mx * my,'x', 2 * mx * my
    
            ! Get ownership ranges for rows
            call MatGetOwnershipRange(A, Istart, Iend, ierr)
            call MatGetOwnershipRangeColumn(A, Jstart, Jend, ierr)
    
            !print *, 'rank ', rank, 'rows', Istart, Iend, 'cols', Jstart, Jend
    
            ! Assemble the matrix A and vector b
            do Ii = Istart, Iend - 1
                i = mod(Ii, mx)
                j = Ii / mx
                
                if (Ii < 2 * mx * my) then
                    if (i == 0) then
                        !!! Setting gradient coefficients
                        rightCoeff = -1/dxRightN(i,j)   !! opposite because -dx in matrix
                        centerCoeff = 1/dxRightN(i,j)
                        leftCoeff = 0
                    elseif (i == mx-1) then
                        rightCoeff = 0
                        centerCoeff = -1/dxLeftN(i,j)
                        leftCoeff =    1/dxLeftN(i,j)
                    else
                        rightCoeff = -1/(dxRightN(i,j) + dxLeftN(i,j))
                        centerCoeff = 0 
                        leftCoeff =  1/(dxRightN(i,j) + dxLeftN(i,j))
                    endif
    
                    if (j == 0 .and. i .ne. 0) then
                        !!! Setting gradient coefficients
                        topCoeff = 1/dyTopN(i,j)
                        centerCoeff = -1/dyTopN(i,j)
                        bottomCoeff = 0
                    elseif (j == my-1 .and. i .ne. mx -1 ) then
                        topCoeff = 0 
                        centerCoeff = 1/dyBottomN(i,j)
                        bottomCoeff = -1/dyBottomN(i,j)
                    else
                        topCoeff = 1/(dyTopN(i,j)+dyBottomN(i,j))
                        centerCoeff = 0
                        bottomCoeff = -1/(dyTopN(i,j)+dyBottomN(i,j))
                    endif
    
                
                    if (Ii < mx *my) then
                        !!! Setting gradient coefficients for RHS u
                        ! if ((Ii + (mx + my) + mx) > (2 *(mx + my)) )then
                        !     print *, 'row ',Ii, 'column', Ii + (mx + my) + mx
                        !     stop 
                        ! endif
    
                        colIndex = Ii - 1
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if (i > 0 .and. inColRange )      call MatSetValue(A, Ii, colIndex , leftCoeff   , ADD_VALUES, ierr)
    
                        colIndex = Ii + 1
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if (i < mx -1  .and. inColRange )  call MatSetValue(A, Ii, colIndex , rightCoeff  , ADD_VALUES, ierr)
    
                        colIndex = Ii + (mx + my) - mx
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if (j > 0 .and. inColRange )      call MatSetValue(A, Ii, colIndex, bottomCoeff , ADD_VALUES, ierr) 
                        
                        colIndex = Ii + (mx + my) + mx
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if (j< my -1 .and. inColRange )  call MatSetValue(A, Ii, colIndex, topCoeff    , ADD_VALUES, ierr)
    
                        colIndex = Ii
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if ((i == 0 .or. i == mx -1) .and. inColRange) call MatSetValue(A, Ii, colIndex, centerCoeff , ADD_VALUES, ierr)
    
                        colIndex = Ii + (mx + my)
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if ((j == 0 .or. j == my -1) .and. inColRange)  call MatSetValue(A, Ii, colIndex, centerCoeff , ADD_VALUES, ierr)
    
                    else
                        ! if ((Ii + (mx + my) + mx) > (2 *(mx + my)) )then
                        !     print *, 'row ',Ii, 'column', Ii + (mx + my) + mx
                        !     stop 
                        ! endif
    
                        !!! Setting gradient coefficients for RHS v
                        colIndex = Ii + (mx + my)- 1
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if (i > 0 .and. inColRange )      call MatSetValue(A, Ii, colIndex , leftCoeff   , ADD_VALUES, ierr)
    
                        colIndex = Ii + (mx + my)+ 1
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if (i < mx -1 .and. inColRange )  call MatSetValue(A, Ii, colIndex , rightCoeff  , ADD_VALUES, ierr)
    
                        colIndex = Ii - mx 
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if (j > 0 .and. inColRange )      call MatSetValue(A, Ii, colIndex          , bottomCoeff , ADD_VALUES, ierr)
    
                        colIndex = Ii + mx
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if ( j< my -1 .and. inColRange )  call MatSetValue(A, Ii, colIndex           , topCoeff    , ADD_VALUES, ierr)
    
                        colIndex = Ii
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if ((i == 0 .or. i == mx -1) .and. inColRange ) call MatSetValue(A, Ii, colIndex, centerCoeff , ADD_VALUES, ierr)
    
                        colIndex = Ii + (mx + my)
                        inColRange = colIndex .GE. Jstart  .and. colIndex < Jend
                        if ((j == 0 .or. j == my -1) .and. inColRange ) call MatSetValue(A, Ii, colIndex, centerCoeff , ADD_VALUES, ierr)
    
                    endif
                    
                    
                else
                    !!!! setting laplacian coefficients
                    if (Ii < 3*mx * my) then
                        !!!! Setting Poission problem for phi
                        colIndex = j * mx + i  ! laplacian coeff for psi is zero
                    else
                        colIndex = (mx + my ) + j * mx + i  ! laplacian coeff for phi is zero
                    endif
                    ! Interior grid points: central difference stencil
                    if ((i > 0) .and. (i < mx - 1) .and. (j > 0) .and. (j < my - 1)) then
                        !!! setting five coefficients for each row based on Finite Volume Descretization in the center points
                        rightCoeff = dyRightE(i,j)/dxRightN(i,j)
                        leftCoeff = dyLeftE(i,j)/dxLeftN(i,j)
                        topCoeff = dxTopE(i,j)/dyTopN(i,j)
                        bottomCoeff = dxBottomE(i,j)/dyBottomN(i,j)
                        centerCoeff = -dyLeftE(i,j)/dxLeftN(i,j) - dyRightE(i,j)/dxRightN(i,j) - dxBottomE(i,j)/dyBottomN(i,j) - dxTopE(i,j)/dyBottomN(i,j)
    
                        
                        colIndex2 = colIndex  - 1
                        inColRange = colIndex2 >= Jstart .and. colIndex2 < Jend 
                        if (inColRange) call MatSetValue(A, Ii, colIndex2, leftCoeff   , ADD_VALUES, ierr)
    
                        colIndex2 = colIndex  + 1
                        inColRange = colIndex2 >= Jstart .and. colIndex2 < Jend 
                        if (inColRange) call MatSetValue(A, Ii,  colIndex2, rightCoeff  , ADD_VALUES, ierr)
    
                        colIndex2 = colIndex  - mx
                        inColRange = colIndex2 >= Jstart .and. colIndex2 < Jend 
                        if (inColRange) call MatSetValue(A, Ii, colIndex2, bottomCoeff , ADD_VALUES, ierr) 
    
                        colIndex2 = colIndex  + mx
                        inColRange = colIndex2 >= Jstart .and. colIndex2 < Jend 
                        if (inColRange) call MatSetValue(A, Ii, colIndex2, topCoeff    , ADD_VALUES, ierr)
    
                        colIndex2 = colIndex
                        inColRange = colIndex2 >= Jstart .and. colIndex2 < Jend  
                        if (inColRange) call MatSetValue(A, Ii,  colIndex2, centerCoeff , ADD_VALUES, ierr)
                    endif
                endif
            end do
    
            ! Finalize assembly
            call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    
            if (rank == 0 ) print *, 'Matrix assembly complete'
    
            ! Create the KSP linear solver
            call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
            call KSPSetOperators(ksp, A, A, ierr)
            call KSPSetFromOptions(ksp, ierr)
    
            ! Set solver tolerances: rel_tol, abs_tol, div_tol, and max iterations
            max_iter = 1000
            rel_tol = 1d-4
            abs_tol = 1d-4
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
            call VecScatterBegin(gather, x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)
            call VecScatterEnd(gather, x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)
    
            if (rank == 0) print *, 'gather complete'
    
            if (rank == 0) then
            allocate(collected_xarray(2*mx*my))
            call VecGetArrayReadF90(sol_globalOnZero, collected_xPointer, ierr)
            collected_xarray = collected_xPointer
            phi = reshape(collected_xarray(1:mx*my), (/mx, my/))
            psi = reshape(collected_xarray(mx*my:2*mx*my), (/mx, my/))
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
    
    