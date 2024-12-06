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

    ! subroutine decomposeHelmholtz(uvel, vvel, psi, phi, dxBottomE, dxTopE, dyLeftE, dyRightE, &
    !                             & dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
    !     real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, & ! For right hand side
    !                                                         &  dxBottomE, dxTopE, dyLeftE, dyRightE, & !at cell edges
    !                                                         &  dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea  ! stencil joining nodes

    !     real(kind=real_kind), intent(inout), dimension(:,:) :: psi, phi
    !     real(kind=real_kind), allocatable :: divU(:,:), curlU(:,:)

    !     PetscErrorCode :: ierr
    !     PetscMPIInt :: rank, size

    !     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    !     call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    !     if (rank == 0) call calcHozDivVertCurl(uvel, vvel, dxBottomE, dxTopE, dyLeftE, dyRightE, cellArea, divU, curlU)

    !     if (rank == 0) print *, 'solving poission equation for toroidal part'
    !     call solvepoission(psi, -curlU, dxBottomE, dxTopE, dyLeftE, dyRightE, &
    !                         &dxLeftN, dxRightN, dyBottomN, dyTopN,cellArea)
    !     if (rank == 0) print *, 'solving poission equation for toroidal part complete'

    !     if (rank == 0) print *, 'solving poission equation for poloidal part'
    !     call solvepoission(phi, -divU, dxBottomE, dxTopE, dyLeftE, dyRightE,&
    !                     &  dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
    !     if (rank == 0) print *, 'solving poission equation for poloidal part complete'
    !     if (rank == 0) deallocate(divU, curlU)
    ! end subroutine

    ! subroutine solvepoission(phi, RHS, dxBottomE, dxTopE, dyLeftE, dyRightE, &
    !                             dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
    !     real(kind=real_kind), intent(in), dimension(:,:) :: RHS, &  ! RHS of the poission equation
    !                                                         dxBottomE, dxTopE, dyLeftE, dyRightE, & !at cell edges
    !                                                         dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea  ! stencil joining nodes
    !     real(kind=real_kind), intent(inout), dimension(:,:) :: phi
    !     logical , optional :: isPoloidal
    !     integer :: shapeArr(2)

    !     PetscErrorCode :: ierr
    !     PetscMPIInt :: rank, size
    !     Vec :: x_globalOnZero, x_local, y_globalOnZero, y_local, sol_globalOnZero
    !     VecScatter :: scatterx, scattery, gather
    !     IS :: xis_localVec, xis_globalVec, yis_localVec, yis_globalVec
    !     Mat :: A, AtA
    !     KSP :: ksp
    !     PetscInt :: Istart, Iend, Ii, i, j, mx, my, localSize, globalSize, low, high
    !     PetscScalar :: leftCoeff, rightCoeff, topCoeff, bottomCoeff, centerCoeff
    !     PetscScalar :: val, norm, rel_tol, abs_tol, div_tol
    !     PetscInt :: max_iter

    !     real(kind=8), pointer :: collected_xPointer(:)
    !     real(kind=real_kind), allocatable :: collected_xarray(:)

    !     ! Initialize PETSc
    !     call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    !     call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    !     call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

    !     shapeArr = shape(cellArea)
    !     mx = shapeArr(1)
    !     my = shapeArr(2)

    !     ! create a distributed vector for phi
    !     ! Create the source vector (fully defined on rank 0)

    !     ! The poission equation is solved in the matrix form as [ A ][x] = [y]

    !     call VecCreate(PETSC_COMM_WORLD, x_globalOnZero,  ierr)
    !     call VecSetFromOptions(x_globalOnZero, ierr)
    !     call VecCreate(PETSC_COMM_WORLD, y_globalOnZero, ierr)
    !     call VecSetFromOptions(y_globalOnZero, ierr)

    !     call VecSetSizes(x_globalOnZero, PETSC_DECIDE, mx * my, ierr)
    !     call VecSetSizes(y_globalOnZero, PETSC_DECIDE, mx * my, ierr)

    !     if (rank == 0) then
    !         print *, taskid
    !         print *, 'creating vectors of size', mx,'x', my
            
    !         do Ii = 0, mx * my - 1
    !             i = mod(Ii, mx)
    !             j = Ii / mx
    !             val = phi(i+1,j+1)
    !             ! print *,'i=',i,'j=', j,'Ii=', Ii,'val =', val
    !             call VecSetValue(x_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                
    !             ! if ((i == 0) .and. (j > 0) .and. (j < my - 1)) then  
    !             !     i = i+1
    !             ! else if ((i == mx - 1) .and. (j > 0) .and. (j < my - 1)) then
    !             !     i = i-1
    !             ! else if ((j == 0) .and. (i > 0) .and. (i < mx - 1)) then
    !             !     j = j+1
    !             ! else if ((j == my - 1) .and. (i > 0) .and. (i < mx - 1)) then
    !             !     j = j-1
    !             ! end if

    !             val = RHS(i+1,j+1) * cellArea(i+1,j+1)
    !             if (i == 0 .or. j==0 .or. i == mx-1 .and. j == my -1 ) val = 0.0
    !             !if (i == 0 .or. i == mx -1 .or. j == 0 .or. j == my -1) val = 0.0
    !             call VecSetValue(y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
    !         end do

    !         print *,'i=',i,'j=', j,'Ii=', Ii,'val =', phi(i+1,j+1), RHS(i+1,j+1) * cellArea(i+1,j+1)

    !     end if

    !     call VecAssemblyBegin(x_globalOnZero, ierr)
    !     call VecAssemblyEnd(x_globalOnZero, ierr)

    !     call VecAssemblyBegin(y_globalOnZero, ierr)
    !     call VecAssemblyEnd(y_globalOnZero, ierr)


        

    !     if (rank == 0) print *, 'End vector create'

        

    !     ! Create the distributed destination vector
    !     if (rank == 0) print *, 'creating vector for distributing'
    !     call VecCreate(PETSC_COMM_WORLD, x_local, ierr)
    !     call VecCreate(PETSC_COMM_WORLD, y_local, ierr)
        
    !     call VecSetSizes(x_local, PETSC_DECIDE, mx * my, ierr)
    !     call VecSetSizes(y_local, PETSC_DECIDE, mx * my, ierr)

    !     call VecSetFromOptions(y_local, ierr)
    !     call VecSetFromOptions(x_local, ierr)
    !     if (rank == 0) print *, 'creating vector for distributing COMPLETE'
        
        
    !     call VecGetOwnershipRange(x_local, low, high, ierr) ! /* low, high are global indices */
    !     call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, xis_localVec, ierr)

    !     call VecGetOwnershipRange(y_local, low, high, ierr) ! /* low, high are global indices */
    !     call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, yis_localVec, ierr)
        
    !     ! call VecGetOwnershipRange(vec_global, low, high, ierr) ! /* low, high are global indices */
    !     ! call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, xis_globalVec, ierr)


    !     ! Create scatter context
    !     if (rank == 0) print *, 'scattering vectors'
    !     call VecScatterCreate(x_globalOnZero, xis_localVec, x_local, xis_localVec, scatterx, ierr)

    !     ! Perform the scatter operation
    !     call VecScatterBegin(scatterx, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !     call VecScatterEnd(scatterx, x_globalOnZero, x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !     if (rank == 0) print *, 'scattering X vector COMPLETE'

    !     ! Create scatter context
    !     call VecScatterCreate(y_globalOnZero, yis_localVec, y_local, yis_localVec, scattery, ierr)

    !     ! Perform the scatter operation
    !     call VecScatterBegin(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !     call VecScatterEnd(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !     if (rank == 0) print *, 'scattering vectors COMPLETE'

    !     ! NOW creating a coefficient distributed matrix A
    !     call MatCreate(PETSC_COMM_WORLD, A, ierr)
    !     call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, mx * my, mx * my, ierr)
    !     call MatSetFromOptions(A, ierr) ! commenting this because we are using the default sparse matrix config
    !     call MatSetUp(A, ierr)

    !     ! Get local ownership range of rows
    !     call MatGetOwnershipRange(A, Istart, Iend, ierr)

    !     if (rank == 0) print * , 'Creating Matrix'

    !     ! Assemble the matrix A and vector b
    !     do Ii = Istart, Iend - 1
    !         i = mod(Ii, mx)
    !         j = Ii / mx

    !         ! Interior grid points: central difference stencil
    !         if ((i > 0) .and. (i < mx - 1) .and. (j > 0) .and. (j < my - 1)) then
    !             rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
    !             leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
    !             topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
    !             bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
    !             centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)    - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
    !                     &     -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)

    !             call MatSetValue(A, Ii, Ii - 1, leftCoeff, ADD_VALUES, ierr)
    !             call MatSetValue(A, Ii, Ii + 1, rightCoeff, ADD_VALUES, ierr)
    !             call MatSetValue(A, Ii, Ii - mx, bottomCoeff, ADD_VALUES, ierr) 
    !             call MatSetValue(A, Ii, Ii + mx, topCoeff, ADD_VALUES, ierr)
    !             call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)


    !         else if ((i == 0)) then
    !             centerCoeff = 1
    !             rightCoeff = -1
    !             call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
    !             call MatSetValue(A, Ii, Ii + 1, rightCoeff, ADD_VALUES, ierr)
    !         else if ((i == mx-1)) then
    !             centerCoeff = 1
    !             leftCoeff = -1
    !             call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
    !             call MatSetValue(A, Ii, Ii - 1, leftCoeff, ADD_VALUES, ierr)
    !         else if ((j == my-1)) then
    !             centerCoeff = 1
    !             bottomCoeff = -1
    !             call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
    !             call MatSetValue(A, Ii, Ii -mx, rightCoeff, ADD_VALUES, ierr)
    !         else if ((j == 0)) then
    !             centerCoeff = 1
    !             topCoeff = -1
    !             call MatSetValue(A, Ii, Ii, centerCoeff, ADD_VALUES, ierr)
    !             call MatSetValue(A, Ii, Ii +mx, topCoeff, ADD_VALUES, ierr)
    !         endif
    !         !endif
    !     end do

    !     ! Finalize assembly
    !     call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    !     call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    !     if (rank == 0) print * , "Mat assembly complete"
        
    !     ! Create the KSP linear solver
    !     call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    !     call KSPSetOperators(ksp, A, A, ierr)
    !     call KSPSetFromOptions(ksp, ierr)

    !     ! Set solver tolerances: rel_tol, abs_tol, div_tol, and max iterations
    !     max_iter = 5000
    !     rel_tol = 1d-4
    !     abs_tol = 1d-4
    !     div_tol = 1d5
    !     call KSPSetTolerances(ksp, rel_tol, abs_tol, div_tol,  max_iter, ierr)

    !     ! Solve the linear system A * u = b
    !     call KSPSolve(ksp, y_local, x_local, ierr)

    !     ! Output the solution norm as a summary

    !     call VecNorm(x_local, NORM_2, norm, ierr)
    !     if (rank == 0) then
    !         write(*,*) 'Solution vector 2-norm:', norm
    !     end if

    !     ! Create a scatter context to gather the global vector to rank 0
    !     call VecScatterCreateToZero(x_local, gather, sol_globalOnZero, ierr)
    !     if (rank == 0) print *, 'gather created'
    !     call VecScatterBegin(gather, x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)
    !     if (rank == 0) print *, 'gather began'
    !     call VecScatterEnd(gather, x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)

    !     if (rank == 0) print *, 'gather complete'

    !     if (rank == 0) then
    !         allocate(collected_xarray(mx*my))
    !         call VecGetArrayReadF90(sol_globalOnZero, collected_xPointer, ierr)
    !         collected_xarray = collected_xPointer
    !         phi = reshape(collected_xarray, (/mx, my/))
    !         deallocate(collected_xarray)
    !     endif

    !     ! Clean up
    !     call KSPDestroy(ksp, ierr)
    !     call VecDestroy(x_local, ierr)
    !     call VecDestroy(x_globalOnZero, ierr)
    !     call VecDestroy(y_local, ierr)
    !     call VecDestroy(y_globalOnZero, ierr)
    !     call VecDestroy(sol_globalOnZero, ierr)
    !     call MatDestroy(A, ierr)
    !     call VecScatterDestroy(scatterx, ierr)
    !     call VecScatterDestroy(scattery, ierr)
    !     call VecScatterDestroy(gather, ierr)
    !     call ISDestroy(xis_localVec, ierr)
    !     call ISDestroy(xis_globalVec, ierr)
    !     call ISDestroy(yis_localVec, ierr)
    !     call ISDestroy(yis_globalVec, ierr)
    !     call PetscFinalize(ierr)
    ! end subroutine


    subroutine decomposeHelmholtz_2(uvel, vvel, psi, phi, dxBottomE, dxTopE, dyLeftE, dyRightE, &
        & dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea)
        real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, & ! For right hand side
                                &  dxBottomE, dxTopE, dyLeftE, dyRightE, & !at cell edges
                                &  dxLeftN, dxRightN, dyBottomN, dyTopN, cellArea  ! stencil joining nodes

        real(kind=real_kind), intent(inout), dimension(:,:) :: psi, phi
        !real(kind=real_kind), allocatable :: divU(:,:), curlU(:,:)

        PetscErrorCode :: ierr
        PetscMPIInt :: rank, size

        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

        if (rank == 0) print *, 'solving poission equation'
        call solvepoissionBig2(psi, phi, uvel, vvel, dxBottomE, dxTopE, dyLeftE, dyRightE, &
                                &dxLeftN, dxRightN, dyBottomN, dyTopN,cellArea)
        if (rank == 0) print *, 'solving poission equation complete'

        !if (rank == 0) deallocate(divU, curlU)
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
        Vec :: x_globalOnZero, x_local, y_globalOnZero, y_local, sol_globalOnZero, Ay
        VecScatter :: scatterx, scattery, gather
        IS :: xis_localVec, xis_globalVec, yis_localVec, yis_globalVec
        Mat :: A, N
        KSP :: ksp
        PetscInt :: Istart, Iend, Jstart, Jend, Ii, Jj, i, j, colIndex, &
                    colIndex2, mx, my, localSize, globalSize, low, high
        PetscScalar :: leftCoeff, rightCoeff, topCoeff, bottomCoeff, centerCoeff, centerXCoeff, centerYCoeff
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
                !j = mod(j,my)
        
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
                
                if (i .ne. 0 .or. i .ne. mx -1 .or. j .ne. 0 .or. j .ne. my -1) then
                    val = -divU(i+1, j+1)  * cellArea(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii + 2*(mx * my), val, INSERT_VALUES, ierr)
                    val = -curlU(i+1, j+1) * cellArea(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii + 3*(mx * my), val, INSERT_VALUES, ierr)
                endif

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
        !call MatSetFromOptions(A, ierr) 
        call MatSetType(A, MATAIJ, ierr) 
        call MatSetUp(A, ierr)

        if (rank == 0) print *, 'Matrix Memory Allocated', 4 * mx * my,'x', 2 * mx * my

        ! Get ownership ranges for rows
        call MatGetOwnershipRange(A, Istart, Iend, ierr)
        call MatGetOwnershipRangeColumn(A, Jstart, Jend, ierr)

        !print *, 'rank ', rank, 'rows', Istart, Iend, 'cols', Jstart, Jend

        ! Assemble the matrix A 
        do Ii = Istart, Iend - 1
            i = mod(Ii, mx)
            j = Ii / mx
            j = mod(j, my)
            
            if (Ii < 2 * mx * my) then
                if (i == 0) then
                    !!! Setting gradient coefficients
                    rightCoeff = 1/dxRightN(i+1,j+1)  
                    centerXCoeff = -1/dxRightN(i+1,j+1)
                    leftCoeff = 0
                elseif (i == mx-1) then
                    rightCoeff = 0
                    centerXCoeff = 1/dxLeftN(i+1,j+1)
                    leftCoeff =   -1/dxLeftN(i+1,j+1)
                else
                    rightCoeff = 1/(dxRightN(i+1,j+1) + dxLeftN(i+1,j+1))
                    centerXCoeff = 0 
                    leftCoeff =  -1/(dxRightN(i+1,j+1) + dxLeftN(i+1,j+1))
                endif

                if (j == 0 ) then
                    !!! Setting gradient coefficients
                    topCoeff = 1/dyTopN(i+1,j+1)
                    centerYCoeff = -1/dyTopN(i+1,j+1)
                    bottomCoeff = 0
                elseif (j == my-1 ) then
                    topCoeff = 0 
                    centerYCoeff = 1/dyBottomN(i+1,j+1)
                    bottomCoeff = -1/dyBottomN(i+1,j+1)
                else
                    topCoeff = 1/(dyTopN(i+1,j+1)+dyBottomN(i+1,j+1))
                    centerYCoeff = 0
                    bottomCoeff = -1/(dyTopN(i+1,j+1)+dyBottomN(i+1,j+1))
                endif

            
                if (Ii < mx *my) then
                    !!! Setting gradient coefficients for RHS u
                    ! if ((Ii + (mx + my) + mx) > (2 *(mx + my)) )then
                    !     print *, 'row ',Ii, 'column', Ii + (mx + my) + mx
                    !     stop 
                    ! endif

                    if (i > 0 )  then
                        colIndex = (mx * j + i) - 1
                        if ( colIndex >= 2 *mx *my ) stop 'Here 01'
                        call MatSetValue(A, Ii, colIndex , -leftCoeff   , ADD_VALUES, ierr)
                    endif

                    if (i < mx -1 )  then
                        colIndex = (mx * j + i) + 1
                        if ( colIndex >= 2 * mx *my ) stop 'Here 02'
                        call MatSetValue(A, Ii, colIndex , -rightCoeff  , ADD_VALUES, ierr)
                    endif

                    colIndex = (mx * j + i)
                    if ( colIndex >= 2 *mx *my ) stop 'Here 03'
                    call MatSetValue(A, Ii, colIndex, -centerXCoeff , ADD_VALUES, ierr)


                    if (j > 0 )   then
                        colIndex = (mx * j + i + mx * my) - mx
                        if ( colIndex >= 2*mx *my ) stop 'Here 04'
                        call MatSetValue(A, Ii, colIndex, bottomCoeff , ADD_VALUES, ierr) 
                    endif
                    
                    if (j< my -1 )   then
                        colIndex = (mx * j + i + mx * my) + mx
                        if ( colIndex >= 2*mx *my ) stop 'Here 05'
                        call MatSetValue(A, Ii, colIndex, topCoeff    , ADD_VALUES, ierr)
                    endif

                    colIndex = (mx * j + i + mx * my)
                    if ( colIndex >= 2*mx *my ) stop 'Here 06'
                    call MatSetValue(A, Ii, colIndex, centerYCoeff , ADD_VALUES, ierr)

                else
                    ! if ((Ii + (mx + my) + mx) > (2 *(mx + my)) )then
                    !     print *, 'row ',Ii, 'column', Ii + (mx + my) + mx
                    !     stop 
                    ! endif

                    !!! Setting gradient coefficients for RHS u
                    if (j > 0 )    then
                        colIndex = (mx * j + i) - mx 
                        if ( colIndex >= 2*mx *my ) stop 'Here 07'
                        call MatSetValue(A, Ii, colIndex , -bottomCoeff , ADD_VALUES, ierr)
                    endif

                    if ( j< my -1 )  then
                        colIndex = (mx * j + i) + mx
                        if ( colIndex >= 2*mx *my ) stop 'Here 08'
                        call MatSetValue(A, Ii, colIndex , -topCoeff    , ADD_VALUES, ierr)
                    endif

                    colIndex = (mx * j + i)
                    if ( colIndex >= 2*mx *my ) stop 'Here 09'
                    call MatSetValue(A, Ii, colIndex, -centerYCoeff , ADD_VALUES, ierr)

                    
                    !!! Setting gradient coefficients for RHS v
                    if (i > 0 )     then
                        colIndex = (mx * j + i) - 1
                        if ( colIndex >= 2*mx *my ) stop 'Here 010'
                        call MatSetValue(A, Ii, colIndex , -leftCoeff   , ADD_VALUES, ierr)
                    end if

                    if (i < mx -1)  then
                        colIndex = (mx * j + i) + 1
                        if ( colIndex >= 2*mx *my ) stop 'Here 011'
                        call MatSetValue(A, Ii, colIndex , -rightCoeff  , ADD_VALUES, ierr)
                    endif

                    colIndex = (mx * j + i) 
                    if ( colIndex >= 2*mx *my ) stop 'Here 012'
                    call MatSetValue(A, Ii, colIndex, -centerXCoeff , ADD_VALUES, ierr)

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
                    rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
                    leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
                    topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
                    bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
                    centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1) - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
                                & -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)

                    
                    colIndex2 = colIndex  - 1
                    if ( colIndex2 >= 2*mx *my ) stop 'Here 13'
                    call MatSetValue(A, Ii, colIndex2, leftCoeff   , ADD_VALUES, ierr)

                    colIndex2 = colIndex  + 1
                    if ( colIndex2 >= 2*mx *my ) stop 'Here 14'
                    call MatSetValue(A, Ii,  colIndex2, rightCoeff  , ADD_VALUES, ierr)

                    colIndex2 = colIndex  - mx
                    if ( colIndex2 >= 2*mx *my ) stop 'Here 15'
                    call MatSetValue(A, Ii, colIndex2, bottomCoeff , ADD_VALUES, ierr) 

                    colIndex2 = colIndex  + mx
                    if ( colIndex2 >= 2*mx *my ) stop 'Here 16'
                    call MatSetValue(A, Ii, colIndex2, topCoeff    , ADD_VALUES, ierr)

                    colIndex2 = colIndex
                    if ( colIndex2 >= 2*mx *my ) stop 'Here 17'
                    call MatSetValue(A, Ii,  colIndex2, centerCoeff , ADD_VALUES, ierr)

                else if ((i == 0) .and. (j > 0) .and. (j < my - 1)) then
                    ! left boundary on physical lat, lon, grid
                    i = i+1
                    rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
                    leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
                    topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
                    bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
                    centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1) - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
                                & -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)
                    
                    call MatSetValue(A, Ii, colIndex +1 - 1, leftCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex +1 + 1, rightCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex +1 - mx, bottomCoeff, ADD_VALUES, ierr) 
                    call MatSetValue(A, Ii, colIndex +1 + mx, topCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex +1, centerCoeff, ADD_VALUES, ierr)


                else if ((i == mx - 1) .and. (j > 0) .and. (j < my - 1)) then
                    i = i-1
                    rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
                    leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
                    topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
                    bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
                    centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1) - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
                                & -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)

                    call MatSetValue(A, Ii, colIndex -1 - 1, leftCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex -1 + 1, rightCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex -1 - mx, bottomCoeff, ADD_VALUES, ierr) 
                    call MatSetValue(A, Ii, colIndex -1 + mx, topCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex -1, centerCoeff, ADD_VALUES, ierr)



                else if ((j == 0) .and. (i > 0) .and. (i < mx - 1)) then
                    ! bottom boundary physical lat, lon, grid
                    j = j+1
                    rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
                    leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
                    topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
                    bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
                    centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1) - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
                                & -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)

                    call MatSetValue(A, Ii, colIndex +mx - 1, leftCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex +mx + 1, rightCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex +mx - mx, bottomCoeff, ADD_VALUES, ierr) 
                    call MatSetValue(A, Ii, colIndex +mx + mx, topCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex +mx, centerCoeff, ADD_VALUES, ierr)


                else if ((j == my - 1) .and. (i > 0) .and. (i < mx - 1)) then
                    !top boundary last column
                    j = j-1
                    rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
                    leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
                    topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
                    bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
                    centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1) - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
                                & -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)

                    call MatSetValue(A, Ii, colIndex -mx - 1, leftCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex -mx + 1, rightCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex -mx - mx, bottomCoeff, ADD_VALUES, ierr) 
                    call MatSetValue(A, Ii, colIndex -mx + mx, topCoeff, ADD_VALUES, ierr)
                    call MatSetValue(A, Ii, colIndex -mx, centerCoeff, ADD_VALUES, ierr)
                endif
            endif
        end do

        ! Finalize assembly
        call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

        !call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

        if (rank == 0 ) print *, 'Matrix assembly complete'

        ! call MatCreateNormalHermitian(A, N, ierr)
        call VecDuplicate(x_local, Ay, ierr) ! just to have same config
        !call MatMultHermitianTranspose(A, y_local, Ay, ierr)
        call MatCreateNormal(A, N, ierr)

        call MatMultTranspose(A, y_local, Ay, ierr)

        ! Create the KSP linear solver
        call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
        call KSPSetOperators(ksp, N, N, ierr)
        
        ! Use LSQR solver for rectangular system (least squares problem)
        ! call KSPSetType(ksp, KSPLSQR, ierr)
        call KSPSetFromOptions(ksp, ierr)

        ! Set solver tolerances: rel_tol, abs_tol, div_tol, and max iterations
        max_iter = 25000
        rel_tol = 1d-10
        abs_tol = 1d-10
        div_tol = 1d5
        call KSPSetTolerances(ksp, rel_tol, abs_tol, div_tol,  max_iter, ierr)

        ! Solve the linear system A * u = b
        call KSPSolve(ksp, Ay, x_local, ierr)

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
            print *, 'array_allocated'
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
        call VecDestroy(Ay, ierr)
        call MatDestroy(A, ierr)
        call MatDestroy(N, ierr)
        call VecScatterDestroy(scatterx, ierr)
        call VecScatterDestroy(scattery, ierr)
        call VecScatterDestroy(gather, ierr)
        call ISDestroy(xis_localVec, ierr)
        call ISDestroy(xis_globalVec, ierr)
        call ISDestroy(yis_localVec, ierr)
        call ISDestroy(yis_globalVec, ierr)
        call PetscFinalize(ierr)
    end subroutine


    subroutine solvepoissionBig2(phi, psi, uvel, vvel, dxBottomE, dxTopE, dyLeftE, dyRightE, &
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
        !         
        !         | laplacian        0      | |  phi  | | -div   | 
        !         |     0        laplacian  | |  psi  | | -curl  |
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
        call VecSetSizes(y_globalOnZero, PETSC_DECIDE, 2*mx * my, ierr)

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
                if (i == 0 .or. i == mx -1 .or. j == 0 .or. j == my -1 ) then
                    val = uvel(i+1, j+1)  + vvel(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                else
                    val = -divU(i+1, j+1)  * cellArea(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                    val = -curlU(i+1, j+1) * cellArea(i+1, j+1)
                    call VecSetValue(y_globalOnZero, Ii + (mx * my), val, INSERT_VALUES, ierr)
                endif

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
        call VecSetSizes(y_local, PETSC_DECIDE, 2*mx * my, ierr)

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
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 2 * mx * my, 2 * mx * my, ierr)
        call MatSetFromOptions(A, ierr) 
        !call MatSetType(A, MATAIJ, ierr) 
        call MatSetUp(A, ierr)

        if (rank == 0) print *, 'Matrix Memory Allocated', 2 * mx * my,'x', 2 * mx * my

        ! Get ownership ranges for rows
        call MatGetOwnershipRange(A, Istart, Iend, ierr)
        call MatGetOwnershipRangeColumn(A, Jstart, Jend, ierr)

        !print *, 'rank ', rank, 'rows', Istart, Iend, 'cols', Jstart, Jend

        ! Assemble the matrix A and vector b
        do Ii = Istart, Iend - 1
            i = mod(Ii, mx)
            j = Ii / mx
            j = mod(j,my)
            
            if (i > 0 .and. i < mx -1 .and. j > 0 .and. j < my -1) then 
                !!!! setting laplacian coefficients for interior grid points
                if (Ii < mx * my) then
                    !!!! Setting Poission problem for phi
                    colIndex = j * mx + i  ! laplacian coeff for psi is zero
                else
                    colIndex = (j * mx + i + mx * my)  ! laplacian coeff for phi is zero
                endif

                !!! setting five coefficients for each row based on Finite Volume Descretization in the center points
                rightCoeff = dyRightE(i+1,j+1)/dxRightN(i+1,j+1)
                leftCoeff = dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1)
                topCoeff = dxTopE(i+1,j+1)/dyTopN(i+1,j+1)
                bottomCoeff = dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1)
                centerCoeff = -dyLeftE(i+1,j+1)/dxLeftN(i+1,j+1) - dyRightE(i+1,j+1)/dxRightN(i+1,j+1) &
                            & -dxBottomE(i+1,j+1)/dyBottomN(i+1,j+1) - dxTopE(i+1,j+1)/dyBottomN(i+1,j+1)

                colIndex2 = colIndex  - 1
                if (colIndex2 > 2 * mx * my ) stop 'Here 01'
                call MatSetValue(A, Ii, colIndex2, leftCoeff, ADD_VALUES, ierr)

                colIndex2 = colIndex  + 1
                if (colIndex2 > 2 * mx * my ) stop 'Here 02'
                call MatSetValue(A, Ii,  colIndex2, rightCoeff  , ADD_VALUES, ierr)

                colIndex2 = colIndex  - mx
                if (colIndex2 > 2 * mx * my ) stop 'Here 03'
                call MatSetValue(A, Ii, colIndex2, bottomCoeff , ADD_VALUES, ierr) 

                colIndex2 = colIndex  + mx
                if (colIndex2 > 2 * mx * my ) stop 'Here 04'
                call MatSetValue(A, Ii, colIndex2, topCoeff    , ADD_VALUES, ierr)

                colIndex2 = colIndex
                if (colIndex2 > 2 * mx * my ) stop 'Here 05'
                call MatSetValue(A, Ii,  colIndex2, centerCoeff , ADD_VALUES, ierr)
            
            else
                if ( j == 0) then
                    !!! Setting dy gradient coefficients phi
                    topCoeff = -1/dyTopN(i+1,j+1)   ! negative because -dy for phi
                    centerCoeff = 1/dyTopN(i+1,j+1)

                    colIndex = Ii + mx
                    if (colIndex > 2 * mx * my ) stop 'Here 06'
                    call MatSetValue(A, Ii, colIndex , topCoeff, ADD_VALUES, ierr)

                    colIndex = Ii
                    if (colIndex > 2 * mx * my ) stop 'Here 07'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)
                    
                    !!! Setting dy gradient coefficients psi
                    topCoeff = 1/dyTopN(i+1,j+1)
                    centerCoeff = -1/dyTopN(i+1,j+1)

                    colIndex = (j * mx + i + mx * my) + mx
                    if (colIndex > 2 * mx * my ) stop 'Here 08'
                    call MatSetValue(A, Ii, colIndex , topCoeff, ADD_VALUES, ierr)

                    colIndex = (j * mx + i + mx * my)
                    if (colIndex > 2 * mx * my ) stop 'Here 09'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)

                elseif ( j == my - 1) then
                    !!! Setting dy gradient coefficients phi
                    centerCoeff = -1/dyTopN(i+1,j+1) ! negative because -dy for phi
                    bottomCoeff = 1/dyTopN(i+1,j+1)

                    colIndex = Ii - mx
                    if (colIndex > 2 * mx * my ) stop 'Here 10'
                    call MatSetValue(A, Ii, colIndex , bottomCoeff, ADD_VALUES, ierr)

                    colIndex = Ii
                    if (colIndex > 2 * mx * my ) stop 'Here 11'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)

                    !!! Setting dy gradient coefficients psi
                    centerCoeff = 1/dyTopN(i+1,j+1)
                    bottomCoeff = -1/dyTopN(i+1,j+1)

                    colIndex = (j * mx + i + mx * my) - mx
                    if (colIndex > 2 * mx * my ) stop 'Here 12'
                    call MatSetValue(A, Ii, colIndex , bottomCoeff, ADD_VALUES, ierr)

                    colIndex = (j * mx + i + mx * my)
                    if (colIndex > 2 * mx * my ) stop 'Here 13'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)

                else
                    !!! Setting dy gradient coefficients phi
                    topCoeff = -1/(dyTopN(i+1,j+1) + dyBottomN(i+1,j+1))! negative because -dy for phi
                    bottomCoeff = 1/(dyTopN(i+1,j+1) + dyBottomN(i+1,j+1))

                    colIndex = Ii + mx
                    if (colIndex > 2 * mx * my ) stop 'Here 14'
                    call MatSetValue(A, Ii, colIndex , topCoeff, ADD_VALUES, ierr)

                    colIndex = Ii - mx
                    if (colIndex > 2 * mx * my ) stop 'Here 15'
                    call MatSetValue(A, Ii, colIndex , bottomCoeff, ADD_VALUES, ierr)


                    !!! Setting dy gradient coefficients psi
                    topCoeff = 1/(dyTopN(i+1,j+1) + dyBottomN(i+1,j+1))
                    bottomCoeff = - 1/(dyTopN(i+1,j+1) + dyBottomN(i+1,j+1))

                    colIndex = (j * mx + i + mx * my) + mx
                    if (colIndex > 2 * mx * my ) then
                        print *, 'i, j, Ii', i, j, Ii
                        stop 'Here 16'
                    endif
                    call MatSetValue(A, Ii, colIndex , topCoeff, ADD_VALUES, ierr)

                    colIndex = (j * mx + i + mx * my) - mx
                    if (colIndex > 2 * mx * my ) stop 'Here 17'
                    call MatSetValue(A, Ii, colIndex , bottomCoeff, ADD_VALUES, ierr)
                
                endif

                if (i == 0) then
                    !!! Setting gradient coefficients phi
                    rightCoeff = -1/dxRightN(i+1,j+1)   !! opposite because -dx in matrix for phi
                    centerCoeff = 1/dxRightN(i+1,j+1)

                    colIndex = Ii +1
                    if (colIndex > 2 * mx * my ) stop 'Here 22'
                    call MatSetValue(A, Ii, colIndex , rightCoeff, ADD_VALUES, ierr)

                    colIndex = Ii
                    if (colIndex > 2 * mx * my ) stop 'Here 23'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)

                    !!! Setting gradient coefficients psi
                    rightCoeff = -1/dxRightN(i+1,j+1)   !! opposite because -dx in matrix for psi
                    centerCoeff = 1/dxRightN(i+1,j+1)

                    colIndex = (j * mx + i + mx * my) +1
                    if (colIndex > 2 * mx * my ) stop 'Here 24'
                    call MatSetValue(A, Ii, colIndex , rightCoeff, ADD_VALUES, ierr)

                    colIndex = (j * mx + i + mx * my)
                    if (colIndex > 2 * mx * my ) stop 'Here 25'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)


                elseif (i == mx-1) then
                    !!! Setting gradient coefficients phi
                    leftCoeff = 1/dxLeftN(i+1,j+1)   !! opposite because -dx in matrix for phi
                    centerCoeff = -1/dxLeftN(i+1,j+1)

                    colIndex = Ii -1
                    if (colIndex > 2 * mx * my ) stop 'Here 26'
                    call MatSetValue(A, Ii, colIndex , leftCoeff, ADD_VALUES, ierr)

                    colIndex = Ii
                    if (colIndex > 2 * mx * my ) stop 'Here 27'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)

                    !!! Setting gradient coefficients psi
                    leftCoeff = 1/dxLeftN(i+1,j+1)   !! opposite because -dx in matrix for psi
                    centerCoeff = -1/dxLeftN(i+1,j+1)

                    colIndex = (j * mx + i + mx * my) -1
                    if (colIndex > 2 * mx * my ) then
                        print *, 'i, j, Ii', i, j, Ii
                        stop 'Here 28'
                    end if
                    call MatSetValue(A, Ii, colIndex , leftCoeff, ADD_VALUES, ierr)

                    colIndex = (j * mx + i + mx * my)
                    if (colIndex > 2 * mx * my ) stop 'Here 29'
                    call MatSetValue(A, Ii, colIndex , centerCoeff, ADD_VALUES, ierr)
                
                else
                    !!! Setting gradient coefficients phi
                    leftCoeff = 1/(dxLeftN(i+1,j+1)  + dxRightN(i+1,j+1)) !! opposite because -dx in matrix for phi
                    rightCoeff = -1/(dxLeftN(i+1,j+1)  + dxRightN(i+1,j+1))

                    colIndex = Ii -1
                    if (colIndex > 2 * mx * my ) stop 'Here 30'
                    call MatSetValue(A, Ii, colIndex , leftCoeff, ADD_VALUES, ierr)

                    colIndex = Ii +1
                    if (colIndex > 2 * mx * my ) stop 'Here 31'
                    call MatSetValue(A, Ii, colIndex , rightCoeff, ADD_VALUES, ierr)

                    !!! Setting gradient coefficients psi
                    leftCoeff = 1/(dxLeftN(i+1,j+1)  + dxRightN(i+1,j+1)) !! opposite because -dx in matrix for psi
                    rightCoeff = -1/(dxLeftN(i+1,j+1)  + dxRightN(i+1,j+1))

                    colIndex = (j * mx + i + mx * my) -1
                    if (colIndex > 2 * mx * my ) stop 'Here 32'
                    call MatSetValue(A, Ii, colIndex , leftCoeff, ADD_VALUES, ierr)

                    colIndex = (j * mx + i + mx * my) +1
                    if (colIndex > 2 * mx * my ) stop 'Here 33'
                    call MatSetValue(A, Ii, colIndex , rightCoeff, ADD_VALUES, ierr)
                endif
            endif

        end do

        ! Finalize assembly
        call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

        !call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

        if (rank == 0 ) print *, 'Matrix assembly complete'

        ! Create the KSP linear solver
        call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
        call KSPSetOperators(ksp, A, A, ierr)
        
        ! Use LSQR solver for rectangular system (least squares problem)
        call KSPSetType(ksp, KSPLSQR, ierr)
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
    
    