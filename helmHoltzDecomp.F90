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

    type :: linearSystemMat
        Mat :: A
        Vec :: x_globalOnZero, x_local, y_globalOnZero, y_local
        integer :: nx, ny
    contains 
            procedure :: setMat => set_mat
            procedure :: delMat => del_mat
            procedure :: setLHS => set_lhs
            procedure :: setRHS => set_rhs
            procedure :: resetLHS => reset_lhs
            procedure :: resetRHS => reset_rhs
            procedure :: getSol => get_lhsOnZero
            procedure :: solve => solve_system
    end type

    PetscMPIInt :: rank, size_nprocs
    PetscErrorCode :: ierr
        

    class(linearSystemMat), allocatable :: multiGridMats(:)

    contains

    subroutine initPETSC()
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, size_nprocs, ierr)
    end subroutine

    subroutine finalizePETSC()
        call PetscFinalize(ierr)
    end subroutine

    subroutine set_mat(self, nx, ny, dx, dy)
        class(linearSystemMat), intent(inout) :: self
        integer, intent(in) :: nx, ny
        real(kind=real_kind), intent(in) :: dx(nx, ny), dy(nx, ny)

        integer, allocatable :: relPosX(:), relPosY(:)
        real(kind=real_kind), allocatable :: fdXCoeffs(:), fdYCoeffs(:)
        integer :: schemeX, schemeY, ncoeffsX, ncoeffsY, count

        PetscInt :: Istart, Iend, Ii, i, j, mx, my, &
                    colIndex
        PetscScalar :: val

        mx = nx
        my = ny

        call MatCreate(PETSC_COMM_WORLD, self%A, ierr)
        call MatSetSizes(self%A, PETSC_DECIDE, PETSC_DECIDE, 4 * mx * my, 2 * mx * my, ierr)
        
        call MatSetType(self%A, MATAIJ, ierr) 
        call MatSetUp(self%A, ierr)

        if (rank == 0) print *, 'Matrix Memory Allocated', 4 * mx * my,'x', 2 * mx * my

        ! Get ownership ranges for rows
        call MatGetOwnershipRange(self%A, Istart, Iend, ierr)
        
        ! Createing the LHS matrix of the following system
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SYSTEM OF EQUATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !          
        !  |    -dx          dy      |              |   u      |
        !  |    -dy         -dx      | |  phi  | =  |   v      |
        !  | laplacian        0      | |  psi  |    | -div  V  | 
        !  |     0        laplacian  |              | -curl V  |
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        

        do Ii = Istart, Iend - 1
            i = mod(Ii, mx)
            j = Ii / mx
            j = mod(j, my)

            if (i < 3) then
                schemeX = 1
            else if (i > mx-4) then
                schemeX = -1
            else
                schemeX = 0
            endif

            if (j < 3) then
                schemeY = 1
            else if (j > my-4) then
                schemeY = -1
            else
                schemeY = 0
            endif

            if (Ii < 2*(mx*my)) then
                call getFDcoefficients(fdXCoeffs, relPosX, derOrder=1, accuracyOrder=6, scheme=schemeX)
                call getFDcoefficients(fdYCoeffs, relPosY, derOrder=1, accuracyOrder=6, scheme=schemeY)

                ncoeffsX = SIZE(fdXCoeffs)
                ncoeffsY = SIZE(fdYCoeffs)                
                if (Ii < (mx*my)) then
                    !!! SETTING THE FIRST ROW TERMS IN THE MATRIX
                    do count =1, ncoeffsX  !setting -dx terms in the matrix
                        val = -fdXCoeffs(count) * 1.0/dx(i+1, j+1)
                        colIndex = Ii + relPosX(count)
                        call MatSetValue(self%A, Ii, colIndex , val, ADD_VALUES, ierr)
                    enddo

                    do count =1, ncoeffsY  !setting dy terms in the matrix
                        val = fdYCoeffs(count) * 1.0/dy(i+1, j+1)
                        colIndex = Ii + (mx * my) + relPosY(count) * mx
                        call MatSetValue(self%A, Ii, colIndex , val, ADD_VALUES, ierr)
                    enddo
                else
                    !!! SETTING THE SECOND ROW TERMS IN THE MATRIX
                    do count =1, ncoeffsY  !setting -dy terms in the matrix
                        val = -fdYCoeffs(count) * 1.0/dy(i+1, j+1)
                        colIndex = Ii - (mx * my) + relPosY(count) * mx
                        call MatSetValue(self%A, Ii, colIndex , val, ADD_VALUES, ierr)
                    enddo

                    do count =1, ncoeffsX  !setting -dx terms in the matrix
                        val = -fdXCoeffs(count) * 1.0/dx(i+1, j+1)
                        colIndex = Ii + relPosX(count)
                        call MatSetValue(self%A, Ii, colIndex , val, ADD_VALUES, ierr)
                    enddo

                endif
            else
                call getFDcoefficients(fdXCoeffs, relPosX, derOrder=2, accuracyOrder=6, scheme=schemeX)
                call getFDcoefficients(fdYCoeffs, relPosY, derOrder=2, accuracyOrder=6, scheme=schemeY)

                ncoeffsX = size(fdXCoeffs)
                ncoeffsY = size(fdYCoeffs)

                do count =1, ncoeffsX  !setting dx^2 terms in the matrix
                    val = fdXCoeffs(count) * 1.0/dx(i+1, j+1)**2
                    colIndex = Ii - 2*(mx * my) + relPosX(count)
                    call MatSetValue(self%A, Ii, colIndex , val, ADD_VALUES, ierr)
                enddo

                do count =1, ncoeffsY  !setting dy^2 terms in the matrix
                    val = fdYCoeffs(count) * 1.0/dy(i+1, j+1)**2
                    colIndex = Ii - 2*(mx * my) + relPosY(count) * mx
                    call MatSetValue(self%A, Ii, colIndex , val, ADD_VALUES, ierr)
                enddo
            endif
        enddo

        call MatAssemblyBegin(self%A, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(self%A, MAT_FINAL_ASSEMBLY, ierr)

        !call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

        if (rank == 0 ) print *, 'Matrix assembly complete'

    end subroutine

    subroutine del_mat(self)
        class(linearSystemMat), intent(inout) :: self
        call MatDestroy(self%A, ierr)
    end subroutine
    

    subroutine set_rhs(self, nx, ny, centerDx, centerDy, uvel, vvel)
        class(linearSystemMat) :: self
        integer, intent(in) :: nx, ny
        real (kind=real_kind), intent(in) :: centerDx(nx, ny), centerDy(nx, ny), &
                                             uvel(:, :), vvel(:, :)

        real (kind=real_kind), allocatable :: divU(:,:), curlU(:,:)

        PetscInt :: mx, my, Ii,  low, high, i, j
        PetscScalar :: val
        IS :: yis_localVec 
        VecScatter :: scattery

        mx = nx
        my = ny

        call VecCreate(PETSC_COMM_WORLD, self%y_globalOnZero, ierr)
        call VecSetFromOptions(self%y_globalOnZero, ierr)
        call VecSetSizes(self%y_globalOnZero, PETSC_DECIDE, 4*mx * my, ierr)

        if (rank == 0) then
            allocate(divU(nx,ny), curlU(nx,ny))
            call calcHozDivVertCurlFD(uvel, vvel, centerDx, centerDy, divU, curlU)
            print *, 'creating RHS vectors grid size size', mx,'x', my

            do Ii = 0, mx * my - 1
                i = mod(Ii, mx)
                j = Ii / mx
                
                !!!!!!!!!!!!!!!!!!!!!!!!!! setting RHS !!!!!!!!!!!!!!!!!!!!!!!
                val = uvel(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                
                val = vvel(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii + (mx * my), val, INSERT_VALUES, ierr)
                
                val = -divU(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii + 2 * (mx * my), val, INSERT_VALUES, ierr)

                val = -curlU(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii + 3 * (mx * my), val, INSERT_VALUES, ierr)
            end do

            deallocate(divU, curlU)

            print *, 'RHS vectors assigned values Set'
        end if

        call VecAssemblyBegin(self%y_globalOnZero, ierr)
        call VecAssemblyEnd(self%y_globalOnZero, ierr)

        ! Create the distributed destination vector
        if (rank == 0) print *, 'creating RHS vector for distributing across processors'

        call VecCreate(PETSC_COMM_WORLD, self%y_local, ierr)
        call VecSetSizes(self%y_local, PETSC_DECIDE, 4*mx * my, ierr)

        call VecSetFromOptions(self%y_local, ierr)
        
        if (rank == 0) print *, 'creating RHS vector for distributing across processors COMPLETE'

        call VecGetOwnershipRange(self%y_local, low, high, ierr) ! /* low, high are global indices */
        call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, yis_localVec, ierr)

        if (rank == 0) print *, 'scattering RHS vector  ... '
        ! Create scatter context
        call VecScatterCreate(self%y_globalOnZero, yis_localVec, self%y_local, yis_localVec, scattery, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scattery, self%y_globalOnZero, self%y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scattery, self%y_globalOnZero, self%y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'scattering RHS vectors COMPLETE'

        call ISDestroy(yis_localVec, ierr)
        call VecScatterDestroy(scattery, ierr)
    end subroutine

    subroutine reset_rhs(self, nx, ny, centerDx, centerDy, uvel, vvel)
        class(linearSystemMat) :: self
        integer, intent(in) :: nx, ny
        real (kind=real_kind), intent(in) :: centerDx(nx, ny), centerDy(nx, ny), &
                                             uvel(:, :), vvel(:, :)

        real (kind=real_kind), allocatable :: divU(:,:), curlU(:,:)

        PetscInt :: mx, my, Ii,  low, high, i, j
        PetscScalar :: val
        IS :: yis_localVec 
        VecScatter :: scattery

        mx = nx
        my = ny

        if (rank == 0) then
            allocate(divU(nx,ny), curlU(nx,ny))
            call calcHozDivVertCurlFD(uvel, vvel, centerDx, centerDy, divU, curlU)
            do Ii = 0, mx * my - 1
                i = mod(Ii, mx)
                j = Ii / mx
                
                !!!!!!!!!!!!!!!!!!!!!!!!!! setting RHS !!!!!!!!!!!!!!!!!!!!!!!
                val = uvel(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                
                val = vvel(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii + (mx * my), val, INSERT_VALUES, ierr)
                
                val = -divU(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii + 2 * (mx * my), val, INSERT_VALUES, ierr)

                val = -curlU(i+1, j+1) 
                call VecSetValue(self%y_globalOnZero, Ii + 3 * (mx * my), val, INSERT_VALUES, ierr)
            end do

            deallocate(divU, curlU)

            print *, 'RHS vectors reassigned values'
        end if

        call VecAssemblyBegin(self%y_globalOnZero, ierr)
        call VecAssemblyEnd(self%y_globalOnZero, ierr)

        call VecGetOwnershipRange(self%y_local, low, high, ierr) ! /* low, high are global indices */
        call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, yis_localVec, ierr)

        if (rank == 0) print *, 'scattering RHS vector  ... '
        ! Create scatter context
        call VecScatterCreate(self%y_globalOnZero, yis_localVec, self%y_local, yis_localVec, scattery, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scattery, self%y_globalOnZero, self%y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scattery, self%y_globalOnZero, self%y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'scattering RHS vectors COMPLETE'

        call ISDestroy(yis_localVec, ierr)
        call VecScatterDestroy(scattery, ierr)
    end subroutine
    

    subroutine set_lhs(self, nx, ny, phi, psi)
        class(linearSystemMat) :: self
        integer, intent(in) :: nx, ny
        real (kind=real_kind), intent(in) :: phi(:, :), psi(:, :)

        PetscInt :: mx, my, Ii,  low, high, i, j
        PetscScalar :: val
        IS :: xis_localVec 
        VecScatter :: scatterx

        mx = nx
        my = ny

        call VecCreate(PETSC_COMM_WORLD, self%x_globalOnZero, ierr)
        call VecSetFromOptions(self%x_globalOnZero, ierr)

        call VecSetSizes(self%x_globalOnZero, PETSC_DECIDE, 2 * mx * my, ierr)

        if (rank == 0) then
            print *, 'creating LHS vectors for grid size', mx,'x', my

            do Ii = 0, mx * my - 1
                i = mod(Ii, mx)
                j = Ii / mx
                
                !!!!!!!!!!!!!!!!!!!!!!!!!! setting LHS !!!!!!!!!!!!!!!!!!!!!!!
                val = phi(i+1, j+1) 
                call VecSetValue(self%x_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                
                val = psi(i+1, j+1) 
                call VecSetValue(self%x_globalOnZero, Ii + (mx * my), val, INSERT_VALUES, ierr)
                
            end do

            print *, 'LHS vectors assigned values Set'
        end if

        call VecAssemblyBegin(self%x_globalOnZero, ierr)
        call VecAssemblyEnd(self%x_globalOnZero, ierr)

        ! Create the distributed destination vector
        if (rank == 0) print *, 'creating LHS vector for distributing across processors'

        call VecCreate(PETSC_COMM_WORLD, self%x_local, ierr)
        call VecSetSizes(self%x_local, PETSC_DECIDE, 2 * mx * my, ierr)

        call VecSetFromOptions(self%x_local, ierr)
        
        if (rank == 0) print *, 'creating LHS vector for distributing across processors COMPLETE'

        call VecGetOwnershipRange(self%x_local, low, high, ierr) ! /* low, high are global indices */
        call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, xis_localVec, ierr)

        if (rank == 0) print *, 'scattering LHS vector  ... '
        ! Create scatter context
        call VecScatterCreate(self%x_globalOnZero, xis_localVec, self%x_local, xis_localVec, scatterx, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scatterx, self%x_globalOnZero, self%x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scatterx, self%x_globalOnZero, self%x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'scattering LHS vectors COMPLETE'

        call ISDestroy(xis_localVec, ierr)
        call VecScatterDestroy(scatterx, ierr)
    end subroutine

    subroutine reset_lhs(self, nx, ny, phi, psi)
        class(linearSystemMat) :: self
        integer, intent(in) :: nx, ny
        real (kind=real_kind), intent(in) :: phi(:, :), psi(:, :)

        PetscInt :: mx, my, Ii,  low, high, i, j
        PetscScalar :: val
        IS :: xis_localVec 
        VecScatter :: scatterx

        mx = nx
        my = ny

        if (rank == 0) then
            do Ii = 0, mx * my - 1
                i = mod(Ii, mx)
                j = Ii / mx
                
                !!!!!!!!!!!!!!!!!!!!!!!!!! setting LHS !!!!!!!!!!!!!!!!!!!!!!!
                val = phi(i+1, j+1) 
                call VecSetValue(self%x_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                
                val = psi(i+1, j+1) 
                call VecSetValue(self%x_globalOnZero, Ii + (mx * my), val, INSERT_VALUES, ierr)
                
            end do

            print *, 'LHS vectors reassigned values'
        end if

        call VecAssemblyBegin(self%x_globalOnZero, ierr)
        call VecAssemblyEnd(self%x_globalOnZero, ierr)

        call VecGetOwnershipRange(self%x_local, low, high, ierr) ! /* low, high are global indices */
        call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, xis_localVec, ierr)

        if (rank == 0) print *, 'scattering LHS vector  ... '
        ! Create scatter context
        call VecScatterCreate(self%x_globalOnZero, xis_localVec, self%x_local, xis_localVec, scatterx, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scatterx, self%x_globalOnZero, self%x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scatterx, self%x_globalOnZero, self%x_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'scattering LHS vectors COMPLETE'

        call ISDestroy(xis_localVec, ierr)
        call VecScatterDestroy(scatterx, ierr)
    end subroutine

    subroutine get_lhsOnZero(self, nx, ny, phi, psi)
        class(linearSystemMat) :: self
        integer, intent(in) :: nx, ny

        Vec :: sol_globalOnZero
        PetscInt :: globalSize, mx, my
        VecScatter :: gather

        real(kind=8), pointer :: collected_xPointer(:)
        real(kind=real_kind), allocatable :: collected_xarray(:)
        real(kind=real_kind), intent(out) :: phi(nx,ny), psi(nx,ny)

        mx = nx
        my = ny

        ! Create a scatter context to gather the global vector to rank 0
        call VecScatterCreateToZero(self%x_local, gather, sol_globalOnZero, ierr)
        call VecScatterBegin(gather, self%x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(gather, self%x_local, sol_globalOnZero, INSERT_VALUES, SCATTER_FORWARD, ierr)

        if (rank == 0) print *, 'gather complete'

        call VecGetSize(sol_globalOnZero, globalSize, ierr)
        
        if (rank == 0) then
            allocate(collected_xarray(2*mx*my))
            print *, 'global size of solution vec', globalSize
            call VecGetArrayReadF90(sol_globalOnZero, collected_xPointer, ierr)
            collected_xarray = collected_xPointer
            phi = reshape(collected_xarray(1:mx*my), (/mx, my/))
            psi = reshape(collected_xarray(mx*my:2*mx*my), (/mx, my/))
            deallocate(collected_xarray)
            nullify(collected_xPointer)
        endif

        call VecDestroy(sol_globalOnZero, ierr)
    end subroutine

    subroutine solve_system(self, maxIter, relTol, absTol, divTol)
        class(linearSystemMat) :: self
        integer :: maxIter
        real(kind=real_kind) :: relTol, absTol, divTol

        PetscScalar :: rel_tol, abs_tol, div_tol, norm
        PetscInt :: max_iter
        Mat :: N
        Vec :: Ay
        KSP :: ksp

        rel_tol = relTol
        abs_tol = absTol
        div_tol = divTol
        max_iter = maxIter
        
        call MatCreateNormalHermitian(self%A, N, ierr)
        call VecDuplicate(self%x_local, Ay, ierr) ! just to have same config
        call MatMultHermitianTranspose(self%A, self%y_local, Ay, ierr)
        
        ! Create the KSP linear solver
        call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
        call KSPSetOperators(ksp, N, N, ierr)
        
        ! Use LSQR solver for rectangular system (least squares problem)
        call KSPSetType(ksp, KSPCG, ierr)
        call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
        call KSPSetFromOptions(ksp, ierr)

        ! Set solver tolerances: rel_tol, abs_tol, div_tol, and max iterations
        call KSPSetTolerances(ksp, rel_tol, abs_tol, div_tol,  max_iter, ierr)

        ! Solve the linear system A * u = b
        if (rank==0) print * ,'now solving'
        call KSPSolve(ksp, Ay, self%x_local, ierr)

        ! Output the solution norm as a summary
        call VecNorm(self%x_local, NORM_2, norm, ierr)
        if (rank == 0) then
            write(*,*) 'Solution vector 2-norm:', norm
        end if

        call VecDestroy(Ay, ierr)
        call MatDestroy(N, ierr)
        call KSPDestroy(ksp, ierr)

    end subroutine

end module helmHoltzDecomp
    
    
