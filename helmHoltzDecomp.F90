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

    type :: linearSystemMat
        MAT :: A
    contains 
            procedure :: setMat => set_mat
            procedure :: delMat => del_mat
    end type

    integer :: ierr, rank, size

    class(linearSystemMat), allocatable :: multiGridMats(:)

    contains

    subroutine initPETSC()
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
    end subroutine


    subroutine set_mat(self, nx, ny, dx, dy)

        class(linearSystemMat), intent(inout) :: self
        integer, intent(in) :: nx, ny
        real(kind=real_kind), intent(in) :: dx(nx, ny), dy(nx, ny)

        integer, allocatable :: relPosX(:), relPosY(:)
        real(kind=real_kind), allocatable :: fdXCoeffs(:), fdYCoeffs(:)
        integer :: schemeX, schmeY, ncoeffsX, ncoeffsY, count


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
                call getFDcoefficients(derOrder=1, accuracyOrder=6, scheme=scheme, fdXCoeffs, relPosX)
                call getFDcoefficients(derOrder=1, accuracyOrder=6, scheme=scheme, fdYCoeffs, relPosY)

                ncoeffsX = size(fdXCoeffs)
                ncoeffsY = size(fdYCoeffs)                
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
                call getFDcoefficients(derOrder=2, accuracyOrder=6, scheme=scheme, fdXCoeffs, relPosX)
                call getFDcoefficients(derOrder=2, accuracyOrder=6, scheme=scheme, fdYCoeffs, relPosY)

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

            call MatAssemblyBegin(self%A, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(self%A, MAT_FINAL_ASSEMBLY, ierr)

            !call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

            if (rank == 0 ) print *, 'Matrix assembly complete'


        enddo

    end subroutine

    subroutine set_rhs(self, uvel, vvel, nx, ny)
        class(linearSystemMat) :: self
        integer, intent(in) :: nx, ny
        real (kind=real_kind) :: uvel(nx, ny), vvel(nx, ny)
        call VecCreate(PETSC_COMM_WORLD, y_globalOnZero, ierr)
        call VecSetFromOptions(y_globalOnZero, ierr)

        call VecSetSizes(y_globalOnZero, PETSC_DECIDE, 4*mx * my, ierr)

        if (rank == 0) then
            print *, 'creating vectors of size', mx,'x', my

            do Ii = 0, mx * my - 1
                i = mod(Ii, mx)
                j = Ii / mx
                
                !!!!!!!!!!!!!!!!!!!!!!!!!! setting RHS !!!!!!!!!!!!!!!!!!!!!!!
                val = uvel(i+1, j+1) 
                call VecSetValue(y_globalOnZero, Ii, val, INSERT_VALUES, ierr)
                
                val = vvel(i+1, j+1) 
                call VecSetValue(y_globalOnZero, Ii + (mx * my), val, INSERT_VALUES, ierr)
                
                val = -divU(i+1, j+1) 
                call VecSetValue(y_globalOnZero, Ii + 2 * (mx * my), val, INSERT_VALUES, ierr)

                val = -curlU(i+1, j+1) 
                call VecSetValue(y_globalOnZero, Ii + 3 * (mx * my), val, INSERT_VALUES, ierr)
            end do

            print *, 'LHS and RHS vectors assigned values Set'
        end if

        call VecAssemblyBegin(y_globalOnZero, ierr)
        call VecAssemblyEnd(y_globalOnZero, ierr)

        ! Create the distributed destination vector
        if (rank == 0) print *, 'creating RHS vector for distributing across processors'

        call VecCreate(PETSC_COMM_WORLD, y_local, ierr)
        call VecSetSizes(y_local, PETSC_DECIDE, 4*mx * my, ierr)

        call VecSetFromOptions(y_local, ierr)
        
        if (rank == 0) print *, 'creating RHS vector for distributing across processors COMPLETE'

        call VecGetOwnershipRange(y_local, low, high, ierr) ! /* low, high are global indices */
        call ISCreateStride(PETSC_COMM_SELF, high - low, low, 1, yis_localVec, ierr)

        if (rank == 0) print *, 'scattering RHS vector  ... '
        ! Create scatter context
        call VecScatterCreate(y_globalOnZero, yis_localVec, y_local, yis_localVec, scattery, ierr)

        ! Perform the scatter operation
        call VecScatterBegin(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(scattery, y_globalOnZero, y_local, INSERT_VALUES, SCATTER_FORWARD, ierr)
        if (rank == 0) print *, 'scattering RHS vectors COMPLETE'

        ! NOW creating a coefficient distributed matrix A
    end subroutine

    subroutine set_lhs(self, LHS, nx, ny)
    end subroutine

end module helmHoltzDecomp
    
    
