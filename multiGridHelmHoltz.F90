module multiGridHelmHoltz
    use kinds   
    use coarsening
    use interpolation
    use mpiwrapper
    use operators
    use helmHoltzDecomp
    implicit none

    type :: grid
        real(kind=real_kind), allocatable, dimension(:,:) :: centerDx, centerDy, &
                                                             lat, lon
        integer :: nx, ny
    contains
            procedure :: setGrid => set_grid     ! method to set grid value at grid level
            procedure :: delGrid => del_grid
    end type

    integer :: nfactors
    class(grid), allocatable :: multiGrid(:)
    
    contains

        subroutine set_grid(self, coarseLevel, lat, lon)
            class(grid), intent(inout) :: self
            integer, intent(in) :: coarseLevel
            real, intent(in) :: lat(:,:), lon(:,:)
            
            integer :: nx, ny, shapeArr

            shapeArr = shape(uvel)
            nx = shapeArr(1)
            ny = shapeArr(2)

            self%nx = nx
            self%ny = ny

            call coarsenLatLon(nx, ny, coarseLevel, lat, lon, self%lat, self%lon)
            call coarsenDXDY(nx, ny, coarseLevel, centerDx, centerDy, self%centerDx, self%centerDy, downCenterUp = 0)
        end subroutine

        subroutine del_grid(self)
            class(grid), intent(inout) :: self
            deallocate(self%lat, self%lon)
            deallocate(self%centerDx, self%centerDy)
        end subroutine


        subroutine setMultiGrid(factorList, lat, lon, centerDx, centerDy)
            integer, intent(in) :: factorList(:)
            real(kind=real_kind), intent(in) :: lat(:,:), lon(:,:), centerDx(:,:), centerDy(:,:)
            integer :: shapeArr(2), i, alloc_err
            

            nfactors = size(factorList)

            allocate(multiGrid(nfactors), list_nx(nfactors), list_ny(nfactors), stat=alloc_err)
            allocate(multiGridMats(nfactors))
            
            do i = 1, nfactors
                multiGrid(i)%setGrid(factorList(i), lat, lon)
                call MPI_Barrier(MPI_COMM_WORLD, i_err)
                multiGridMats(i)%setMat(multiGrid(i)%nx, multiGrid(i)%ny, multiGrid(i)%centerDx, multiGrid(i)%centerDy)
            end do
        end subroutine
    
    

end module