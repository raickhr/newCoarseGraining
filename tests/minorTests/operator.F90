module operators
    use kinds
    implicit none

    contains

    subroutine fromFaceCenterToEdge(phi, leftEdgePhi, rightEdgePhi, bottomEdgePhi, topEdgePhi)
        real(kind = real_kind), intent(in) :: phi(:,:)
        real(kind = real_kind), allocatable, dimension(:,:), intent(out) :: leftEdgePhi, rightEdgePhi, &
                                                                         &  bottomEdgePhi, topEdgePhi

        real(kind = real_kind), allocatable, dimension(:,:) :: dummy

        integer :: nx, ny, ierr, shapeArr(2)

        shapeArr = shape(phi)
        nx = shapeArr(1)
        ny = shapeArr(2)

        if (.NOT. allocated(leftEdgePhi)) then
            allocate(leftEdgePhi(nx, ny), stat= ierr)
        endif
        if (.NOT. allocated(rightEdgePhi)) then
            allocate(rightEdgePhi(nx, ny), stat= ierr)
        endif
        if (.NOT. allocated(bottomEdgePhi)) then
            allocate(bottomEdgePhi(nx, ny), stat= ierr)
        endif
        if (.NOT. allocated(topEdgePhi)) then
            allocate(topEdgePhi(nx, ny), stat= ierr)
        endif

        allocate(dummy(nx, ny), stat= ierr)

        dummy = cshift(phi, shift=-1, dim=1) ! left values
        call setBoundary(dummy, boundary='left', boundaryType='linear')
        leftEdgePhi = 0.5 * (dummy + phi)

        dummy = cshift(phi, shift=1, dim=1) ! right values
        call setBoundary(dummy, boundary='right', boundaryType='linear')
        rightEdgePhi = 0.5 * (dummy + phi)

        dummy = cshift(phi, shift=-1, dim=2) ! bottom values
        call setBoundary(dummy, boundary='left', boundaryType='linear')
        bottomEdgePhi = 0.5 * (dummy + phi)

        dummy = cshift(phi, shift=1, dim=2) ! top values
        call setBoundary(dummy, boundary='right', boundaryType='linear')
        topEdgePhi = 0.5 * (dummy + phi)
        deallocate(dummy)
    end subroutine

    subroutine setBoundary(phi, boundary, boundaryType)
        real(kind=real_kind), intent(inout) :: phi(:,:)
        character(len=*), optional, intent(in) :: boundary, boundaryType
        character(len=50) :: bdry, bdryType

        integer :: nx, ny, shapeArr(2)
        
        shapeArr = shape(phi)
        nx = shapeArr(1)
        ny = shapeArr(2)

        if (present(boundary)) then
            bdry=boundary
        else
            bdry='all'
        endif

        if (present(boundaryType)) then
            bdryType=boundaryType
        else
            bdryType='linear'
        endif
    
        if (trim(adjustl(bdryType)) == 'linear') then
            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'left') then
                phi(1,:) = 2 * phi(2,:) - phi(1,:)
            endif

            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'right') then
                phi(nx,:) = 2 * phi(nx-1,:) - phi(nx-2,:)
            endif

            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'bottom') then
                phi(:, 1) = 2 * phi(:, 2) - phi(:, 1)
            endif

            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'top') then
                phi(:, ny) = 2 * phi(:, ny-1) - phi(:, ny-2)
            endif

        elseif (trim(adjustl(bdryType)) == 'zerograd') then
            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'left') then
                phi(1,:) = phi(2,:)
            endif

            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'right') then
                phi(nx,:) = phi(nx-1,:)
            endif

            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'bottom') then
                phi(:, 1) = phi(:, 2)
            endif

            if (trim(adjustl(bdry)) == 'all' .OR. trim(adjustl(bdry)) == 'top') then
                phi(:, ny) = phi(:, ny-1)
            endif
        endif
        
        
    end subroutine

    subroutine calcGradFV(phi, dxBottom, dxTop, dyLeft, dyRight)
        
    end subroutine
end module