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

        if (allocated(leftEdgePhi)) then
            deallocate(leftEdgePhi)    
            deallocate(rightEdgePhi)    
            deallocate(bottomEdgePhi)    
            deallocate(topEdgePhi)    
        endif

        allocate(leftEdgePhi(nx, ny), stat= ierr)
        allocate(rightEdgePhi(nx, ny), stat= ierr)
        allocate(bottomEdgePhi(nx, ny), stat= ierr)
        allocate(topEdgePhi(nx, ny), stat= ierr)

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

    ! subroutine fromFaceCenterToEdgeHermite(phi, dxLeftN, dxRightN, dyBottomN, dyTopN, latArr, lonArr, leftEdgePhi, rightEdgePhi, bottomEdgePhi, topEdgePhi)
    !     real(kind = real_kind), dimension(:,:), intent(in) :: phi, dxLeftN, dxRightN, dyBottomN, dyTopN, &
    !                                                        &  latArr, lonArr
    !     real(kind = real_kind), allocatable, dimension(:,:), intent(out) :: leftEdgePhi, rightEdgePhi, &
    !                                                                      &  bottomEdgePhi, topEdgePhi, 

    !     real(kind = real_kind), allocatable, dimension(:,:) :: x, x1, x2, fx1, fx2, dfx1, dfx2, &
    !                                                            h0, h1, h2, h3, t, dummy

    !     integer :: nx, ny, ierr, shapeArr(2)

    !     shapeArr = shape(phi)
    !     nx = shapeArr(1)
    !     ny = shapeArr(2)

        
    !     allocate(x(nx, ny), x1(nx, ny), x2(nx, ny), t(nx, ny), &
    !             fx1(nx, ny), fx2(nx, ny), dfx1(nx, ny), dfx2(nx, ny), &
    !             h0(nx, ny), h1(nx, ny), h2(nx, ny), h3(nx, ny), &
    !             dummy1(nx, ny), dummy2(nx, ny), stat= ierr)

    !     if (allocated(leftEdgePhi)) then
    !         deallocate(leftEdgePhi)    
    !         deallocate(rightEdgePhi)    
    !         deallocate(bottomEdgePhi)    
    !         deallocate(topEdgePhi)    
    !     endif

    !     allocate(leftEdgePhi(nx, ny), stat= ierr)
    !     allocate(rightEdgePhi(nx, ny), stat= ierr)
    !     allocate(bottomEdgePhi(nx, ny), stat= ierr)
    !     allocate(topEdgePhi(nx, ny), stat= ierr)


    !     ! interpolate for left Edge Value
    !     x1 = cshift(lonArr, shift=-1, dim = 1)
    !     x2 = lonArr
    !     x = (x1 + x2)/2.0
    !     x1(1,:) = lonArr(1,:)
    !     x2(1,:) = lonArr(2,:)
    !     x(1,:) = x1(1,:) - 0.5 *(x2(1,:) - x1(1,:))

    !     fx1 = cshift(phi, shift=-1, dim = 1)
    !     fx2 = phi
    !     fx1(1,:) = phi(1,:)
    !     fx2(1,:) = phi(2,:)

    !     dummy = (cshift(phi, shift=1, dim =1) - cshift(phi, shift=-1, dim =1))/(dxLeftN + dxRightN)
    !     dfx1 = cshift(dummy, shift=-1, dim =1)
    !     dfx2 = dummy

    !     t = (x - x1)/(x2 - x1)

    !     h0 = 2*t**3 - 3*t**2 + 1
    !     h1 = -2*t**3 + 3*t**2
    !     h2 = t**3 - 2*t**2 + t
    !     h3 = t**3 - t**2

    !     leftEdgePhi = h0 * fx1  &
    !                 + h1 * fx2  &
    !                 + h2 * dfx1 &
    !                 + h3 * dfx2



        

        

    !     \
    !     dummy = cshift(phi, shift=-1, dim=1) ! left values
    !     call setBoundary(dummy, boundary='left', boundaryType='linear')
    !     leftEdgePhi = 0.5 * (dummy + phi)

    !     dummy = cshift(phi, shift=1, dim=1) ! right values
    !     call setBoundary(dummy, boundary='right', boundaryType='linear')
    !     rightEdgePhi = 0.5 * (dummy + phi)

    !     dummy = cshift(phi, shift=-1, dim=2) ! bottom values
    !     call setBoundary(dummy, boundary='left', boundaryType='linear')
    !     bottomEdgePhi = 0.5 * (dummy + phi)

    !     dummy = cshift(phi, shift=1, dim=2) ! top values
    !     call setBoundary(dummy, boundary='right', boundaryType='linear')
    !     topEdgePhi = 0.5 * (dummy + phi)
    !     deallocate(dummy)
    ! end subroutine

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

    subroutine calcGradFV(phi, dxBottom, dxTop, dyLeft, dyRight, cellArea, gradX, grady)
        real(kind=real_kind), intent(in), dimension(:,:) :: phi, dxBottom, dxTop, dyLeft, dyRight, cellArea
        real(kind=real_kind), intent(out), allocatable, dimension(:,:) :: gradX, grady

        integer :: nx, ny, ierr, shapeArr(2)
        real(kind=real_kind), allocatable, dimension(:,:) :: leftEdgePhi, rightEdgePhi, bottomEdgePhi, topEdgePhi

        shapeArr = shape(phi)

        nx = shapeArr(1)
        ny = shapeArr(2)

        !Make sure the shape of the arrays match

        if (.NOT. all(shape(dxBottom) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dxTop) .EQ. shapeArr) .OR. & 
            .NOT. all(shape(dyLeft) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dyRight) .EQ. shapeArr) .OR. &
            .NOT. all(shape(cellArea) .EQ.shapeArr)) then

                stop "calcGradFV:: shape of array inconsistent"
        endif

        if (allocated(gradX)) then
            deallocate(gradX)
            deallocate(gradY)
        endif

        allocate(gradX(nx, ny), stat = ierr)
        allocate(gradY(nx, ny), stat = ierr)

        call fromFaceCenterToEdge(phi, leftEdgePhi, rightEdgePhi, bottomEdgePhi, topEdgePhi)

        gradX = (rightEdgePhi * dyRight - leftEdgePhi * dyLeft)/cellArea
        gradY = (topEdgePhi * dxTop - bottomEdgePhi * dxBottom)/cellArea

        deallocate (leftEdgePhi, rightEdgePhi, bottomEdgePhi, topEdgePhi)
    end subroutine

    subroutine calcGradFD(phi, dx, dy, gradX, grady)
        real(kind=real_kind), intent(in), dimension(:,:) :: phi, dx, dy
        real(kind=real_kind), intent(out), allocatable, dimension(:,:) :: gradX, grady

        real(kind=real_kind), allocatable, dimension(:,:) :: dummy

        integer :: nx, ny, ierr, shapeArr(2)


        shapeArr = shape(phi)

        nx = shapeArr(1)
        ny = shapeArr(2)

        !Make sure the shape of the arrays match

        if (.NOT. all(shape(dx) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dy) .EQ. shapeArr)) then
                stop "calcGradFD:: shape of array inconsistent"
        endif

        if (allocated(gradX)) then
            deallocate(gradY)
        endif

        if (allocated(gradY)) then
            deallocate(gradY)
        endif

        allocate(gradX(nx, ny), stat = ierr)
        allocate(gradY(nx, ny), stat = ierr)
        
        allocate(dummy(nx, ny), stat = ierr)
        
        gradX = 1./12. * cshift(phi, shift = -2, dim=1)
        gradX = gradX - 2./3.  * cshift(phi, shift = -1, dim=1)
        gradX = gradX + 2./3.  * cshift(phi, shift =  1, dim=1)
        gradX = gradX - 1./12. * cshift(phi, shift =  2, dim=1) 
        
        dummy = -1./2. * cshift(phi, shift = -1, dim=1)
        dummy = dummy +1./2.  * cshift(phi, shift =  1, dim=1) 

        ! gradX = 1./12. * cshift(phi, shift = -2, dim=1) &
        !      & -2./3.  * cshift(phi, shift = -1, dim=1) &
        !      & +2./3.  * cshift(phi, shift =  1, dim=1) &
        !      & -1./12. * cshift(phi, shift =  2, dim=1) 
        
        ! dummy = -1./2. * cshift(phi, shift = -1, dim=1) &
        !        & +1./2.  * cshift(phi, shift =  1, dim=1) 

        !print *, 'aaaa'
        
        gradX(2,:) = dummy(2,:)
        gradX(nx-2, :) = dummy(nx-2,:)

        dummy = phi - cshift(phi, shift = -1, dim=1) 
        gradX(nx,:) = dummy(nx,:)

        dummy = cshift(phi, shift = 1, dim=1) - phi
        gradX(1,:) = dummy(1,:)

        gradX = gradX/dx

        
        gradY = 1./12. * cshift(phi, shift = -2, dim=2) &
               -2./3.  * cshift(phi, shift = -1, dim=2) &
               +2./3.  * cshift(phi, shift =  1, dim=2) &
               -1./12. * cshift(phi, shift =  2, dim=2) 
        
        dummy = -1./2. * cshift(phi, shift = -1, dim=2) &
                +1./2.  * cshift(phi, shift =  1, dim=2) 
       
        gradY(:, 2) = dummy(:, 2)
        gradY(:, ny-2) = dummy(:, ny-2)

        dummy = phi - cshift(phi, shift = -1, dim=2) 
        gradY(:, ny) = dummy(:, ny)

        dummy = cshift(phi, shift = 1, dim=2) - phi
        gradY(:, 1) = dummy(:, 1)

        gradY = gradY/dy

        deallocate(dummy)
    end subroutine

    subroutine calcHozDivVertCurl(uvel, vvel, dxBottom, dxTop, dyLeft, dyRight, cellArea, horzDiv, vertCurl)
        real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, dxBottom, dxTop, dyLeft, dyRight, cellArea
        real(kind=real_kind), intent(out), allocatable, dimension(:,:) :: horzDiv, vertCurl

        real(kind=real_kind), allocatable, dimension(:,:) :: gradX_uvel, gradY_uvel, gradX_vvel, gradY_vvel
        integer :: nx, ny, ierr, shapeArr(2)

        shapeArr = shape(uvel)

        nx = shapeArr(1)
        ny = shapeArr(2)

        !Make sure the shape of the arrays match

        if (.NOT. all(shape(vvel) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dxBottom) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dxTop) .EQ. shapeArr) .OR. & 
            .NOT. all(shape(dyLeft) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dyRight) .EQ. shapeArr) .OR. &
            .NOT. all(shape(cellArea) .EQ.shapeArr)) then

                stop "calcHozDivVertCurl:: shape of array inconsistent"
        endif

        if (allocated(horzDiv)) then
            deallocate(horzDiv)
        endif

        if (allocated(vertCurl)) then
            deallocate(vertCurl)
        endif

        allocate(horzDiv(nx, ny), stat = ierr)
        allocate(vertCurl(nx, ny), stat = ierr)

        call calcGradFV(uvel, dxBottom, dxTop, dyLeft, dyRight, cellArea, gradX_uvel, gradY_uvel)
        call calcGradFV(vvel, dxBottom, dxTop, dyLeft, dyRight, cellArea, gradX_vvel, gradY_vvel)

        horzDiv = gradX_uvel + gradY_vvel
        vertCurl = gradX_vvel - gradY_uvel
        
        deallocate(gradX_uvel, gradY_uvel, gradX_vvel, gradY_vvel, stat=ierr)
    end subroutine


    subroutine calcHozDivVertCurlFD(uvel, vvel, dx, dy, horzDiv, vertCurl)
        real(kind=real_kind), intent(in), dimension(:,:) :: uvel, vvel, dx, dy
        real(kind=real_kind), intent(out), allocatable, dimension(:,:) :: horzDiv, vertCurl

        real(kind=real_kind), allocatable, dimension(:,:) :: gradX_uvel, gradY_uvel, gradX_vvel, gradY_vvel
        integer :: nx, ny, ierr, shapeArr(2)

        shapeArr = shape(uvel)

        nx = shapeArr(1)
        ny = shapeArr(2)

        !Make sure the shape of the arrays match

        if (.NOT. all(shape(vvel) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dx) .EQ. shapeArr) .OR. &
            .NOT. all(shape(dy) .EQ.shapeArr)) then
                stop "calcHozDivVertCurl:: shape of array inconsistent"
        endif

        if (allocated(horzDiv)) then
            deallocate(horzDiv)
        endif

        if (allocated(vertCurl)) then
            deallocate(vertCurl)
        endif

        allocate(horzDiv(nx, ny), stat = ierr)
        allocate(vertCurl(nx, ny), stat = ierr)

        

        call calcGradFD(uvel, dx, dy, gradX_uvel, gradY_uvel)
        call calcGradFD(vvel, dx, dy, gradX_vvel, gradY_vvel)
        print *, 'okay here'

        horzDiv = gradX_uvel + gradY_vvel
        vertCurl = gradX_vvel - gradY_uvel
        
        deallocate(gradX_uvel, gradY_uvel, gradX_vvel, gradY_vvel, stat=ierr)
    end subroutine

    subroutine getPolTorVel(psi, phi, centerDX, centerDy, dxBottom, dxTop, dyLeft, dyRight, cellArea, polUvel, torUvel, polVvel, torVvel)
        real(kind=real_kind), intent(in), dimension(:,:) :: psi, phi, centerDX, centerDy, dxBottom, dxTop, dyLeft, dyRight, cellArea
        real(kind=real_kind), intent(out), dimension(:,:) :: polUvel, torUvel, polVvel, torVvel

        real(kind=real_kind), allocatable, dimension(:,:) :: gradX_psi, gradY_psi, gradX_phi, gradY_phi

        !call calcGradFV(psi, dxBottom, dxTop, dyLeft, dyRight, cellArea, gradX_psi, gradY_psi)

        call calcGradFD(psi, centerDx, centerDy, gradX_psi, gradY_psi)
        print*, 'gradients psi calculation complete'
        
        !call calcGradFV(phi, dxBottom, dxTop, dyLeft, dyRight, cellArea, gradX_phi, gradY_phi)

        call calcGradFD(phi, centerDx, centerDy, gradX_phi, gradY_phi)

        print*, 'gradients calculation complete'

        polUvel = -gradX_phi
        polVvel = -gradY_phi

        torUvel = gradY_psi
        torVvel = -gradX_psi
    end subroutine
end module