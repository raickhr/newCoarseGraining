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
        dummy = dummy + 1./2.  * cshift(phi, shift =  1, dim=1) 
        
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

    subroutine calcGradFD2(phi, dx, dy, gradX, grady, accuracyOrder)
        real(kind=real_kind), intent(in), dimension(:,:) :: phi, dx, dy
        real(kind=real_kind), intent(out), allocatable, dimension(:,:) :: gradX, grady
        integer, optional :: accuracyOrder 

        integer :: nx, ny, ierr, shapeArr(2), i, derorder

        integer, allocatable ::  relPos(:)
        real(kind=real_kind), allocatable :: coefficients(:), dummy(:,:) 

        integer, parameter :: CENTRAL = 0, FORWARD = 1, BACKWARD = -1

        derorder = 1
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
        
        if (present(accuracyOrder)) then
            call getFDcoefficients(coefficients, relPos, derOrder, accuracyOrder, scheme=CENTRAL)
        else
            call getFDcoefficients(coefficients, relPos, derOrder, 4, scheme=CENTRAL)
        endif

        gradX = 0
        do i = 1, size(coefficients)
            gradX = gradX + coefficients(i) * cshift(phi, shift = relPos(i), dim=1)
        end do

        gradY = 0
        do i = 1, size(coefficients)
            gradY = gradY + coefficients(i) * cshift(phi, shift = relPos(i), dim=2)
        end do

        if (present(accuracyOrder)) then
            call getFDcoefficients(coefficients, relPos, derOrder, accuracyOrder, scheme=FORWARD)
        else
            call getFDcoefficients(coefficients, relPos, derOrder, 4, scheme=FORWARD)
        endif

        dummy = 0
        do i = 1, size(coefficients)
            dummy = dummy + coefficients(i) * cshift(phi, shift = relPos(i), dim=1)
        end do
        gradX(1:size(coefficients),:) = dummy(1:size(coefficients),:)

        dummy = 0
        do i = 1, size(coefficients)
            dummy = dummy + coefficients(i) * cshift(phi, shift = relPos(i), dim=2)
        end do
        gradY(:,1:size(coefficients)) = dummy(:,1:size(coefficients))

        if (present(accuracyOrder)) then
            call getFDcoefficients(coefficients, relPos, derOrder, accuracyOrder, scheme=BACKWARD)
        else
            call getFDcoefficients(coefficients, relPos, derOrder, 4, scheme=BACKWARD)
        endif

        dummy = 0
        do i = 1, size(coefficients)
            dummy = dummy + coefficients(i) * cshift(phi, shift = relPos(i), dim=1)
        end do
        gradX(nx-size(coefficients)+1:nx,:) = dummy(nx-size(coefficients)+1:nx,:)
        
        dummy = 0
        do i = 1, size(coefficients)
            dummy = dummy + coefficients(i) * cshift(phi, shift = relPos(i), dim=2)
        end do
        gradY(:,ny-size(coefficients)+1:ny) = dummy(:,ny-size(coefficients)+1:ny)

        deallocate(dummy)
        gradX = gradX/dx
        gradY = gradY/dy

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
                print *, shapeArr, shape(vvel), shape(dx), shape(dy)
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

    subroutine getPolTorVelFD(psi, phi, centerDX, centerDy, polUvel, torUvel, polVvel, torVvel)
        real(kind=real_kind), intent(in), dimension(:,:) :: psi, phi, centerDX, centerDy
        real(kind=real_kind), intent(out), dimension(:,:) :: polUvel, torUvel, polVvel, torVvel

        real(kind=real_kind), allocatable, dimension(:,:) :: gradX_psi, gradY_psi, gradX_phi, gradY_phi

        !call calcGradFV(psi, dxBottom, dxTop, dyLeft, dyRight, cellArea, gradX_psi, gradY_psi)

        call calcGradFD(psi, centerDx, centerDy, gradX_psi, gradY_psi)
        
        !call calcGradFV(phi, dxBottom, dxTop, dyLeft, dyRight, cellArea, gradX_phi, gradY_phi)

        call calcGradFD(phi, centerDx, centerDy, gradX_phi, gradY_phi)

        polUvel = -gradX_phi
        polVvel = -gradY_phi

        torUvel = gradY_psi
        torVvel = -gradX_psi

        deallocate(gradX_psi, gradY_psi, gradX_phi, gradY_phi)
    end subroutine


    subroutine getFDcoefficients(coefficients, relPos, derOrder, accuracyOrder, scheme)
        integer, intent(in) :: derOrder, accuracyOrder, scheme
        ! scheme are for the edge of the grid, 
        ! if 1 provides forward, if -1 backward FD scheme, if 0  use centeral FD scheme

        real(kind=real_kind), intent(out), allocatable :: coefficients(:)
        integer, intent(out), allocatable :: relPos(:)

        integer :: ncoeffs, ierr

        ncoeffs = derOrder + accuracyOrder
        if (allocated(coefficients)) then
            deallocate(coefficients, relPos)
        endif
        
        allocate(coefficients(ncoeffs), relPos(ncoeffs), stat=ierr)

        coefficients(:) = 0.0

        select case (derOrder) ! first or second derivative
        case (1) ! First derivative
            select case(scheme) ! forward central or backward scheme
            case(0) ! central scheme
                select case (accuracyOrder)
                case (2)
                    coefficients(1) = -0.5
                    coefficients(2) = 0.0
                    coefficients(3) = 0.5

                    relPos(1) = -1
                    relPos(2) = 0
                    relPos(3) = 1

                case (4)
                    coefficients(1) = 1.0/12.0
                    coefficients(2) = -2.0/3.0
                    coefficients(3) = 0.0
                    coefficients(4) = 2.0/3.0
                    coefficients(5) = -1.0/12.0

                    relPos(1) = -2
                    relPos(2) = -1
                    relPos(3) = 0
                    relPos(4) = 1
                    relPos(5) = 2

                case (6)
                    coefficients(1) = -1.0/60.0
                    coefficients(2) = 3.0/20.0
                    coefficients(3) = -3.0/4.0
                    coefficients(4) = 0.0
                    coefficients(5) = 3.0/4.0
                    coefficients(6) = -3.0/20.0
                    coefficients(7) = 1.0/60.0


                    relPos(1) = -3
                    relPos(2) = -2
                    relPos(3) = -1
                    relPos(4) = 0
                    relPos(5) = 1
                    relPos(6) = 2
                    relPos(7) = 3
                    
                end select

            case(1) ! Forward Scheme
                select case (accuracyOrder)
                case (2)
                    coefficients(1) = -3.0/2.0
                    coefficients(2) = 2.0
                    coefficients(3) = -1.0/2.0

                    relPos(1) = 0
                    relPos(2) = 1
                    relPos(3) = 2

                case (4)
                    coefficients(1) = -25/12
                    coefficients(2) = 4.0
                    coefficients(3) = -3.0
                    coefficients(4) = 4.0/3.0
                    coefficients(5) = -1.0/4.0

                    relPos(1) = 0
                    relPos(2) = 1
                    relPos(3) = 2
                    relPos(4) = 3
                    relPos(5) = 4

                case (6)
                    coefficients(1) = -49.0/20.0
                    coefficients(2) = 6.0
                    coefficients(3) = -15.0/2.0
                    coefficients(4) = 20.0/3.0
                    coefficients(5) = -15.0/4.0
                    coefficients(6) = 6.0/5.0
                    coefficients(7) = -1.0/6.0

                    relPos(1) = 0
                    relPos(2) = 1
                    relPos(3) = 2
                    relPos(4) = 3
                    relPos(5) = 4
                    relPos(6) = 5
                    relPos(7) = 6
                end select


            case(-1)! Backward Scheme
                select case (accuracyOrder)
                case (2)
                    coefficients(3) = 3.0/2.0
                    coefficients(2) = -2.0
                    coefficients(1) = 1.0/2.0

                    relPos(1) = -2
                    relPos(2) = -1
                    relPos(3) = 0

                case (4)
                    coefficients(5) = 25/12
                    coefficients(4) = -4.0
                    coefficients(3) = 3.0
                    coefficients(2) = -4.0/3.0
                    coefficients(1) = 1.0/4.0

                    relPos(1) = -4
                    relPos(2) = -3
                    relPos(3) = -2
                    relPos(4) = -1
                    relPos(5) = 0

                case (6)
                    coefficients(7) = 49.0/20.0
                    coefficients(6) = -6.0
                    coefficients(5) = 15.0/2.0
                    coefficients(4) = -20.0/3.0
                    coefficients(3) = 15.0/4.0
                    coefficients(2) = -6.0/5.0
                    coefficients(1) = 1.0/6.0

                    relPos(1) = -6
                    relPos(2) = -5
                    relPos(3) = -4
                    relPos(4) = -3
                    relPos(5) = -2
                    relPos(6) = -1
                    relPos(7) = 0
                end select
            end select

        case (2) ! Second Derivative
            select case(scheme) ! forward central or backward scheme
            case(0) ! Cenbtral Scheme
                deallocate(coefficients)
                allocate(coefficients(derOrder + accuracyOrder - 1))
                deallocate(relPos)
                allocate(relPos(derOrder + accuracyOrder - 1))
                select case (accuracyOrder)
                case (2)
                    coefficients(1) = 1.0
                    coefficients(2) = -2.0
                    coefficients(3) = 1.0

                    relPos(1) = -1
                    relPos(2) = 0
                    relPos(3) = 1

                case (4)
                    coefficients(1) = -1.0/12.0
                    coefficients(2) = 4.0/3.0
                    coefficients(3) = -5.0/2.0
                    coefficients(4) = 4.0/3.0
                    coefficients(5) = -1.0/12.0

                    relPos(1) = -2
                    relPos(2) = -1
                    relPos(3) = 0
                    relPos(4) = 1
                    relPos(5) = 2

                case (6)
                    coefficients(1) = 1.0/90.0
                    coefficients(2) = -3.0/20.0
                    coefficients(3) = 3.0/2.0
                    coefficients(4) = -49.0/18.0
                    coefficients(5) = 3.0/2.0
                    coefficients(6) = -3.0/20.0
                    coefficients(7) = 1.0/90.0

                    relPos(1) = -3
                    relPos(2) = -2
                    relPos(3) = -1
                    relPos(4) = 0
                    relPos(5) = 1
                    relPos(6) = 2
                    relPos(7) = 3
                end select
            case(1) ! Forward Scheme
                select case (accuracyOrder)
                case (2)
                    coefficients(1) = 2.0
                    coefficients(2) = -5.0
                    coefficients(3) = 4.0
                    coefficients(4) = -1.0

                    relPos(1) = 0
                    relPos(2) = 1
                    relPos(3) = 2
                    relPos(4) = 3

                case (4)
                    coefficients(1) = 15.0/4.0
                    coefficients(2) = -77.0/6.0
                    coefficients(3) = 107.0/6.0
                    coefficients(4) = -13
                    coefficients(5) = 61.0/12.0
                    coefficients(6) = -5.0/6.0

                    relPos(1) = 0
                    relPos(2) = 1
                    relPos(3) = 2
                    relPos(4) = 3
                    relPos(5) = 4
                    relPos(6) = 5

                case (6)
                    coefficients(1) = 469.0/90.0
                    coefficients(2) = -223.0/10.0
                    coefficients(3) = 879.0/20.0
                    coefficients(4) = -949.0/18.0
                    coefficients(5) = 41.0
                    coefficients(6) = -201/10.0
                    coefficients(7) = 1019.0/180.0
                    coefficients(8) = -7.0/10.0

                    relPos(1) = 0
                    relPos(2) = 1
                    relPos(3) = 2
                    relPos(4) = 3
                    relPos(5) = 4
                    relPos(6) = 5
                    relPos(7) = 6
                    relPos(8) = 7
                end select
            case(-1)! Backward Scheme
                select case (accuracyOrder)
                case (2)
                    coefficients(4) = 2.0
                    coefficients(3) = -5.0
                    coefficients(2) = 4.0
                    coefficients(1) = -1.0

                    relPos(1) = -3
                    relPos(2) = -2
                    relPos(3) = -1
                    relPos(4) = 0

                case (4)
                    coefficients(6) = 15.0/4.0
                    coefficients(5) = -77.0/6.0
                    coefficients(4) = 107.0/6.0
                    coefficients(3) = -13
                    coefficients(2) = 61.0/12.0
                    coefficients(1) = -5.0/6.0

                    relPos(1) = -5
                    relPos(2) = -4
                    relPos(3) = -3
                    relPos(4) = -2
                    relPos(5) = -1
                    relPos(6) = 0

                case (6)
                    coefficients(8) = 469.0/90.0
                    coefficients(7) = -223.0/10.0
                    coefficients(6) = 879.0/20.0
                    coefficients(5) = -949.0/18.0
                    coefficients(4) = 41.0
                    coefficients(3) = -201/10.0
                    coefficients(2) = 1019.0/180.0
                    coefficients(1) = -7.0/10.0

                    relPos(1) = -7
                    relPos(2) = -6
                    relPos(3) = -5
                    relPos(4) = -4
                    relPos(5) = -3
                    relPos(6) = -2
                    relPos(7) = -1
                    relPos(8) = 0
                end select
            end select
        
        end select
    end subroutine
end module