module interpolation
    use kinds
    implicit none
    contains

    subroutine blinearInterpolation(coarse_LAT, coarse_LON, fine_LAT, fine_LON, coarse_field, fine_field)
        real(kind = real_kind), intent(in) :: coarse_LAT(:,:), coarse_LON(:,:), &
                                            & fine_LAT(:,:), fine_LON(:,:), &
                                            & coarse_field(:,:)
        
        real(kind=real_kind), intent(out) :: fine_field(:,:)

        integer :: cnx, cny, fnx, fny, coarseShape(2), fineShape(2), ierr

        real(kind=real_kind), allocatable, dimension(:,:) :: bottomLeft_lat, bottomRight_lat, topRight_lat, topLeft_lat, &
                                                           & bottomLeft_lon, bottomRight_lon, topRight_lon, topLeft_lon, &
                                                           & bottomLeft_field, bottomRight_field, topRight_field, topLeft_field 

        logical(kind=log_kind), allocatable, dimension(:,:) :: insideBox, nearestBox, leftOfEdge1, &
                                                               leftOfEdge2, leftOfEdge3, leftOfEdge4
        real(kind=real_kind) :: x1, y1, x2, y2       
        integer :: i, j
        real(kind=real_kind) :: counter

        coarseShape = shape(coarse_field)
        cnx = coarseShape(1)
        cny = coarseShape(2)

        fineShape = shape(coarse_field)
        fnx = fineShape(1)
        fny = fineShape(2)

        allocate(bottomLeft_lat(fnx, fny), bottomRight_lat(fnx, fny), topRight_lat(fnx, fny), topLeft_lat(fnx, fny), &
               & bottomLeft_lon(fnx, fny), bottomRight_lon(fnx, fny), topRight_lon(fnx, fny), topLeft_lon(fnx, fny), &
               & bottomLeft_field(fnx, fny), bottomRight_field(fnx, fny), topRight_field(fnx, fny), &
               & topLeft_field(fnx, fny), stat=ierr)

        allocate(leftOfEdge1(fnx, fny), leftOfEdge2(fnx, fny), leftOfEdge3(fnx, fny), leftOfEdge4(fnx, fny), &
               & insideBox(fnx, fny), nearestBox(fnx, fny), stat=ierr)

        counter = 0.0

        do i = 1, cnx-1
            do j = 1, cny-1
                counter = counter + 1.0
                ! check with in box

                ! bottom edge
                x1 = coarse_LON(i, j)
                y1 = coarse_LAT(i, j)

                x2 = coarse_LON(i+1, j)
                y2 = coarse_LAT(i+1, j)

                leftOfEdge1 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0

                ! right edge
                x1 = coarse_LON(i+1, j)
                y1 = coarse_LAT(i+1, j)

                x2 = coarse_LON(i+1, j+1)
                y2 = coarse_LAT(i+1, j+1)

                leftOfEdge2 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0

                ! top edge
                x1 = coarse_LON(i+1, j+1)
                y1 = coarse_LAT(i+1, j+1)

                x2 = coarse_LON(i, j+1)
                y2 = coarse_LAT(i, j+1)

                leftOfEdge3 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0

                ! left edge
                x1 = coarse_LON(i, j+1)
                y1 = coarse_LAT(i, j+1)

                x2 = coarse_LON(i, j)
                y2 = coarse_LAT(i, j)

                leftOfEdge4 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0


                if (i == 1) then
                    if ( j == 1 ) then
                        insideBox =  leftOfEdge2 .AND. leftOfEdge3
                    elseif (j ==  cny -1 ) then
                        insideBox = leftOfEdge1 .AND. leftOfEdge2
                    else 
                        insideBox = leftOfEdge1 .AND. leftOfEdge2 .AND. leftOfEdge3
                    endif
                elseif (i == cnx -1 ) then
                    if ( j == 1 ) then
                        insideBox =  leftOfEdge3 .AND. leftOfEdge4
                    elseif (j ==  cny - 1 ) then
                        insideBox = leftOfEdge1 .AND. leftOfEdge4
                    else 
                        insideBox = leftOfEdge1 .AND. leftOfEdge3 .AND. leftOfEdge4
                    endif
                else
                    if (j == 1) then
                        insideBox = leftOfEdge2 .AND. leftOfEdge3 .AND. leftOfEdge4
                    elseif ( j == cny - 1) then
                        insideBox = leftOfEdge1 .AND. leftOfEdge2 .AND. leftOfEdge4
                    else
                        insideBox = leftOfEdge1 .AND. leftOfEdge2 .AND. leftOfEdge3 .AND. leftOfEdge4
                    endif
                endif

                !insideBox = leftOfEdge1 .AND. leftOfEdge2 .AND. leftOfEdge3 .AND. leftOfEdge4
                
                where (insideBox .EQV. .TRUE. )
                    fine_field = i+j
                end where

                ! where (insideBox)
                !     bottomLeft_lat = coarse_LAT(i, j)
                !     bottomLeft_lon = coarse_LON(i, j)
                !     bottomRight_lat = coarse_LAT(i+1, j)
                !     bottomRight_lon = coarse_LON(i+1, j)
                !     topLeft_lat = coarse_LAT(i, j+1)
                !     topLeft_lon = coarse_LON(i, j+1)
                !     topRight_lat = coarse_LAT(i+1, j+1)
                !     topRight_lon = coarse_LON(i+1, j+1)

                !     bottomLeft_field = coarse_field(i, j)
                !     bottomRight_field = coarse_field(i+1, j)
                !     topLeft_field = coarse_field(i, j+1)
                !     topRight_field = coarse_field(i+1, j+1)

                !     fine_field = counter
                ! end where

                ! where (insideBox)
                    
                ! end where

                ! !!!!!!!!!!!
                ! where (insideBox)
                    
                ! end where

                ! where (insideBox)
                    
                ! end where
                ! !!!!!!!!!!!

                ! where (insideBox)
                    
                ! end where

                ! where (insideBox)
                    
                ! end where

                ! !!!!!!!!!!!

                ! where (insideBox)
                    
                ! end where

                ! where (insideBox)
                    
                ! end where
                
            end do
        end do
        
    end subroutine

    subroutine blinearInterpolationLatLon(coarse_LAT1d, coarse_LON1d, fine_LAT1d, fine_LON1d, coarse_field, fine_field)
        real(kind=real_kind), intent(in) :: coarse_LAT1d(:), coarse_LON1d(:), fine_LAT1d(:), fine_LON1d(:), coarse_field(:,:)
        real(kind=real_kind), intent(out) :: fine_field(:,:)

        integer, allocatable :: leftRightIndex(:,:), bottomTopIndex(:,:) 
        
        real(kind=real_kind), allocatable  ::   x(:,:), y(:,:), &
                                            &   x1(:,:), x2(:,:), y1(:,:), y2(:,:), &
                                            &   f11(:,:), f21(:,:), f22(:,:), f12(:,:)

        integer :: shapeFineField(2), nearestIndex, ierr, i, j

        shapeFineField = shape(fine_field)

        allocate(leftRightIndex(shapeFineField(1), 2), bottomTopIndex(shapeFineField(2),2), stat=ierr)

        allocate(x(shapeFineField(1),shapeFineField(2)), y(shapeFineField(1),shapeFineField(2)), stat = ierr)

        allocate(x1(shapeFineField(1),shapeFineField(2)), x2(shapeFineField(1),shapeFineField(2)), stat = ierr)

        allocate(y1(shapeFineField(1),shapeFineField(2)), y2(shapeFineField(1),shapeFineField(2)), stat = ierr)

        allocate(f11(shapeFineField(1),shapeFineField(2)), f21(shapeFineField(1),shapeFineField(2)), &
        &        f22(shapeFineField(1),shapeFineField(2)), f12(shapeFineField(1),shapeFineField(2)), stat = ierr)

        
        ! co-ordiante values at the corners
        do i =1, shapeFineField(1)
            x(i,:) = fine_LON1d(i)
            leftRightIndex(i,:) = enclosingTwoIndices(coarse_LON1d, fine_LON1d(i))
            x1(i, :) = coarse_LON1d(leftRightIndex(i,1))
            x2(i, :) = coarse_LON1d(leftRightIndex(i,2))
            
            if ((x2(i,1) - x1(i,1)) < 0.0000001 ) then
                print *, 'check lon vals'
                stop
            endif
        enddo

        do i =1, shapeFineField(2)
            y(:,i) = fine_LAT1d(i)
            bottomTopIndex(i,:) = enclosingTwoIndices(coarse_LAT1d, fine_LAT1d(i))
            y1(:, i) = coarse_LAT1d(bottomTopIndex(i,1))
            y2(:, i) = coarse_LAT1d(bottomTopIndex(i,2))
            if ((y2(1,i) - y1(1,i)) < 0.0000001 ) then
                print *, coarse_LAT1d
                print *, 'size(coarse_LAT1d)', size(coarse_LAT1d)
                do j = bottomTopIndex(i,1) -2,  bottomTopIndex(i,2) + 2
                    print *, j, coarse_LAT1d(j)
                end do
                print *, 'i, bottomTopIndex(i,1), bottomTopIndex(i,2)', i, bottomTopIndex(i,1), bottomTopIndex(i,2)
                print *, 'y1(1, i), y2(1, i)', y1(1, i), y2(1, i)
                print *, 'check lat vals'
                stop
            endif
        enddo

        ! field values at the corners
        do i = 1, shapeFineField(1)
            do j = 1, shapeFineField(2)
                f11(i,j) = coarse_field(leftRightIndex(i,1), bottomTopIndex(j,1))
                f21(i,j) = coarse_field(leftRightIndex(i,2), bottomTopIndex(j,1))
                f22(i,j) = coarse_field(leftRightIndex(i,2), bottomTopIndex(j,2))
                f12(i,j) = coarse_field(leftRightIndex(i,1), bottomTopIndex(j,2))
                
            end do
        end do

        fine_field = 1/((x2 - x1)*(y2 - y1)) * &
                    & (f11 * (x2 - x  ) * (y2 - y ) + &
                    &  f21 * (x  - x1 ) * (y2 - y ) + &
                    &  f22 * (x  - x1 ) * (y  - y1) + &
                    &  f12 * (x2-  x  ) * (y  - y1))


        deallocate(leftRightIndex, bottomTopIndex)
        deallocate(x, y, x1, x2, y1, y2)
        deallocate(f11, f21, f22, f12)


    end subroutine

    function enclosingTwoIndices(coordArray, coordVal)
        real(kind=real_kind), intent(in) :: coordArray(:), coordVal
        ! coordinate should be in ascending order
        integer :: enclosingTwoIndices(2), nearestIndex

        integer :: lenCoord, ierr, i
        real(kind=real_kind), allocatable :: diffArr(:)
        real(kind=real_kind) :: minVal

        lenCoord = size(coordArray)

        allocate(diffArr(lenCoord), stat=ierr)

        diffArr = abs(coordArray - coordVal)

        minVal = diffArr(1)
        nearestIndex = 1

        do i =1, lenCoord
            if (diffArr(i) < minVal) then
                nearestIndex = i
                minVal = diffArr(i)
            endif
        end do

        if (nearestIndex == 1 .AND. coordArray(nearestIndex) > coordVal ) then
            ! have to extraplolate 
            enclosingTwoIndices(1) = 1
            enclosingTwoIndices(2) = 2
        elseif (nearestIndex == lenCoord .AND. coordArray(nearestIndex) < coordVal ) then
            ! have to extrapolate
            enclosingTwoIndices(1) = lenCoord-1
            enclosingTwoIndices(2) = lenCoord
        else
            if (coordArray(nearestIndex) > coordVal) then
                enclosingTwoIndices(1) = nearestIndex -1
                enclosingTwoIndices(2) = nearestIndex
            else
                enclosingTwoIndices(1) = nearestIndex
                enclosingTwoIndices(2) = nearestIndex+1
            endif
        endif

        if (enclosingTwoIndices(1) < 1) then
            enclosingTwoIndices(1) = 1
            enclosingTwoIndices(2) = 2
        endif

        if (enclosingTwoIndices(2) > lenCoord) then
            enclosingTwoIndices(1) = lenCoord-1
            enclosingTwoIndices(2) = lenCoord
        endif




    end function

    subroutine biliearInterpolationLatLonResidual(coarse_LAT1d, coarse_LON1d, &
                        fine_LAT1d, fine_LON1d, coarse_residual, fine_residual)

        real(kind=real_kind), intent(in) :: coarse_LAT1d(:), coarse_LON1d(:), &
                                            fine_LAT1d(:), fine_LON1d(:), coarse_residual(:)

        real(kind=real_kind), allocatable, intent(out) :: fine_residual(:)

        real(kind=real_kind), allocatable :: crs_dummy(:,:), fine_dummy(:, :)


        integer :: nx, ny, cnx, cny

        print *, 'interpolating residual'

        cny = size(coarse_LAT1d)
        cnx = size(coarse_LON1d)

        ny = size(fine_LAT1d)
        nx = size(fine_LON1d)

        print *, 'from ', cnx, 'x', cny, ' to ', nx, 'x', ny

        allocate(fine_dummy(nx, ny))
        allocate(crs_dummy(cnx, cny))
        if (allocated(fine_residual)) deallocate(fine_residual)
        allocate(fine_residual(4*nx*ny))

        print *, 'shape of coarse residual ', size(coarse_residual)

        crs_dummy = reshape(coarse_residual(1 : cnx* cny), (/cnx, cny/))
        print *, 'reshaped uvel'
        call blinearInterpolationLatLon(coarse_LAT1d, coarse_LON1d, &
                                        fine_LAT1d, fine_LON1d, &
                                        crs_dummy, fine_dummy)
        print *, 'interpolation uvel complete'
        print *, 'shape uvel', shape(fine_dummy)
        fine_residual(1:nx*ny) = pack(fine_dummy, .TRUE.)
        print *, 'flattening uvel complete'

        crs_dummy = reshape(coarse_residual(cnx* cny+1:2*cnx* cny), (/cnx, cny/))
        call blinearInterpolationLatLon(coarse_LAT1d, coarse_LON1d, &
                                        fine_LAT1d, fine_LON1d, &
                                        crs_dummy, fine_dummy)
        fine_residual(nx*ny+1: 2*nx*ny) = pack(fine_dummy, .TRUE.)

        crs_dummy = reshape(coarse_residual(2*cnx* cny+1:3*cnx* cny), (/cnx, cny/))
        call blinearInterpolationLatLon(coarse_LAT1d, coarse_LON1d, &
                                        fine_LAT1d, fine_LON1d, &
                                        crs_dummy, fine_dummy)
        fine_residual(2*nx*ny+1: 3*nx*ny) = pack(fine_dummy, .TRUE.)

        crs_dummy = reshape(coarse_residual(3*cnx* cny+1:4*cnx* cny), (/cnx, cny/))
        call blinearInterpolationLatLon(coarse_LAT1d, coarse_LON1d, &
                                        fine_LAT1d, fine_LON1d, &
                                        crs_dummy, fine_dummy)
        fine_residual(3*nx*ny+1: 4*nx*ny) = pack(fine_dummy, .TRUE.)

        deallocate(fine_dummy, crs_dummy)

    end subroutine

    subroutine biliearInterpolationLatLonLHS(coarse_LAT1d, coarse_LON1d, &
        fine_LAT1d, fine_LON1d, coarse_LHS, fine_LHS)

        real(kind=real_kind), intent(in) :: coarse_LAT1d(:), coarse_LON1d(:), &
                                    fine_LAT1d(:), fine_LON1d(:), coarse_LHS(:)

        real(kind=real_kind), allocatable, intent(out) :: fine_LHS(:)

        real(kind=real_kind), allocatable :: crs_dummy(:,:), fine_dummy(:, :)


        integer :: nx, ny, cnx, cny

        print *, 'interpolating residual'

        cny = size(coarse_LAT1d)
        cnx = size(coarse_LON1d)

        ny = size(fine_LAT1d)
        nx = size(fine_LON1d)

        print *, 'from ', cnx, 'x', cny, ' to ', nx, 'x', ny

        allocate(fine_dummy(nx, ny))
        allocate(crs_dummy(cnx, cny))
        if (allocated(fine_LHS)) deallocate(fine_LHS)
        allocate(fine_LHS(2*nx*ny))

        print *, 'shape of coarse LHS ', size(coarse_LHS)

        crs_dummy = reshape(coarse_residual(1 : cnx* cny), (/cnx, cny/))
        print *, 'reshaped phi'
        call blinearInterpolationLatLon(coarse_LAT1d, coarse_LON1d, &
                                fine_LAT1d, fine_LON1d, &
                                crs_dummy, fine_dummy)
        print *, 'interpolation phi complete'
        print *, 'shape phi', shape(fine_dummy)
        fine_residual(1:nx*ny) = pack(fine_dummy, .TRUE.)
        print *, 'flattening phi complete'

        crs_dummy = reshape(coarse_residual(cnx* cny+1:2*cnx* cny), (/cnx, cny/))
        call blinearInterpolationLatLon(coarse_LAT1d, coarse_LON1d, &
                                fine_LAT1d, fine_LON1d, &
                                crs_dummy, fine_dummy)
        fine_LHS(nx*ny+1: 2*nx*ny) = pack(fine_dummy, .TRUE.)

        deallocate(fine_dummy, crs_dummy)

        end subroutine


end module