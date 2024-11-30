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
                                                             bottomLeft_lon, bottomRight_lon, topRight_lon, topLeft_lon

        logical(kind=log_kind), allocatable, dimension(:,:) :: insideBox, leftOfEdge1, leftOfEdge2, leftOfEdge3, leftOfEdge4
        real(kind=real_kind) :: x1, y1, x2, y2       
        integer :: i, j

        coarseShape = shape(coarse_field)
        cnx = coarseShape(1)
        cny = coarseShape(2)

        fineShape = shape(coarse_field)
        fnx = fineShape(1)
        fny = fineShape(2)

        allocate(bottomLeft_lat(fnx, fny), bottomRight_lat(fnx, fny), topRight_lat(fnx, fny), topLeft_lat(fnx, fny), &
               & bottomLeft_lon(fnx, fny), bottomRight_lon(fnx, fny), topRight_lon(fnx, fny), topLeft_lon(fnx, fny), stat=ierr)
        allocate(leftOfEdge1(fnx, fny), leftOfEdge2(fnx, fny), leftOfEdge3(fnx, fny), leftOfEdge4(fnx, fny), &
               & insideBox(fnx, fny), stat=ierr)

        do i = 1, cnx-1
            do j = 1, cny-1
                ! check with in box

                ! bottom edge
                x1 = coarse_LON(i, j)
                x2 = coarse_LON(i+1, j)
                y1 = coarse_LAT(i, j)
                y2 = coarse_LAT(i+1, j)
                leftOfEdge1 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0

                ! right edge
                x1 = coarse_LON(i+1, j)
                x2 = coarse_LON(i+1, j+1)
                y1 = coarse_LAT(i+1, j)
                y2 = coarse_LAT(i+1, j+1)
                leftOfEdge2 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0

                ! top edge
                x1 = coarse_LON(i+1, j+1)
                x2 = coarse_LON(i, j+1)
                y1 = coarse_LAT(i+1, j+1)
                y2 = coarse_LAT(i, j+1)
                leftOfEdge3 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0

                ! left edge
                x1 = coarse_LON(i, j+1)
                x2 = coarse_LON(i, j)
                y1 = coarse_LAT(i, j+1)
                y2 = coarse_LAT(i, j)
                leftOfEdge4 = ((x2 - x1) * (fine_LAT - y1) - (y2 - y1) *(fine_LON - x1)) .GE. 0

                insideBox = leftOfEdge1 .AND. leftOfEdge2 .AND. leftOfEdge3 .AND. leftOfEdge4

                where (insideBox)
                    bottomLeft_lat = coarse_LAT(i, j)
                end where

                where (insideBox)
                    bottomLeft_lon = coarse_LON(i, j)
                end where

                !!!!!!!!!!!
                where (insideBox)
                    bottomRight_lat = coarse_LAT(i+1, j)
                end where

                where (insideBox)
                    bottomRight_lon = coarse_LON(i+1, j)
                end where
                !!!!!!!!!!!

                where (insideBox)
                    topLeft_lat = coarse_LAT(i, j+1)
                end where

                where (insideBox)
                    topLeft_lon = coarse_LON(i, j+1)
                end where

                !!!!!!!!!!!

                where (insideBox)
                    topRight_lat = coarse_LAT(i+1, j+1)
                end where

                where (insideBox)
                    topRight_lon = coarse_LON(i+1, j+1)
                end where
                
            end do
        end do






    end subroutine
end module