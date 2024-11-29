module coarsening
    use kinds
    implicit none

    contains
    subroutine coarsenField(nx, ny, factor, field, cellArea, outField, fieldWithPadding, cellAreaWithPadding)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: field(nx,  ny), cellArea(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: outField(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, intent(inout) :: fieldWithPadding(:,:) , cellAreaWithPadding(:,:)

        print *, 'Original field size', nx, ny

        print *, 'Coarsening by factor', factor

        coarse_nx = nx/factor
        coarse_ny = ny/factor

        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        endif

        print *, 'pad_x1, pad_x2, pad_y1, pad_y2', pad_x1, pad_x2, pad_y1, pad_y2
        print *, 'Coarsening to', coarse_nx, coarse_ny
        allocate(outField(coarse_nx, coarse_ny), stat=ierr)
        

        !print *,' padded field size' , nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2
        allocate(fieldWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
                &cellAreaWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2),  stat=ierr)

        fieldWithPadding(pad_x1+1:pad_x1 +nx, pad_y1+1: pad_y1 +ny) = field
        cellAreaWithPadding(pad_x1+1:pad_x1 +nx, pad_y1+1: pad_y1 +ny) = cellArea

        do i =1, pad_x1
            fieldWithPadding(i, pad_y1+1:ny+pad_y1) = field(1,:) 
            cellAreaWithPadding(i, pad_y1+1:ny+pad_y1) = cellArea(1,:)
        end do

        do i =pad_x1 + nx + 1, pad_x1 + nx + pad_x2
            fieldWithPadding(i,pad_y1+1:ny+pad_y1) = field(nx,:) 
            cellAreaWithPadding(i,pad_y1+1:ny+pad_y1) = cellArea(nx,:) 
        end do

        do i =1, pad_y1
            fieldWithPadding(:, i) = fieldWithPadding(:, pad_y1+1)
            cellAreaWithPadding(:, i) = cellAreaWithPadding(:, pad_y1+1)
        end do

        do i =pad_y1 + ny + 1, pad_y1 + ny + pad_y2
            fieldWithPadding(:, i) = fieldWithPadding(:, pad_y1+ny) 
            cellAreaWithPadding(:, i) = cellAreaWithPadding(:, pad_y1+ny) 
        end do

        is = 1
        do i = 1, coarse_nx
            ie = is + factor-1
            js = 1
            do j = 1, coarse_ny
                je = js + factor-1
                outField(i,j) = sum(fieldWithPadding(is:ie, js:je) * cellAreaWithPadding(is:ie, js:je))/ &
                                & sum(cellAreaWithPadding(is:ie, js:je))
                js = je + 1
            end do
            is = ie + 1
        end do

    end subroutine

    subroutine coarsenLatLon(nx, ny, factor, LAT, LON, cLAT, cLON, latWithPadding, lonWithPadding)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: LAT(nx,  ny), LON(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: cLAT(:,:), cLON(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, intent(inout) :: latWithPadding(:,:) , lonWithPadding(:,:)
        real(kind=real_kind), allocatable ::  weightsWithPadding(:,:), deltax_lon(:), deltay_lon(:), deltax_lat(:), deltay_lat(:)

        print *, 'Original field size', nx, ny

        print *, 'Coarsening by factor', factor

        coarse_nx = nx/factor
        coarse_ny = ny/factor

        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        endif

        !print *, 'pad_x1, pad_x2, pad_y1, pad_y2', pad_x1, pad_x2, pad_y1, pad_y2
        print *, 'Coarsening to', coarse_nx, coarse_ny
        allocate(cLAT(coarse_nx, coarse_ny),cLON(coarse_nx, coarse_ny), stat=ierr)
        

        !print *,' padded field size' , nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2
        allocate(latWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
            &    lonWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
            &    weightsWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), stat=ierr)

        latWithPadding(pad_x1+1:pad_x1 +nx, pad_y1+1: pad_y1 +ny) = LAT
        lonWithPadding(pad_x1+1:pad_x1 +nx, pad_y1+1: pad_y1 +ny) = LON
        weightsWithPadding(:,:) = 1

        allocate(deltax_lon(ny), deltax_lat(ny), stat=ierr)

        ! LEFT PADDING  ... DOES NOT FILL TOP AND BOTTOM PADDING!
        deltax_lat = LAT(1,:) - LAT(2,:)
        deltax_lon = LON(1,:) - LON(2,:)

        print *, 'deltax_lon', deltax_lon(1:10)
        
        do i =1, pad_x1
            latWithPadding(i, pad_y1+1:ny+pad_y1) = LAT(1,:) + (pad_x1 -i + 1) * deltax_lat
            lonWithPadding(i, pad_y1+1:ny+pad_y1) = LON(1,:) + (pad_x1 -i + 1) * deltax_lon
        end do

        ! RIGHT PADDING  ... DOES NOT FILL TOP AND BOTTOM PADDING!
        deltax_lat = LAT(nx,:) - LAT(nx-1,:)
        deltax_lon = LON(nx,:) - LON(nx-1,:)
        
        do i =pad_x1 + nx + 1, pad_x1 + nx + pad_x2
            latWithPadding(i,pad_y1+1:ny+pad_y1) = LAT(nx,:) + (i - pad_x1 - nx) * deltax_lat
            lonWithPadding(i,pad_y1+1:ny+pad_y1) = LON(nx,:) + (i - pad_x1 - nx) * deltax_lon
        end do

        deallocate(deltax_lon, deltax_lat)

        allocate(deltay_lon(pad_x1 + nx + pad_x2), deltay_lat(pad_x1 + nx + pad_x2), stat=ierr)

        ! DOWN PADDING  ... FILL LEFT AND RIGHT PADDING!
        deltay_lat = latWithPadding(:,pad_y1 + 1) - latWithPadding(:,pad_y1 + 2)
        deltay_lon = lonWithPadding(:,pad_y1 + 1) - lonWithPadding(:,pad_y1 + 2)

        print *, 'deltay_lat', deltay_lat(1:10)
        
        do i =1, pad_y1
            latWithPadding(:, i) = latWithPadding(:, pad_y1+1) + (pad_y1 - i + 1) * deltay_lat
            lonWithPadding(:, i) = lonWithPadding(:, pad_y1+1) + (pad_y1 - i + 1) * deltay_lon
        end do

        ! UP PADDING  ... FILL LEFT AND RIGHT PADDING!
        deltay_lat = latWithPadding(:,pad_y1 + ny) - latWithPadding(:,pad_y1 + ny - 1)
        deltay_lon = lonWithPadding(:,pad_y1 + ny) - lonWithPadding(:,pad_y1 + ny - 1)

        do i =pad_y1 + ny + 1, pad_y1 + ny + pad_y2
            latWithPadding(:, i) = latWithPadding(:, pad_y1+ny)  + (i - pad_y1 - ny) * deltay_lat
            lonWithPadding(:, i) = lonWithPadding(:, pad_y1+ny)  + (i - pad_y1 - ny) * deltay_lon
        end do


        deallocate(deltay_lon, deltay_lat)

        is = 1
        do i = 1, coarse_nx
            ie = is + factor-1
            js = 1
            do j = 1, coarse_ny
                je = js + factor-1
                cLAT(i,j) = sum(latWithPadding(is:ie, js:je) * weightsWithPadding(is:ie, js:je))/ &
                                & sum(weightsWithPadding(is:ie, js:je))
                cLON(i,j) = sum(lonWithPadding(is:ie, js:je) * weightsWithPadding(is:ie, js:je))/ &
                                & sum(weightsWithPadding(is:ie, js:je))
                js = je + 1
            end do
            is = ie + 1
        end do

    end subroutine

    subroutine coarsenDXDY(nx, ny, factor, DX, DY, cDX, cDY, dxWithPadding, dyWithPadding, downCenterUp)
        integer, intent(in):: nx, ny, factor, &
                              & downCenterUp       ! this variable is for location of the grid Distances
                                                   ! up(1) will sum the top edges, down(-1) will sum with bottom edges, center(0) will average the edges.
                                                   ! same variable is used for left center right location 
                                                   ! dx dy combination is (left, down) , (right, up), (center, center) 

        real(kind=real_kind), intent(in) :: DX(nx,  ny), DY(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: cDX(:,:), cDY(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, intent(inout) :: dyWithPadding(:,:) , dxWithPadding(:,:)
        real(kind=real_kind), allocatable ::  weightsWithPadding(:,:)

        print *, 'Original field size', nx, ny

        print *, 'Coarsening by factor', factor

        coarse_nx = nx/factor
        coarse_ny = ny/factor

        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        endif

        !print *, 'pad_x1, pad_x2, pad_y1, pad_y2', pad_x1, pad_x2, pad_y1, pad_y2
        print *, 'Coarsening to', coarse_nx, coarse_ny
        allocate(cDX(coarse_nx, coarse_ny),cDY(coarse_nx, coarse_ny), stat=ierr)
        

        !print *,' padded field size' , nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2
        allocate(dyWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
            &    dxWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
            &    weightsWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), stat=ierr)

        dxWithPadding(pad_x1+1:pad_x1 +nx, pad_y1+1: pad_y1 +ny) = DX
        dyWithPadding(pad_x1+1:pad_x1 +nx, pad_y1+1: pad_y1 +ny) = DY
        weightsWithPadding(:,:) = 1

        ! LEFT PADDING  ... DOES NOT FILL TOP AND BOTTOM PADDING!
        
        do i =1, pad_x1
            dxWithPadding(i, pad_y1+1:ny+pad_y1) = DX(1,:) 
            dyWithPadding(i, pad_y1+1:ny+pad_y1) = DY(1,:) 
        end do

        ! RIGHT PADDING  ... DOES NOT FILL TOP AND BOTTOM PADDING!
    
        do i =pad_x1 + nx + 1, pad_x1 + nx + pad_x2
            dxWithPadding(i,pad_y1+1:ny+pad_y1) = DX(nx,:) 
            dyWithPadding(i,pad_y1+1:ny+pad_y1) = DY(nx,:) 
        end do

        ! DOWN PADDING  ... FILL LEFT AND RIGHT PADDING!
        
        do i =1, pad_y1
            dxWithPadding(:, i) = dxWithPadding(:, pad_y1+1) 
            dyWithPadding(:, i) = dyWithPadding(:, pad_y1+1) 
        end do

        ! UP PADDING  ... FILL LEFT AND RIGHT PADDING!

        do i =pad_y1 + ny + 1, pad_y1 + ny + pad_y2
            dxWithPadding(:, i) = dxWithPadding(:, pad_y1+ny)  
            dyWithPadding(:, i) = dyWithPadding(:, pad_y1+ny)  
        end do

        is = 1
        do i = 1, coarse_nx
            ie = is + factor-1
            js = 1
            do j = 1, coarse_ny
                je = js + factor-1
                if (downCenterUp .EQ. -1) then 
                    cDX(i,j) = sum(dxWithPadding(is:ie, js))
                    cDY(i,j) = sum(dyWithPadding(is, js:je))

                else if (downCenterUp .EQ. 1) then 
                    cDX(i,j) = sum(dxWithPadding(is:ie, je))
                    cDY(i,j) = sum(dyWithPadding(ie, js:je))

                else if (downCenterUp .EQ. 0) then 
                    cDX(i,j) = 0.5 *( sum(dxWithPadding(is:ie, js)) + sum(dxWithPadding(is:ie, je)))
                    cDY(i,j) = 0.5 *( sum(dyWithPadding(is, js:je)) + sum(dyWithPadding(ie, js:je)))
                
                else
                    stop 'incorrect value for variable downCenterUp'
                end if
                js = je + 1
            end do
            is = ie + 1
        end do

    end subroutine

    subroutine coarsenAREA(nx, ny, factor, cellArea, coarseCellArea, cellAreaWithPadding)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: cellArea(nx,  ny)
        real(kind=real_kind), allocatable, intent(out) :: coarseCellArea(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, intent(inout) :: cellAreaWithPadding(:,:)
    
        print *, 'Original cellArea size', nx, ny
    
        print *, 'Coarsening by factor', factor
    
        coarse_nx = nx/factor
        coarse_ny = ny/factor
    
        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        endif
    
        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        endif
    
        print *, 'pad_x1, pad_x2, pad_y1, pad_y2', pad_x1, pad_x2, pad_y1, pad_y2
        print *, 'Coarsening to', coarse_nx, coarse_ny
        allocate(coarseCellArea(coarse_nx, coarse_ny), stat=ierr)
        
    
        !print *,' padded cellArea size' , nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2
        allocate(cellAreaWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2),  stat=ierr)
    
        cellAreaWithPadding(pad_x1+1:pad_x1 +nx, pad_y1+1: pad_y1 +ny) = cellArea
    
        do i =1, pad_x1
            cellAreaWithPadding(i, pad_y1+1:ny+pad_y1) = cellArea(1,:) 
        end do
    
        do i =pad_x1 + nx + 1, pad_x1 + nx + pad_x2
            cellAreaWithPadding(i,pad_y1+1:ny+pad_y1) = cellArea(nx,:) 
        end do
    
        do i =1, pad_y1
            cellAreaWithPadding(:, i) = cellAreaWithPadding(:, pad_y1+1)
        end do
    
        do i =pad_y1 + ny + 1, pad_y1 + ny + pad_y2
            cellAreaWithPadding(:, i) = cellAreaWithPadding(:, pad_y1+ny) 
        end do
    
        is = 1
        do i = 1, coarse_nx
            ie = is + factor-1
            js = 1
            do j = 1, coarse_ny
                je = js + factor-1
                coarseCellArea(i,j) = sum(cellAreaWithPadding(is:ie, js:je))
                js = je + 1
            end do
            is = ie + 1
        end do
    
    end subroutine

end module