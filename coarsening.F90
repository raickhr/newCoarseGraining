module coarsening
    use kinds
    implicit none
    logical:: verbose

#ifdef VERBOSE
    parameter (verbose = .TRUE.)
#else
    parameter (verbose = .FALSE.)
#endif

    contains
    subroutine coarsenField(nx, ny, factor, field, cellArea, outField, fieldWithPadding_out, cellAreaWithPadding_out)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: field(nx,  ny), cellArea(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: outField(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, intent(inout), optional :: fieldWithPadding_out(:,:) , &
                                                                      cellAreaWithPadding_out(:,:)

        real(kind=real_kind), allocatable :: fieldWithPadding(:,:) , &
                                             cellAreaWithPadding(:,:)

        if (verbose) then
            print *, 'Coarsening field '
            print *, 'Original field size', nx, ny
            print *, 'Coarsening by factor', factor
        endif 

        coarse_nx = nx/factor
        coarse_ny = ny/factor

        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        else
            pad_x1 = 0
            pad_x2 = 0
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        else
            pad_y1 = 0
            pad_y2 = 0
        endif

        if (verbose) then
            print *, 'Padding for x and y direction'
            print *, 'left pad x1, right pad x2, bottom pad y1, top pad y2', pad_x1, pad_x2, pad_y1, pad_y2
            print *, 'Grid Size with padding', coarse_nx, coarse_ny
        end if

        allocate(outField(coarse_nx, coarse_ny), stat=ierr)
        
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

        if (present(fieldWithPadding_out) .and. present(cellAreaWithPadding_out)) then
            if (allocated(fieldWithPadding_out)) deallocate(fieldWithPadding_out)
            if (allocated(cellAreaWithPadding_out)) deallocate(cellAreaWithPadding_out)
            allocate(fieldWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
                    &cellAreaWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2),  stat=ierr)
            fieldWithPadding_out = fieldWithPadding
            cellAreaWithPadding_out = cellAreaWithPadding
        endif

        deallocate(fieldWithPadding, cellAreaWithPadding)

        if (verbose) then
            print *, 'Finished coarsening fields'
        endif

    end subroutine

    subroutine coarsenFieldWeighted(nx, ny, factor, field, cellArea, outField, fieldWithPadding_out, cellAreaWithPadding_out)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: field(nx,  ny), cellArea(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: outField(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, intent(inout), optional :: fieldWithPadding_out(:,:) , &
                                                                      cellAreaWithPadding_out(:,:)

        real(kind=real_kind), allocatable :: fieldWithPadding(:,:) , &
                                             cellAreaWithPadding(:,:), &
                                             kernelShapeX(:,:), kernelShapeY(:,:), &
                                             kernelShapeR(:,:)

        print *, 'Original field size', nx, ny
        print *, 'Coarsening by factor', factor

        coarse_nx = nx/factor
        coarse_ny = ny/factor

        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        else
            pad_x1 = 0
            pad_x2 = 0
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        else
            pad_y1 = 0
            pad_y2 = 0
        endif

        print *, 'pad_x1, pad_x2, pad_y1, pad_y2', pad_x1, pad_x2, pad_y1, pad_y2
        print *, 'Coarsening to', coarse_nx, coarse_ny
        if (allocated(outField)) deallocate(outField)
        allocate(outField(coarse_nx, coarse_ny), stat=ierr)
        

        !print *,' padded field size' , nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2
        allocate(fieldWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
                &cellAreaWithPadding(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), stat=ierr)

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

        
        ! Gaussian shape for kernel
        allocate(kernelShapeX(factor, factor), kernelShapeY(factor, factor), &
                 kernelShapeR(factor, factor), stat= ierr )

        kernelShapeX(:,1) = (/(i, i=0,factor-1)/)
        kernelShapeX(:,1) = kernelShapeX(:,1)/factor
        kernelShapeX(:,1) = kernelShapeX(:,1) - sum(kernelShapeX(:,1))/factor

        do i=1, factor
            kernelShapeX(:,i) = kernelShapeX(:,1)
            kernelShapeY(i,:) = kernelShapeX(:,1)
        end do

        kernelShapeR = sqrt(kernelShapeX**2 + kernelShapeY**2)
        kernelShapeR = exp(-(kernelShapeR**2/0.04)/2)
        kernelShapeR = kernelShapeR/sum(kernelShapeR)

        is = 1
        do i = 1, coarse_nx
            ie = is + factor-1
            js = 1
            do j = 1, coarse_ny
                je = js + factor-1
                outField(i,j) = sum(fieldWithPadding(is:ie, js:je) * cellAreaWithPadding(is:ie, js:je) * kernelShapeR)/ &
                                & sum(cellAreaWithPadding(is:ie, js:je) * kernelShapeR)
                js = je + 1
            end do
            is = ie + 1
        end do

        if (present(fieldWithPadding_out) .and. present(cellAreaWithPadding_out)) then
            if (allocated(fieldWithPadding_out)) deallocate(fieldWithPadding_out)
            if (allocated(cellAreaWithPadding_out)) deallocate(cellAreaWithPadding_out)
            allocate(fieldWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
                    &cellAreaWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2),  stat=ierr)
            fieldWithPadding_out = fieldWithPadding
            cellAreaWithPadding_out = cellAreaWithPadding
        endif

        deallocate(fieldWithPadding, cellAreaWithPadding)

    end subroutine

    subroutine coarsenLatLon(nx, ny, factor, LAT, LON, cLAT, cLON, latWithPadding_out, lonWithPadding_out)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: LAT(nx,  ny), LON(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: cLAT(:,:), cLON(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, optional, intent(inout) :: latWithPadding_out(:,:) , lonWithPadding_out(:,:)
        real(kind=real_kind), allocatable ::  weightsWithPadding(:,:), deltax_lon(:), &
                                              deltay_lon(:), deltax_lat(:), deltay_lat(:), &
                                              latWithPadding(:,:) , lonWithPadding(:,:)

        if (verbose) then
            print *, 'Coarsening lat lon by factor', factor
            print *, 'Original field size', nx, ny
        endif 

        coarse_nx = nx/factor
        coarse_ny = ny/factor

        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        else
            pad_x1 = 0
            pad_x2 = 0
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        else
            pad_y1 = 0
            pad_y2 = 0
        endif

        if (verbose) then
            print *, 'Padding for x and y direction'
            print *, 'left pad x1, right pad x2, bottom pad y1, top pad y2', pad_x1, pad_x2, pad_y1, pad_y2
            print *, 'Grid Size with padding', coarse_nx, coarse_ny
        end if

        allocate(cLAT(coarse_nx, coarse_ny),cLON(coarse_nx, coarse_ny), stat=ierr)
        

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

        if (present(latWithPadding_out).and. present(lonWithPadding_out)) then
            if (allocated(latWithPadding_out)) deallocate(latWithPadding_out)
            if (allocated(lonWithPadding_out)) deallocate(lonWithPadding_out)
            allocate(latWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
            &    lonWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2))
            latWithPadding_out = latWithPadding
            lonWithPadding_out = lonWithPadding
        endif

        deallocate(latWithPadding, lonWithPadding, weightsWithPadding)

        if (verbose) then
            print *, 'coarsening latitude and longitude complete'
        endif 

    end subroutine

    subroutine coarsenDXDY(nx, ny, factor, DX, DY, cDX, cDY, dxWithPadding_out, dyWithPadding_out, downCenterUp)
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
        real(kind=real_kind), allocatable, optional, intent(inout) :: dyWithPadding_out(:,:) , dxWithPadding_out(:,:)
        real(kind=real_kind), allocatable :: dyWithPadding(:,:) , dxWithPadding(:,:)
        real(kind=real_kind), allocatable ::  weightsWithPadding(:,:)


        if (verbose) then
            print *, 'Coarsening dx, dy by factor', factor
            print *, 'Original field size', nx, ny
            print *, 'DownCenterUP value', downCenterUp
            print *, 'Down   (-1):   Bottom Edge and Left Edge of Grid cell'
            print *, 'Center ( 0):   Cell Width and Cell Height of Grid cell'
            print *, 'Up     ( 1):   Top Edge and Right Edge of Grid cell'
        endif

        coarse_nx = nx/factor
        coarse_ny = ny/factor

        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        else
            pad_x1 = 0
            pad_x2 = 0
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        else
            pad_y1 = 0
            pad_y2 = 0
        endif

        if (verbose) then
            print *, 'Padding for x and y direction'
            print *, 'left pad x1, right pad x2, bottom pad y1, top pad y2', pad_x1, pad_x2, pad_y1, pad_y2
            print *, 'Grid Size with padding', coarse_nx, coarse_ny
        end if

        allocate(cDX(coarse_nx, coarse_ny), cDY(coarse_nx, coarse_ny), stat=ierr)
        

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

        if (present(dyWithPadding_out) .and. present(dxWithPadding_out)) then
            if (allocated(dyWithPadding_out)) deallocate(dyWithPadding_out)
            if (allocated(dxWithPadding_out)) deallocate(dxWithPadding_out)

            allocate(dyWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), &
            &    dxWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2), stat=ierr)
            
            dxWithPadding_out = dxWithPadding
            dyWithPadding_out = dyWithPadding
        endif
        deallocate(dxWithPadding, dyWithPadding, weightsWithPadding)

    end subroutine

    subroutine coarsenAREA(nx, ny, factor, cellArea, coarseCellArea, cellAreaWithPadding_out)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: cellArea(nx,  ny)
        real(kind=real_kind), allocatable, intent(out) :: coarseCellArea(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable, optional, intent(inout) :: cellAreaWithPadding_out(:,:)
        real(kind=real_kind), allocatable :: cellAreaWithPadding(:,:)
    
        print *, 'Original cellArea size', nx, ny
    
        print *, 'Coarsening by factor', factor
    
        coarse_nx = nx/factor
        coarse_ny = ny/factor
    
        if (mod(nx,factor) > 0) then
            coarse_nx = coarse_nx + 1
            dummy = coarse_nx * factor - nx
            pad_x1 = dummy/2 
            pad_x2 = dummy - pad_x1
        else
            pad_x1 = 0
            pad_x2 = 0
        endif

        if (mod(ny,factor) > 0) then
            coarse_ny = coarse_ny + 1
            dummy = coarse_ny * factor - ny
            pad_y1 = dummy/2
            pad_y2 = dummy - pad_y1
        else
            pad_y1 = 0
            pad_y2 = 0
        endif
    
        print *, 'pad_x1, pad_x2, pad_y1, pad_y2', pad_x1, pad_x2, pad_y1, pad_y2
        print *, 'Coarsening to', coarse_nx, coarse_ny
        if (allocated(coarseCellArea)) deallocate(coarseCellArea)
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

        if (present(cellAreaWithPadding_out)) then
            if (allocated(cellAreaWithPadding_out)) deallocate(cellAreaWithPadding_out)
            allocate(cellAreaWithPadding_out(nx + pad_x1 + pad_x2, ny+ pad_y1 + pad_y2),  stat=ierr)
            cellAreaWithPadding_out = cellAreaWithPadding
        endif
        deallocate(cellAreaWithPadding)
    
    end subroutine


    subroutine coarsen_residual(nx, ny, factor, residual, cellArea, crs_residual)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: residual(:), cellArea(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: crs_residual(:)

        real(kind=real_kind), allocatable :: crs_dummy(:,:), fine_dummy(:,:)
        integer :: crsArrShape(2), cnx, cny

        allocate(fine_dummy(nx,ny))

        fine_dummy = reshape(residual(1:nx*ny), shape=(/nx, ny/))
        call coarsenFieldWeighted(nx, ny, factor, fine_dummy, cellArea, crs_dummy)
        crsArrShape = shape(crs_dummy)
        cnx = crsArrShape(1)
        cny = crsArrShape(2)

        if (allocated(crs_residual)) deallocate(crs_residual)
        allocate(crs_residual(4*cnx*cny))
        crs_residual(1:cnx*cny) = pack(crs_dummy, .true.)

        fine_dummy = reshape(residual(nx*ny+1:2*nx*ny), shape=(/nx, ny/))
        call coarsenFieldWeighted(nx, ny, factor, fine_dummy, cellArea, crs_dummy)
        crs_residual(cnx*cny+1 : 2*cnx*cny) = pack(crs_dummy, .true.)

        fine_dummy = reshape(residual(2*nx*ny+1:3*nx*ny), shape=(/nx, ny/))
        call coarsenFieldWeighted(nx, ny, factor, fine_dummy, cellArea, crs_dummy)
        crs_residual(2*cnx*cny+1 : 3*cnx*cny) = pack(crs_dummy, .true.)


        fine_dummy = reshape(residual(3*nx*ny+1:4*nx*ny), shape=(/nx, ny/))
        call coarsenFieldWeighted(nx, ny, factor, fine_dummy, cellArea, crs_dummy)
        crs_residual(3*cnx*cny+1 : 4*cnx*cny) = pack(crs_dummy, .true.)

        deallocate(crs_dummy, fine_dummy)

    end subroutine

    ! subroutine coarsenGridAndUVforHelmHoltz(uvel, vvel, &
    !                                         lat, lon, &
    !                                         centerDx, centerDy, &
    !                                         crs_uvel, crs_vvel, &
    !                                         crs_lat, crs_lon, &
    !                                         crs_centerDx, crs_centerDy)

    !     real (kind=real_kind), intent(in), dimension(:,:) ::uvel, vvel, &
    !                                                         lat, lon, &
    !                                                         centerDx, centerDy

    !     real(kind=real_kind), intent(out), allocatable, dimension(:,:) ::crs_uvel, crs_vvel, &
    !                                                         crs_lat, crs_lon, &
    !                                                         crs_centerDx, crs_centerDy


    !     integer :: shapeArr(2), nx, ny
    
    !     shapeArr = shape(uvel)
    !     nx = shapeArr(1)
    !     ny = shapeArr(2)
        
    !     call coarsenLatLon(nx, ny, factor, lat, lon, wrk_lat, wrk_lon)
    !     call coarsenDXDY(nx, ny, factor, centerDx, centerDy, wrk_centerDx, wrk_centerDy, downCenterUp = 0)
    !     call coarsenField(nx, ny, factor, uvel, cellArea, wrk_uvel)
    !     call coarsenField(nx, ny, factor, vvel, cellArea, wrk_vvel)
    ! end subroutine

    
end module