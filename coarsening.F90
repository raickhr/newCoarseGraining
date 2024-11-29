module coarsening
    use kinds
    implicit none

    contains
    subroutine coarsenField(nx, ny, factor, field, cellArea, outField)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: field(nx,  ny), cellArea(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: outField(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable :: fieldWithPadding(:,:) , cellAreaWithPadding(:,:)

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

    subroutine coarsenXYpos(nx, ny, factor, field, DX, DY)
        integer, intent(in):: nx, ny, factor
        real(kind=real_kind), intent(in) :: field(nx,  ny), cellArea(nx, ny)
        real(kind=real_kind), allocatable, intent(out) :: outField(:,:)
        
        ! local variables
        integer :: coarse_nx, coarse_ny, pad_x1, pad_x2, pad_y1, pad_y2, ierr
        integer :: i,j, dummy, is, ie, js, je
        real(kind=real_kind), allocatable :: fieldWithPadding(:,:) , cellAreaWithPadding(:,:)

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

end module