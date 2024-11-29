program forTestMain
    use kinds
    use coarsening
    use forTestReadWrite


    implicit none

    integer, parameter :: nx = 100, ny = 100                 ! Dimensions of the array
    real :: random_array(nx, ny), cellArea(nx, ny)              ! Declare the array

    integer:: factor, cnx, cny 
    real(kind=real_kind), allocatable :: outField(:,:), fieldWithPadding(:,:), cellAreaWithPadding(:,:)

    ! Fill the array with random numbers between 0 and 1
    call random_number(random_array)
    cellArea(:, :) = 1.0

    factor = 3

    call coarsenField(nx, ny, factor, random_array, cellArea, outField, fieldWithPadding, cellAreaWithPadding)

    call write2dVar('coaseningTest.nc','test', outField)

    deallocate(cellAreaWithPadding)
    deallocate(fieldWithPadding)
    deallocate(outField)



end program