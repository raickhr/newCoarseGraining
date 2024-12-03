program forTestMain
    use kinds
    use coarsening
    use forTestReadWrite
    use interpolation


    implicit none

    integer, parameter :: nx = 2597, ny = 597                 ! Dimensions of the array

    real :: LNDX_RHO(nx, ny), &
          & DNDY_RHO(nx, ny), &
          & RNDX_RHO(nx, ny), &
          & UNDY_RHO(nx, ny), &
          & DEDX_RHO(nx, ny), &
          & LEDY_RHO(nx, ny), &
          & UEDX_RHO(nx, ny), &
          & REDY_RHO(nx, ny), &
          & DX_RHO(nx, ny), &
          & DY_RHO(nx, ny), &
          & AREA(nx, ny), &
          & LAT_RHO(nx, ny), &
          & LON_RHO(nx, ny), &
          & uvel(nx, ny), &
          & vvel(nx, ny)      

    integer:: factor, cnx, cny

    real(kind=real_kind), allocatable :: Crs_LNDX_RHO(:,:), &
                                      &  Crs_DNDY_RHO(:,:), &
                                      &  Crs_RNDX_RHO(:,:), &
                                      &  Crs_UNDY_RHO(:,:), &
                                      &  Crs_DEDX_RHO(:,:), &
                                      &  Crs_LEDY_RHO(:,:), &
                                      &  Crs_UEDX_RHO(:,:), &
                                      &  Crs_REDY_RHO(:,:), &
                                      &  Crs_DX_RHO(:,:), &
                                      &  Crs_DY_RHO(:,:), &
                                      &  Crs_AREA(:,:), &
                                      &  Crs_LAT_RHO(:, :), &
                                      &  Crs_LON_RHO(:, :), &
                                      &  Crs_uvel(:, :), &  
                                      &  Crs_vvel(:, :), &
                                      &  Crs_curl(:, :), &
                                      &  Crs_div(:,:), &
                                      &  Crs_psi(:,:), &
                                      &  Crs_phi(:,:), &
                                      &  Crs_uvel_pol(:,:), &
                                      &  Crs_uvel_tor(:,:)


    real(kind=real_kind), allocatable :: padded_LNDX_RHO(:,:), &
                                      &  padded_DNDY_RHO(:,:), &
                                      &  padded_RNDX_RHO(:,:), &
                                      &  padded_UNDY_RHO(:,:), &
                                      &  padded_DEDX_RHO(:,:), &
                                      &  padded_LEDY_RHO(:,:), &
                                      &  padded_UEDX_RHO(:,:), &
                                      &  padded_REDY_RHO(:,:), &
                                      &  padded_DX_RHO(:,:), &
                                      &  padded_DY_RHO(:,:), &
                                      &  padded_AREA(:,:), &
                                      &  padded_LAT_RHO(:, :), &
                                      &  padded_LON_RHO(:, :), &
                                      &  padded_uvel(:, :)   &
                                      &  padded_vvel(:, :)    

    call read2Dvar('WPE_ROMS_grid.nc', 'LNDX_RHO', nx, ny, LNDX_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'DNDY_RHO', nx, ny, DNDY_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'RNDX_RHO', nx, ny, RNDX_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'UNDY_RHO', nx, ny, UNDY_RHO)

    call read2Dvar('WPE_ROMS_grid.nc', 'DEDX_RHO', nx, ny, DEDX_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'LEDY_RHO', nx, ny, LEDY_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'UEDX_RHO', nx, ny, UEDX_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'REDY_RHO', nx, ny, REDY_RHO)

    call read2Dvar('WPE_ROMS_grid.nc', 'DX_RHO', nx, ny, DX_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'DY_RHO', nx, ny, DY_RHO)
    
    call read2Dvar('WPE_ROMS_grid.nc', 'AREA', nx, ny, AREA)

    call read2Dvar('WPE_ROMS_grid.nc', 'LAT', nx, ny, LAT_RHO)
    call read2Dvar('WPE_ROMS_grid.nc', 'LON', nx, ny, LON_RHO)

    call read2Dvar('test.nc', 'uvel', nx, ny, uvel)
    call read2Dvar('test.nc', 'vvel', nx, ny, vvel)
    
    factor = 8

    call coarsenLatLon(nx, ny, factor, LAT_RHO, LON_RHO, crs_LAT_RHO, crs_LON_RHO, padded_LAT_RHO, padded_LON_RHO)
    call coarsenDXDY(nx, ny, factor, DX_RHO, DY_RHO, crs_DX_RHO, crs_DY_RHO, padded_DX_RHO, padded_DY_RHO, 0)
    call coarsenAREA(nx, ny, factor, AREA, Crs_AREA, padded_AREA)
    call coarsenDXDY(nx, ny, factor, LNDX_RHO, DNDY_RHO, crs_LNDX_RHO, crs_DNDY_RHO, padded_LNDX_RHO, padded_DNDY_RHO, 0)
    call coarsenDXDY(nx, ny, factor, RNDX_RHO, UNDY_RHO, crs_RNDX_RHO, crs_UNDY_RHO, padded_RNDX_RHO, padded_UNDY_RHO, 0)
    call coarsenDXDY(nx, ny, factor, DEDX_RHO, LEDY_RHO, crs_DEDX_RHO, crs_LEDY_RHO, padded_DEDX_RHO, padded_LEDY_RHO, -1)
    call coarsenDXDY(nx, ny, factor, UEDX_RHO, REDY_RHO, crs_UEDX_RHO, crs_REDY_RHO, padded_UEDX_RHO, padded_REDY_RHO, 1)

    deallocate(padded_AREA)
    call coarsenField(nx, ny, factor, uvel, AREA, Crs_uvel, padded_uvel, padded_AREA)
    deallocate(padded_AREA)
    call coarsenField(nx, ny, factor, vvel, AREA, Crs_vvel, padded_vvel, padded_AREA)

    !call blinearInterpolationLatLon(Crs_LAT_RHO(1,:), Crs_LON_RHO(:,1), LAT_RHO(1,:), LON_RHO(:,1), Crs_uvel, interpolated_uvel)

    
    call write2dVar('int_uvel.nc','int_uvel', interpolated_uvel)


    deallocate(padded_LNDX_RHO, &
            &  padded_DNDY_RHO, &
            &  padded_RNDX_RHO, &
            &  padded_UNDY_RHO, &
            &  padded_DEDX_RHO, &
            &  padded_LEDY_RHO, &
            &  padded_UEDX_RHO, &
            &  padded_REDY_RHO, &
            &  padded_DX_RHO, &
            &  padded_DY_RHO, &
            &  padded_AREA, &
            &  padded_LAT_RHO, &
            &  padded_LON_RHO )


    deallocate(Crs_LNDX_RHO, &
            &  Crs_DNDY_RHO, &
            &  Crs_RNDX_RHO, &
            &  Crs_UNDY_RHO, &
            &  Crs_DEDX_RHO, &
            &  Crs_LEDY_RHO, &
            &  Crs_UEDX_RHO, &
            &  Crs_REDY_RHO, &
            &  Crs_DX_RHO, &
            &  Crs_DY_RHO, &
            &  Crs_AREA, &
            &  Crs_LAT_RHO, &
            &  Crs_LON_RHO )




end program