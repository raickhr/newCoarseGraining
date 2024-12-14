program forTestMain
        use kinds
        use operators
        use coarsening
        use interpolation
        use mpiwrapper
        use forTestReadWrite
        use helmHoltzDecomp
        use multiGridHelmHoltz

        implicit none
        integer, parameter :: nx = 2597, ny = 597                 ! Dimensions of the array

        real, allocatable, dimension(:,:) :: LNDX_RHO, RNDX_RHO, DNDY_RHO, UNDY_RHO, &
                        & UEDX_RHO, DEDX_RHO, REDY_RHO, LEDY_RHO, &
                        & DX_RHO, DY_RHO, AREA, LAT_RHO, LON_RHO, &
                        & uvel, vvel, phi, psi, &
                        & uvel_pol, uvel_tor, vvel_pol, vvel_tor, &  
                        & uvel_tot, vvel_tot      

        integer:: factorList(3), cnx, cny, shapeArr(2), ierr

        call startMPI()

        if (taskid == 0) print *,'factorList', factorList
        !factorList = (/ 8, 4, 2/)
        factorList = (/ 3, 5, 9/)
        

        if (taskid == 0) then
                allocate(LNDX_RHO(nx, ny), &
                & RNDX_RHO(nx, ny), &
                & DNDY_RHO(nx, ny), &
                & UNDY_RHO(nx, ny), &
                & UEDX_RHO(nx, ny), &
                & DEDX_RHO(nx, ny), &
                & REDY_RHO(nx, ny), &
                & LEDY_RHO(nx, ny), &
                & DX_RHO(nx, ny), &
                & DY_RHO(nx, ny), &
                & AREA(nx, ny), &
                & LAT_RHO(nx, ny), &
                & LON_RHO(nx, ny), &
                & uvel(nx, ny), &
                & vvel(nx, ny), &
                & phi(nx, ny), &
                & psi(nx, ny), &
                & uvel_pol(nx, ny), &
                & uvel_tor(nx, ny), &  
                & vvel_pol(nx, ny), &
                & vvel_tor(nx, ny), &  
                & uvel_tot(nx, ny), &
                & vvel_tot(nx, ny) )
                print *, 'READING VARIABLES ...'
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
        endif

        call doMultiGridHelmHoltz(uvel, vvel, &
                                  LAT_RHO, LON_RHO,&
                                  DX_RHO, DY_RHO, &
                                  UEDX_RHO, DEDX_RHO, &
                                  LEDY_RHO, REDY_RHO, &
                                  UNDY_RHO, DNDY_RHO, &
                                  LNDX_RHO, RNDX_RHO, &
                                  AREA, &
                                  factorList, psi, phi)

        ! call doMultiGridHelmHoltzRes(uvel, vvel, &
        !                         LAT_RHO, LON_RHO,&
        !                         DX_RHO, DY_RHO, &
        !                         UEDX_RHO, DEDX_RHO, &
        !                         LEDY_RHO, REDY_RHO, &
        !                         UNDY_RHO, DNDY_RHO, &
        !                         LNDX_RHO, RNDX_RHO, &
        !                         AREA, &
        !                         factorList, psi, phi)

        if (taskid == 0) then
                call getPolTorVel(psi, phi, DEDX_RHO, UEDX_RHO, LEDY_RHO, REDY_RHO, AREA, &
                                 uvel_pol, uvel_tor, vvel_pol, vvel_tor)

                call write2dVar('uvel_pol.nc', 'uvel_pol', uvel_pol )
                call write2dVar('uvel_tor.nc', 'uvel_tor', uvel_tor )
                
                call write2dVar('vvel_pol.nc', 'vvel_pol', vvel_pol )
                call write2dVar('vvel_tor.nc', 'vvel_tor', vvel_tor )

                call write2dVar('uvel_tot.nc', 'uvel_tot', uvel_tor + uvel_pol )
                call write2dVar('vvel_tot.nc', 'vvel_tot', vvel_tor + vvel_pol )

                call write2dVar('psi.nc', 'psi', psi )
                call write2dVar('phi.nc', 'phi', phi )
                deallocate(LNDX_RHO, RNDX_RHO, DNDY_RHO, UNDY_RHO, &
                & UEDX_RHO, DEDX_RHO, REDY_RHO, LEDY_RHO, &
                & DX_RHO, DY_RHO, AREA, LAT_RHO, LON_RHO, &
                & uvel, vvel, phi, psi, &
                & uvel_pol, uvel_tor, vvel_pol, vvel_tor, &  
                & uvel_tot, vvel_tot )
        endif

        call MPI_Barrier(MPI_COMM_WORLD, i_err)

        call stopMPI()

end program
