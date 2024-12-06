program forTestMain
        use kinds
        use operators
        use coarsening
        use interpolation
        use mpiwrapper
        use forTestReadWrite
        use helmHoltzDecomp

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

        integer:: factor, cnx, cny, shapeArr(2), ierr

        real(kind=real_kind), allocatable :: Crs_LNDX_RHO(:,:), &
                &  crs_DNDY_RHO(:,:), &
                &  crs_RNDX_RHO(:,:), &
                &  crs_UNDY_RHO(:,:), &
                &  crs_DEDX_RHO(:,:), &
                &  crs_LEDY_RHO(:,:), &
                &  crs_UEDX_RHO(:,:), &
                &  crs_REDY_RHO(:,:), &
                &  crs_DX_RHO(:,:), &
                &  crs_DY_RHO(:,:), &
                &  crs_AREA(:,:), &
                &  crs_LAT_RHO(:, :), &
                &  crs_LON_RHO(:, :), &
                &  crs_uvel(:, :), &  
                &  crs_vvel(:, :), &
                &  crs_curl(:, :), &
                &  crs_div(:,:), &
                &  crs_psi(:,:), &
                &  crs_phi(:,:), &
                &  crs_uvel_pol(:,:), &
                &  crs_uvel_tor(:,:), &
                &  crs_vvel_pol(:,:), &
                &  crs_vvel_tor(:,:)


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
                                        &  padded_uvel(:, :),   &
                                        &  padded_vvel(:, :) 

        call startMPI()

        if (taskid == 0) then
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

                factor = 2
                print *, 'COARSENING THE GRID VARIABLES'
                call coarsenLatLon(nx, ny, factor, LAT_RHO, LON_RHO, crs_LAT_RHO, crs_LON_RHO, padded_LAT_RHO, padded_LON_RHO)
                call coarsenDXDY(nx, ny, factor, DX_RHO, DY_RHO, crs_DX_RHO, crs_DY_RHO, padded_DX_RHO, padded_DY_RHO, 0)
                call coarsenAREA(nx, ny, factor, AREA, crs_AREA, padded_AREA)
                call coarsenDXDY(nx, ny, factor, LNDX_RHO, DNDY_RHO, crs_LNDX_RHO, crs_DNDY_RHO, padded_LNDX_RHO, padded_DNDY_RHO, 0)
                call coarsenDXDY(nx, ny, factor, RNDX_RHO, UNDY_RHO, crs_RNDX_RHO, crs_UNDY_RHO, padded_RNDX_RHO, padded_UNDY_RHO, 0)
                call coarsenDXDY(nx, ny, factor, DEDX_RHO, LEDY_RHO, crs_DEDX_RHO, crs_LEDY_RHO, padded_DEDX_RHO, padded_LEDY_RHO, -1)
                call coarsenDXDY(nx, ny, factor, UEDX_RHO, REDY_RHO, crs_UEDX_RHO, crs_REDY_RHO, padded_UEDX_RHO, padded_REDY_RHO, 1)


                print *, 'COARSENING THE FIELD VARIABLES'
                deallocate(padded_AREA)
                call coarsenField(nx, ny, factor, uvel, AREA, Crs_uvel, padded_uvel, padded_AREA)
                deallocate(padded_AREA)
                call coarsenField(nx, ny, factor, vvel, AREA, Crs_vvel, padded_vvel, padded_AREA)

                print *, 'COARSENING COMPLETE'

                shapeArr = shape(crs_AREA)
                cnx = shapeArr(1)
                cny = shapeArr(2)
                allocate(crs_psi(cnx, cny), crs_phi(cnx, cny), &
                &        crs_uvel_pol(cnx, cny), crs_uvel_tor(cnx, cny), &
                &        crs_vvel_pol(cnx, cny), crs_vvel_tor(cnx, cny), stat=ierr)

                crs_psi(:,:) = 0.0
                crs_phi(:,:) = 0.0
        endif


        call MPI_BCAST(cnx, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(cny, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

        call MPI_Barrier(MPI_COMM_WORLD, i_err)
        if (taskid .NE. MASTER) then 
                allocate( crs_LNDX_RHO(cnx, cny), stat=ierr)
                allocate( crs_DNDY_RHO(cnx, cny), stat=ierr)
                allocate( crs_RNDX_RHO(cnx, cny), stat=ierr)
                allocate( crs_UNDY_RHO(cnx, cny), stat=ierr)
                allocate( crs_DEDX_RHO(cnx, cny), stat=ierr)
                allocate( crs_LEDY_RHO(cnx, cny), stat=ierr)
                allocate( crs_UEDX_RHO(cnx, cny), stat=ierr)
                allocate( crs_REDY_RHO(cnx, cny), stat=ierr)
                allocate( crs_AREA(cnx, cny), stat=ierr)
                allocate( crs_uvel(cnx, cny),stat =ierr) 
                allocate( crs_vvel(cnx, cny),stat =ierr)
                allocate( crs_psi(cnx, cny),stat =ierr) 
                allocate( crs_phi(cnx, cny),stat =ierr)
        endif

        call MPI_BCAST(crs_LNDX_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_DNDY_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_RNDX_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_UNDY_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_DEDX_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_LEDY_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_UEDX_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_REDY_RHO, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(crs_AREA, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

        call MPI_Barrier(MPI_COMM_WORLD, i_err)

        call decomposeHelmholtz_2(crs_uvel, crs_vvel, crs_psi, crs_phi, &
                                & crs_DEDX_RHO, crs_UEDX_RHO, crs_LEDY_RHO, crs_REDY_RHO, &
                                & crs_LNDX_RHO, crs_RNDX_RHO, crs_DNDY_RHO, crs_UNDY_RHO, crs_AREA)
        
        if (taskid == 0 ) print *, 'Decomposition complete'

        if (taskid == 0) then
                call getPolTorVel(crs_psi, crs_phi, crs_DEDX_RHO, crs_UEDX_RHO, crs_LEDY_RHO, crs_REDY_RHO, crs_AREA, &
                                 crs_uvel_pol, crs_uvel_tor, crs_vvel_pol, crs_vvel_tor)

                call write2dVar('crs_uvel.nc', 'crs_uvel',crs_uvel )
                call write2dVar('crs_vvel.nc', 'crs_vvel',crs_vvel )

                call write2dVar('crs_uvel_pol.nc', 'crs_uvel_pol',crs_uvel_pol )
                call write2dVar('crs_uvel_tor.nc', 'crs_uvel_tor',crs_uvel_tor )
                
                call write2dVar('crs_vvel_pol.nc', 'crs_vvel_pol',crs_vvel_pol )
                call write2dVar('crs_vvel_tor.nc', 'crs_vvel_tor',crs_vvel_tor )

                call write2dVar('crs_uvel_tol.nc', 'crs_uvel_tot',crs_uvel_tor + crs_uvel_pol )
                call write2dVar('crs_vvel_tol.nc', 'crs_vvel_tot',crs_vvel_tor + crs_vvel_pol )

                call write2dVar('crs_psi.nc', 'crs_psi', crs_psi )
                call write2dVar('crs_phi.nc', 'crs_phi', crs_phi )
        endif

        call MPI_Barrier(MPI_COMM_WORLD, i_err)

        if (allocated(padded_LNDX_RHO)) deallocate(padded_LNDX_RHO)
        if (allocated(padded_DNDY_RHO)) deallocate(padded_DNDY_RHO)
        if (allocated(padded_RNDX_RHO)) deallocate(padded_RNDX_RHO)
        if (allocated(padded_UNDY_RHO)) deallocate(padded_UNDY_RHO)
        if (allocated(padded_DEDX_RHO)) deallocate(padded_DEDX_RHO)
        if (allocated(padded_LEDY_RHO)) deallocate(padded_LEDY_RHO)
        if (allocated(padded_UEDX_RHO)) deallocate(padded_UEDX_RHO)
        if (allocated(padded_REDY_RHO)) deallocate(padded_REDY_RHO)
        if (allocated(padded_DX_RHO)) deallocate(padded_DX_RHO)
        if (allocated(padded_DY_RHO)) deallocate(padded_DY_RHO)
        if (allocated(padded_AREA)) deallocate(padded_AREA)
        if (allocated(padded_LAT_RHO)) deallocate(padded_LAT_RHO)
        if (allocated(padded_LON_RHO)) deallocate(padded_LON_RHO )


        if (allocated(Crs_LNDX_RHO)) deallocate(Crs_LNDX_RHO)
        if (allocated(Crs_DNDY_RHO)) deallocate(Crs_DNDY_RHO)
        if (allocated(Crs_RNDX_RHO)) deallocate(Crs_RNDX_RHO)
        if (allocated(Crs_UNDY_RHO)) deallocate(Crs_UNDY_RHO)
        if (allocated(Crs_DEDX_RHO)) deallocate(Crs_DEDX_RHO)
        if (allocated(Crs_LEDY_RHO)) deallocate(Crs_LEDY_RHO)
        if (allocated(Crs_UEDX_RHO)) deallocate(Crs_UEDX_RHO)
        if (allocated(Crs_REDY_RHO)) deallocate(Crs_REDY_RHO)
        if (allocated(Crs_DX_RHO)) deallocate(Crs_DX_RHO)
        if (allocated(Crs_DY_RHO)) deallocate(Crs_DY_RHO)
        if (allocated(Crs_AREA)) deallocate(Crs_AREA)
        if (allocated(Crs_LAT_RHO)) deallocate(Crs_LAT_RHO)
        if (allocated(Crs_LON_RHO)) deallocate(Crs_LON_RHO)
        if (allocated(Crs_uvel)) deallocate(Crs_uvel)  
        if (allocated(Crs_vvel)) deallocate(Crs_vvel)
        if (allocated(Crs_curl)) deallocate(Crs_curl)
        if (allocated(Crs_div)) deallocate(Crs_div)
        if (allocated(Crs_psi)) deallocate(Crs_psi)
        if (allocated(Crs_phi)) deallocate(Crs_phi)
        if (allocated(Crs_uvel_pol)) deallocate(Crs_uvel_pol)
        if (allocated(Crs_uvel_tor)) deallocate(Crs_uvel_tor)

        call stopMPI()

end program
