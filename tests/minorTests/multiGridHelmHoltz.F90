module multiGridHelmHoltz
    use kinds
    use coarsening
    use interpolation
    use helmHoltzDecomp
    use mpiwrapper
    implicit none

    contains
    subroutine doMultiGridHelmHoltz(uvel, vvel, &
                                  lat, lon,&
                                  centerDx, centerDy, &
                                  topEdx, bottomEdx, &
                                  leftEdy, rightEdy, &
                                  topNdy, bottomNdy, &
                                  leftNdx, rightNdx, &
                                  cellArea, &
                                  coarsenList, psi, phi)

        real (kind=real_kind), intent(in), dimension(:,:) ::uvel, vvel, &
                                                            lat, lon, &
                                                            centerDx, centerDy, &
                                                            topEdx, bottomEdx, &
                                                            leftEdy, rightEdy, &
                                                            topNdy, bottomNdy, &
                                                            leftNdx, rightNdx, &
                                                            cellArea
        
        integer, intent(in), dimension(:) :: coarsenList

        real (kind=real_kind), intent(out), dimension(:,:) ::psi, phi

        real(kind=real_kind), allocatable, dimension(:,:) ::crs_lat, crs_lon, &
                                                            crs_psi, crs_phi, &
                                                            wrk_uvel, wrk_vvel, &
                                                            wrk_lat, wrk_lon, &
                                                            wrk_centerDx, wrk_centerDy, &
                                                            wrk_topEdx, wrk_bottomEdx, &
                                                            wrk_leftEdy, wrk_rightEdy, &
                                                            wrk_topNdy, wrk_bottomNdy, &
                                                            wrk_leftNdx, wrk_rightNdx, &
                                                            wrk_cellArea, &
                                                            wrk_psi, wrk_phi

        integer :: i, factor, nfactors, nx, ny, cnx, cny, shapeArr(2), ierr

        shapeArr = shape(uvel)
        nx = shapeArr(1)
        ny = shapeArr(2)

        nfactors = size(coarsenList)

        do i = 1, nfactors+1
            factor = coarsenList(i)
            if (taskid == MASTER) then
                if (i > 1) then
                    if (allocated(crs_lat)) deallocate(crs_lat, crs_lon, crs_psi, crs_phi)
                    allocate(crs_lat(cnx, cny), &
                             crs_lon(cnx, cny), &
                             crs_phi(cnx, cny), &
                             crs_psi(cnx, cny))
                    crs_lat = wrk_lat
                    crs_lon = wrk_lon
                    crs_phi = wrk_phi
                    crs_psi = wrk_psi
                    print *, 'allocated size ', cnx, cny
                endif

                if (i <= nfactors) then
                    print *, '' 
                    print *, '' 
                    print *, '' 
                    print *, 'COARSENING THE GRID VARIABLES'
                    call coarsenLatLon(nx, ny, factor, lat, lon, wrk_lat, wrk_lon)
                    !call coarsenDXDY(nx, ny, factor, centerDx, centerDy, wrk_centerDx, wrk_centerDy, downCenterUp = 0)
                    call coarsenAREA(nx, ny, factor, cellArea, wrk_cellArea)
                    call coarsenDXDY(nx, ny, factor, leftNdx, bottomNdy, wrk_leftNdx, wrk_bottomNdy, downCenterUp = 0)
                    call coarsenDXDY(nx, ny, factor, rightNdx, topNdy, wrk_rightNdx, wrk_topNdy, downCenterUp = 0)
                    call coarsenDXDY(nx, ny, factor, bottomEdx, leftEdy, wrk_bottomEdx, wrk_leftEdy, downCenterUp = -1)
                    call coarsenDXDY(nx, ny, factor, topEdx, rightEdy, wrk_topEdx, wrk_rightEdy, downCenterUp =  1)

                    call coarsenField(nx, ny, factor, uvel, cellArea, wrk_uvel)
                    call coarsenField(nx, ny, factor, vvel, cellArea, wrk_vvel)

                    print *, 'COARSENING COMPLETE'

                    shapeArr = shape(wrk_cellArea)
                    cnx = shapeArr(1)
                    cny = shapeArr(2)
                    if (allocated(wrk_psi)) then 
                        deallocate(wrk_psi, wrk_phi)
                    endif
                    allocate(wrk_psi(cnx, cny), wrk_phi(cnx, cny), stat=ierr)
                else
                    print *, 'copyng original grid for final time decomposition'
                    if (allocated(wrk_leftNdx)) then
                        deallocate(wrk_leftNdx, wrk_bottomNdy, wrk_rightNdx, wrk_topNdy, &
                        wrk_bottomEdx, wrk_leftEdy, wrk_topEdx, wrk_rightEdy, &
                        wrk_cellArea, wrk_uvel, wrk_vvel, wrk_psi, wrk_phi)
                    endif
                    cnx = nx
                    cny = ny
                    allocate(wrk_leftNdx(cnx, cny))
                    allocate(wrk_bottomNdy(cnx, cny))
                    allocate(wrk_rightNdx(cnx, cny))
                    allocate(wrk_topNdy(cnx, cny))
                    allocate(wrk_bottomEdx(cnx, cny))
                    allocate(wrk_leftEdy(cnx, cny))
                    allocate(wrk_topEdx(cnx, cny))
                    allocate(wrk_rightEdy(cnx, cny))
                    allocate(wrk_cellArea(cnx, cny))
                    allocate(wrk_uvel(cnx, cny))
                    allocate(wrk_vvel(cnx, cny))
                    allocate(wrk_psi(cnx, cny))
                    allocate(wrk_phi(cnx, cny))

                    wrk_leftNdx= leftNdx
                    wrk_bottomNdy= bottomNdy
                    wrk_rightNdx= rightNdx
                    wrk_topNdy= topNdy
                    wrk_bottomEdx= bottomEdx
                    wrk_leftEdy= leftEdy
                    wrk_topEdx= topEdx
                    wrk_rightEdy= rightEdy
                    wrk_cellArea= cellArea
                    wrk_uvel= uvel
                    wrk_vvel= vvel
                endif

                if (i == 1) then
                    wrk_psi(:,:) = 0.0
                    wrk_phi(:,:) = 0.0
                else
                    print *, 'Interpolating '
                    !print *, 'crs_lat', crs_lat(1,:)
                    !print *, 'crs_lon', crs_lon(:,1)
                    call blinearInterpolationLatLon(crs_lat(1,:), crs_lon(:,1), wrk_lat(1,:), wrk_lon(:,1), crs_psi, wrk_psi)
                    call blinearInterpolationLatLon(crs_lat(1,:), crs_lon(:,1), wrk_lat(1,:), wrk_lon(:,1), crs_phi, wrk_phi)
                    print *, 'Interpolation complete '
                    print *, ' '
                    print *, ' '
                endif
            endif



            call MPI_BCAST(cnx, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(cny, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

            call MPI_Barrier(MPI_COMM_WORLD, i_err)
            if (taskid .NE. MASTER) then 
                if (allocated(wrk_leftEdy)) then
                    deallocate(wrk_uvel, wrk_vvel, &
                   wrk_topEdx, wrk_bottomEdx, &
                   wrk_leftEdy, wrk_rightEdy, &
                   wrk_topNdy, wrk_bottomNdy, &
                   wrk_leftNdx, wrk_rightNdx, &
                   wrk_cellArea, &
                   wrk_psi, wrk_phi)
                endif
                    allocate( wrk_leftNdx(cnx, cny), stat=ierr)
                    allocate( wrk_bottomNdy(cnx, cny), stat=ierr)
                    allocate( wrk_rightNdx(cnx, cny), stat=ierr)
                    allocate( wrk_topNdy(cnx, cny), stat=ierr)
                    allocate( wrk_bottomEdx(cnx, cny), stat=ierr)
                    allocate( wrk_leftEdy(cnx, cny), stat=ierr)
                    allocate( wrk_topEdx(cnx, cny), stat=ierr)
                    allocate( wrk_rightEdy(cnx, cny), stat=ierr)
                    allocate( wrk_cellArea(cnx, cny), stat=ierr)
                    allocate( wrk_uvel(cnx, cny),stat =ierr) 
                    allocate( wrk_vvel(cnx, cny),stat =ierr)
                    allocate( wrk_psi(cnx, cny),stat =ierr) 
                    allocate( wrk_phi(cnx, cny),stat =ierr)
            endif

            call MPI_BCAST(wrk_leftNdx, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_bottomNdy, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_rightNdx, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_topNdy, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_bottomEdx, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_leftEdy, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_topEdx, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_rightEdy, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_cellArea, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

            call MPI_BCAST(wrk_uvel, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(wrk_vvel, cnx*cny, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

            call MPI_Barrier(MPI_COMM_WORLD, i_err)


            if (taskid == 0 ) print *, 'Starting Helmholtz Decomposition'
            ! call decomposeHelmholtz_2(wrk_uvel, wrk_vvel, wrk_psi, wrk_phi, &
            !                         & wrk_bottomEdx, wrk_topEdx, wrk_leftEdy, wrk_rightEdy, &
            !                         & wrk_leftNdx, wrk_rightNdx, wrk_bottomNdy, wrk_topNdy, &
            !                         & wrk_cellArea)
            
            ! if (taskid == 0 ) print *, 'Decomposition complete'

        enddo
        if (taskid == 0) then
            psi = wrk_psi
            phi = wrk_phi
        endif

        deallocate(crs_lat, crs_lon, &
                   crs_psi, crs_phi, &
                   wrk_uvel, wrk_vvel, &
                   wrk_lat, wrk_lon, &
                   wrk_topEdx, wrk_bottomEdx, &
                   wrk_leftEdy, wrk_rightEdy, &
                   wrk_topNdy, wrk_bottomNdy, &
                   wrk_leftNdx, wrk_rightNdx, &
                   wrk_cellArea, &
                   wrk_psi, wrk_phi)
        
    end subroutine

end module