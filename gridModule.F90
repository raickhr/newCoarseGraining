module gridModule
    use kinds
    use constants
    use mpiwrapper
    use configurationMod

    implicit none

    ! 2D grid info

    real (kind=real_kind), allocatable, dimension(:,:) :: &
       &   DXU, DYU,         &  ! {x,y} spacing centered at U points
       &   ULAT, ULONG,      &  ! {latitude,longitude} of U points (radians)
       &   UAREA,            &  ! area of U cells
       &   FCORU,             &   ! Coriolis parameter at U  points
       &   KMU,               &   ! mask at U points  
       &   HU                     !Bathymetry at U points      
       
    integer (kind = int_kind) :: nxu, nyu, nzu
       
       
    contains
        subroutine set_size(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            nxu = nx
            nyu = ny
            nzu = nz
        end subroutine

        subroutine init_grid()
            integer(kind=int_kind) :: ierr

            if (taskid = MASTER) call set_size(config%nx, config%ny, config%nz)
            
            call MPI_BCAST(nxu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(nyu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(nzu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            allocate(DXU(nxu, nyu),  DYU(nxu, nyu),     &
                &   ULAT(nxu, nyu), ULONG(nxu, nyu), &
                &   UAREA(nxu, nyu), FCORU(nxu, nyu), &
                &   KMU(nxu, nyu), HU(nxu, nyu), stat=ierr)

            if ( ierr /= 0 ) then
                print *,'ERROR : could not allocate xy grid at rank', taskid
                stop 999
            endif
        
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            ! reading grid file
            if (taskid = MASTER)call get_grid_nc(config%InputPath, config%gridFile)

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call MPI_BCAST(DXU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(DYU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(ULAT, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(ULONG, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(KMU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(UAREA, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(HU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(FCORU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
            
        end subroutine

        subroutine deallocate_gridVars()
            deallocate( DXU, DYU,         &  ! {x,y} spacing centered at U points
                    &   ULAT, ULONG,      &  ! {latitude,longitude} of U points (radians)
                    &   UAREA,            &  ! area of U cells
                    &   FCORU,             &   ! Coriolis parameter at U  points
                    &   KMU,               &   ! mask at U points  
                    &   HU)                     !Bathymetry at U points      )
        end subroutine

end module


