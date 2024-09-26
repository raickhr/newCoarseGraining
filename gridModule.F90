module gridModule
    use kinds
    use constants
    use mpiwrapper

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

        subroutine init_grid(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            integer(kind=int_kind) :: ierr

            if (.not. allocated(DXU)) then
                allocate(DXU(nx, ny),  DYU(nx, ny),     &
                &   ULAT(nx, ny), ULONG(nx, ny), &
                &   UAREA(nx, ny), FCORU(nx, ny), &
                &   KMU(nx, ny), HU(nx, ny), stat=ierr)
            endif

            if ( ierr /= 0 ) then
                print *,'ERROR : could not allocate xy grid'
                stop 999
            endif
        end subroutine

        subroutine broadCastGridInfo()
            call MPI_BCAST(nxu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(nyu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(nzu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
            if (taskid .NE. MASTER) call init_grid(nxu, nyu, nzu)
    
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


