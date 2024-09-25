module gridModule
    use kinds
    use constants
    use ncdf_wrapper
    use netcdf_io

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

        subroutine get_grid_nc(path, grid_filename)

            character(len=*) :: path
            character(len=*) :: grid_filename
    
            ! local variables
            character(len=600) :: file_netcdf

            integer :: ncerr, fileid

            !---- Open file
            file_netcdf = trim(path) // '/' // trim(grid_filename)
            print *,'Opening grid file : ',trim(file_netcdf)
    
            ncerr = nf_open(trim(file_netcdf), nf_nowrite, fileid)
            if ( ncerr /= nf_noerr )  call handle_err(ncerr, 'nf_open')
    
            print *,'File ',trim(file_netcdf), ' opened ..'
            print *,' '
    
            print *,'Reading file now ..'
            
            call getVar2D_real(fileid, 'DXU', DXU,  'error reading DXU')
            call getVar2D_real(fileid, 'DYU', DYU,  'error reading DYU')
            call getVar2D_real(fileid, 'ULAT', ULAT,  'error reading ULAT')
            call getVar2D_real(fileid, 'ULONG', ULONG,  'error reading ULONG')
            call getVar2D_real(fileid, 'UAREA', UAREA,  'error reading UAREA')
            call getVar2D_real(fileid, 'KMU', KMU,  'error reading KMU')
            call getVar2D_real(fileid, 'HU', HU,  'error reading HU')
            call getVar2D_real(fileid, 'FCORU', FCORU,  'error reading FCORU')

            ! closing the grid file
            ncerr = nf_close(fileid)

            print *,'Read grid file complete ..'
            print *,' '

            
        end subroutine get_grid_nc

end module


