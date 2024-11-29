module forTestReadWrite
    use netcdf
    use kinds
    implicit none
    contains

    subroutine read2Dvar(fileName, varName, nx, ny, varArr)
        character(len=*) :: fileName
        character(len=*) :: varName
        integer, intent(in):: nx, ny

        integer :: file_id, varid, ncerr
        real, allocatable, intent(out) :: varArr(:,:)

        allocate(varArr(nx, ny))

        ncerr = nf90_open(trim(adjustl(fileName)), nf90_nowrite, file_id)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, 'nf90_open')

        ncerr = nf90_inq_varid(file_id, varName, varid)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, 'inq variable')

        ncerr = nf90_get_var(file_id, varid, varArr)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, 'get variable')
        print *,'complete reading ', varName

        ! closing the grid file
        ncerr = nf90_close(file_id)
    end subroutine

    subroutine write2dVar(fileName, varName, varArr)
        character(len=*) , intent(in) :: fileName, varName
        integer :: arrshape(2), nx, ny
        real(kind=real_kind), intent(in) :: varArr(:, :)

        integer :: file_id, xdim_id, ydim_id, var_id, ncerr
    
        arrshape = shape(varArr)
        nx = arrshape(1)
        ny = arrshape(2)

        ncerr = nf90_create(fileName, nf90_clobber, file_id)
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'nf90_open to write')
    
        !-------------------------------------------------------------------
        !  define dimensions
        !-------------------------------------------------------------------
        print *, 'nx ny', nx, ny
        ncerr = nf90_def_dim(file_id, 'X' , nx, xdim_id)
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'dimension X')

        ncerr = nf90_def_dim(file_id, 'Y', ny, ydim_id)
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'dimension Y')
    
        !-------------------------------------------------------------------
        !-- define variable 

        ncerr = nf90_def_var(file_id, varName , nf90_double, (/xdim_id, ydim_id/), var_id)
    
        ncerr = nf90_enddef(file_id)
        if (ncerr /= nf90_noerr) stop 'at enddef'

        !-------------------------------------------------------------------
        !-- write array to the variable

        ncerr = nf90_put_var(file_id, var_id, varArr,       &
                        start = (/1, 1/), &
                        count = (/nx, ny /))

        ncerr = nf90_close(file_id)

        if (ncerr /= nf90_noerr) stop 'at close'
    end subroutine

    subroutine handle_err (status, string)

        implicit none

        integer status
        character*(*) string

        write (*,*) nf90_strerror(status),': ',string
        print *, ''
        stop 'ERROR: STOPPING'

        return

    end subroutine handle_err

end module
