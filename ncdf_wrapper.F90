module ncdf_wrapper

    !include 'netcdf.inc'
    use netcdf
    contains

    subroutine handle_err (status, string)

        implicit none

        integer status
        character*(*) string

        write (*,*) nf90_strerror(status),': ',string
        print *, ''
        stop 'ERROR: STOPPING'

        return

    end subroutine handle_err

end module ncdf_wrapper
