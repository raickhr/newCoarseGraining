module netcdf_io
    
    use ncdf_wrapper
    use kinds
    implicit none

    contains

    subroutine getVar2D_real(fileid, varname_in_ncfile, varArr,  error_string )
        integer, intent(in) :: fileid
        character (*), intent(in) :: varname_in_ncfile, error_string
        real (kind=real_kind), intent(out):: varArr(:,:)

        integer :: ncerr, varid

        ncerr = nf_inq_varid(fileid, varname_in_ncfile, varid)
        if ( ncerr /= nf_noerr )  call handle_err(ncerr, error_string)
        ncerr = nf_get_var_real(fileid, varid, varArr)
        if ( ncerr /= nf_noerr )  call handle_err(ncerr, error_string)
        print *,'complete reading ', varname_in_ncfile
    end subroutine


    ! subroutine getVar2D_real_inLevel(fileid, varname_in_ncfile, varArr,  error_string )
    !     integer, intent(in) :: fileid
    !     character (len = 100), intent(in) :: varname_in_ncfile, error_string
    !     real (kind=real_kind), intent(out):: varArr(:,:)

    !     integer :: ncerr, varid

    !     ncerr = nf_inq_varid(fileid, varname_in_ncfile, varid)
    !     if ( ncerr /= nf_noerr )  call handle_err(ncerr, error_string)
    !     ncerr = nf_get_var_real(fileid, varid, varArr)
    !     if ( ncerr /= nf_noerr )  call handle_err(ncerr, error_string)
    !     print *,'complete reading ', varname_in_ncfile
    ! end subroutine


end module
