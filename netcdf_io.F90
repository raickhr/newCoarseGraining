module netcdf_io
    use gridModule
    use ncdf_wrapper
    use kinds
    use fields
    implicit none

    contains

    subroutine getTimeVar(fileid, varname_in_ncfile, units, calendar, time_index, value, error_string)
        integer, intent(in) :: fileid, time_index
        character(*), intent(in) :: varname_in_ncfile, error_string
        character(len=units_len), intent(out) :: units, calendar
        real(kind=real_kind), intent(out) :: value(1)

        integer :: ncerr, varid, startcount(1), endcount(1)

        startcount = (/time_index/)
        endcount = (/1/)
        ncerr = nf90_inq_varid(fileid, varname_in_ncfile, varid)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)

        ncerr = nf90_get_att(fileid, varid, "units", units)
        if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining units of "//varname_in_ncfile)

        ncerr = nf90_get_att(fileid, varid, "calendar", calendar)
        if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining calendar of "//varname_in_ncfile)

        ncerr = nf90_get_var(fileid, varid, value, startcount, endcount)
        if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining value of "//varname_in_ncfile)

    end subroutine

    subroutine getVertDimVals(fileid, nz, varname_in_ncfile, vertdim_indices, vertdim_vals, error_string)
        integer, intent(in) :: fileid, nz
        character(*), intent(in) :: varname_in_ncfile, error_string
        integer, intent(in) :: vertdim_indices(nz)
        real(kind=real_kind), intent(out) :: vertdim_vals(nz)

        integer :: ncerr, varid, vert_index, vert_index_count, startcount(1), endcount(1)
        real(kind=real_kind) :: readval(1)
        
        endcount = (/1/)
        do vert_index_count=1, nz
            vert_index = vertdim_indices(vert_index_count)
            startcount = (/vert_index/)

            ncerr = nf90_inq_varid(fileid, varname_in_ncfile, varid)
            if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)

            ncerr = nf90_get_var(fileid, varid, readval, startcount, endcount)
            if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining value of "//varname_in_ncfile)
            vertdim_vals(vert_index_count) = readval(1)

        end do


    end subroutine

    subroutine getVar2D_real(fileid, varname_in_ncfile, varArr,  error_string )
        integer, intent(in) :: fileid
        character (*), intent(in) :: varname_in_ncfile, error_string
        real (kind=real_kind), intent(out):: varArr(:,:)

        integer :: ncerr, varid

        ncerr = nf90_inq_varid(fileid, varname_in_ncfile, varid)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)
        ncerr = nf90_get_var(fileid, varid, varArr)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)
        print *,'complete reading ', varname_in_ncfile
    end subroutine

    subroutine getVar3DatZlevel_real(fileid, varname_in_ncfile, z_index, time_index, varArr, units, long_name, error_string )
        integer, intent(in) :: fileid, z_index, time_index
        character (*), intent(in) :: varname_in_ncfile, error_string
        character (len=units_len), intent(out) :: units
        character (len=longname_len), intent(out) :: long_name
        real (kind=real_kind), intent(out):: varArr(:,:)

        integer :: ncerr, varid, endcount(4),startcount(4), unitlength, longnamelength

        ! check to see if variable is present
        ncerr = nf90_inq_varid(fileid, varname_in_ncfile, varid)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)

        ! check to see if variable attribute units is present and get units string
        ncerr = nf90_get_att(fileid, varid, "units", units)
        if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining units of "//varname_in_ncfile)

        ! check to see if variable attribute long_name is present and get long_name string
        ncerr = nf90_get_att(fileid, varid, "long_name", long_name)
        if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining long_name of "//varname_in_ncfile)

        startcount = (/1, 1, z_index, time_index/) !x, y, z, time
        endcount = (/nxu, nyu, 1 , 1/)

        ncerr = nf90_get_var(fileid, varid, varArr, startcount, endcount)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)
        print *,'complete reading ', varname_in_ncfile
    end subroutine

    subroutine getVar2D_WithAttrs_real(fileid, varname_in_ncfile, time_index, varArr, units, long_name, error_string )
        integer, intent(in) :: fileid, time_index
        character (*), intent(in) :: varname_in_ncfile, error_string
        character (len=units_len), intent(out) :: units
        character (len=longname_len), intent(out) :: long_name
        real (kind=real_kind), intent(out):: varArr(:,:)

        integer :: ncerr, varid, endcount(3),startcount(3), unitlength, longnamelength

        ! check to see if variable is present
        ncerr = nf90_inq_varid(fileid, varname_in_ncfile, varid)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)

        ! check to see if variable attribute units is present and get units string
        ! ncerr = nf90_inq_att(fileid, varid, "units", &
        !                         len = unitlength)
        ! if ( ncerr .NE. nf90_noerr) call handle_err(ncerr, "error inq units of "//varname_in_ncfile)
        ncerr = nf90_get_att(fileid, varid, "units", units)
        if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining units of "//varname_in_ncfile)

        ! check to see if variable attribute long_name is present and get long_name string
        ! ncerr = nf90_inq_att(fileid, varid, "long_name", &
        !                         len = longnamelength)
        ! if ( ncerr .NE. nf90_noerr) call handle_err(ncerr, "error inq long_name of "//varname_in_ncfile)
        ncerr = nf90_get_att(fileid, varid, "long_name", long_name)
        if (ncerr .NE. nf90_noerr) call handle_err(ncerr, "error obtaining long_name of "//varname_in_ncfile)

        startcount = (/1, 1, time_index/) !x, y, z, time
        endcount = (/nxu, nyu, 1 /)

        ncerr = nf90_get_var(fileid, varid, varArr, startcount, endcount)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, error_string)
        print *,'complete reading ', varname_in_ncfile
    end subroutine


    subroutine defineDimension(file_id, dimlen, dim_id, dim_name, error_string)
        integer, intent(in) :: file_id, dimlen
        integer, intent(out):: dim_id 
        character(len=*), intent(in) :: dim_name, error_string

        integer :: ncerr

        ncerr = nf90_def_dim(file_id, trim(adjustl(dim_name)), dimlen, dim_id)
        if (ncerr /= nf90_noerr) call handle_err(ncerr, error_string)
    end subroutine

    subroutine defineVariables(file_id, varname, vartype, coords_ids, var_id, att_names, att_values )

        integer(kind=int_kind), intent(in) :: file_id, vartype
        character(len=*), intent(in) :: varname
        character(len=longname_len), intent(in) :: att_names(:), att_values(:)

        integer(kind=int_kind), intent(in) :: coords_ids(:)
        integer(kind=int_kind), intent(out) :: var_id

        integer(kind=int_kind) :: num_atts, num_coords, ncerr, counter

        select case (vartype)
        case (1) ! int var type
            ncerr = nf90_def_var(file_id, varname , nf90_int,&
            coords_ids, var_id)
        case (2) ! float var type
            ncerr = nf90_def_var(file_id, varname , nf90_float,&
            coords_ids, var_id)
        case (3) ! double var type
            ncerr = nf90_def_var(file_id, varname , nf90_double,&
            coords_ids, var_id)
        end select
        
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'defining variable '//varname)

        num_atts = size(att_names)
        
        do counter=1, num_atts
    
            ncerr = nf90_put_att(file_id, var_id, trim(adjustl(att_names(counter))), &
                                & trim(adjustl(att_values(counter))))
            if (ncerr /= nf90_noerr)  call handle_err(ncerr, 'defining '//varname//' '//trim(adjustl(att_names(counter))))
    
        end do
    
        
    end subroutine


end module
