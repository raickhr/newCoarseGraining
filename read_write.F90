module read_write
    use kinds
    use ncdf_wrapper
    use netcdf_io
    use fields
    use input_data_info
    use gridModule
    implicit none

    contains

    subroutine get_grid_nc(path, grid_filename)

        character(len=*) :: path
        character(len=*) :: grid_filename

        ! local variables
        character(len=fullfilename_len) :: file_netcdf

        integer :: ncerr, file_id

        !---- Open file
        file_netcdf = trim(adjustl(path)) // '/' // trim(adjustl(grid_filename))
        print *,'Opening grid file : ',trim(adjustl(file_netcdf))

        ncerr = nf90_open(trim(adjustl(file_netcdf)), nf90_nowrite, file_id)
        if ( ncerr .NE. nf90_noerr )  call handle_err(ncerr, 'nf90_open')

        print *,'File ',trim(adjustl(file_netcdf)), ' opened ..'
        print *,' '

        print *,'Reading file now ..'
        
        call getVar2D_real(file_id, 'DXU', DXU,  'error reading DXU')
        call getVar2D_real(file_id, 'DYU', DYU,  'error reading DYU')
        call getVar2D_real(file_id, 'ULAT', ULAT,  'error reading ULAT')
        call getVar2D_real(file_id, 'ULONG', ULONG,  'error reading ULONG')
        call getVar2D_real(file_id, 'UAREA', UAREA,  'error reading UAREA')
        call getVar2D_real(file_id, 'KMU', KMU,  'error reading KMU')
        call getVar2D_real(file_id, 'HU', HU,  'error reading HU')
        call getVar2D_real(file_id, 'FCORU', FCORU,  'error reading FCORU')

        ! closing the grid file
        ncerr = nf90_close(file_id)

        print *,'Read grid file complete ..'
        print *,' '

        
    end subroutine get_grid_nc

    subroutine read_fields(filename, &
                           time_index)

        character(len=*) , intent(in) :: filename
        integer(kind=int_kind) :: time_index

        integer(kind=int_kind) :: f_error, ierr, file_id, ncerr, nfields, nzlevels, &
                                  z_index, field_count, z_count, field_count3, fieldIndex

        print *,''
        print *, 'Opening file  ', trim(adjustl(filename)), ' ... '

        ierr = nf90_open(trim(adjustl(filename)), nf90_nowrite, file_id)
        if ( ierr .NE. nf90_noerr )  call handle_err(ierr, 'nf90_open')

        print *,''
        print *, 'Opened  ', trim(adjustl(filename))

        do field_count =1, num_scalar_fields
            do z_count =1, num_zlevels
                z_index = arr_z_index(z_count)
                call getVar3DatZlevel_real(file_id, trim(adjustl(scalar_field_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       scalar_fields(:,:,z_count, field_count), &
                                    &       scalar_field_info(field_count)%units, &
                                    &       scalar_field_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real" )
            end do
        end do

        do field_count=1, num_2Dvector_fields
            call getVar2D_WithAttrs_real(file_id,  trim(adjustl(vector2DX_field_info(field_count)%varname)), &
                                  &       time_index, &
                                  &       vector2DX_fields(:,:,field_count), &
                                  &       vector2DX_field_info(field_count)%units, &
                                  &       vector2DX_field_info(field_count)%long_name, &
                                  &       "getVar2D_withAttrs_real vector2dx" )

            call getVar2D_WithAttrs_real(file_id,  trim(adjustl(vector2DX_field_info(field_count)%varname)), &
                                  &       time_index, &
                                  &       vector2DY_fields(:,:,field_count), &
                                  &       vector2DY_field_info(field_count)%units, &
                                  &       vector2DY_field_info(field_count)%long_name, &
                                  &       "getVar2D_withAttrs_real vector2dy" )
        end do

        do field_count=1, num_3Dvector_fields
            do z_count =1, num_zlevels
                z_index = arr_z_index(z_count)
                call getVar3DatZlevel_real(file_id, trim(adjustl(vector3DX_field_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       vector3DX_fields(:,:,z_count, field_count), &
                                    &       vector3DX_field_info(field_count)%units, &
                                    &       vector3DX_field_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real vector3dx" )

                call getVar3DatZlevel_real(file_id, trim(adjustl(vector3DY_field_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       vector3DY_fields(:,:,z_count, field_count), &
                                    &       vector3DY_field_info(field_count)%units, &
                                    &       vector3DY_field_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real vector3dz" )

                call getVar3DatZlevel_real(file_id, trim(adjustl(vector3DZ_field_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       vector3DZ_fields(:,:,z_count, field_count), &
                                    &       vector3DZ_field_info(field_count)%units, &
                                    &       vector3DZ_field_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real vector3dz" )
            end do
        end do

        ncerr = nf90_close(file_id)
        print *, 'closed file  ', trim(adjustl(filename)), ' ... '
        
    end subroutine

end module