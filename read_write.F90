module read_write
    use kinds
    use ncdf_wrapper
    use netcdf_io
    use fields
    use input_data_info
    use gridModule
    use filterparallel
    use constants
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

        ULONG = ULONG * 3.141592653589793/180.0
        ULAT = ULAT * 3.141592653589793/180.0

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

        if (taskid .EQ. MASTER) then
            call getTimeVar(file_id, trim(adjustl(timevar_name)), &
                                                  timevar_units, timevar_calendar, time_index, &
                                                  timevar_val, "error getting time variable info "//timevar_name)
            call getVertDimVals(file_id, nzu, trim(adjustl(vertdim_name)), arr_z_index, vertdim_vals, "error getting vertical dim ")

        endif

        do field_count =1, num_scalar2D_fields
            call getVar2D_WithAttrs_real(file_id, trim(adjustl(scalar2D_fields_info(field_count)%varname)), &
                        &       time_index, &
                        &       scalar2D_fields(:,:, field_count), &
                        &       scalar2D_fields_info(field_count)%units, &
                        &       scalar2D_fields_info(field_count)%long_name, &
                        &       "getVar2D_real" )
            where (abs(scalar2D_fields(:,:,field_count)) > 1d20)
                scalar2D_fields(:,:, field_count) = 0
            end where
        end do

        do field_count =1, num_scalar3D_fields
            do z_count =1, num_zlevels
                z_index = arr_z_index(z_count)
                call getVar3DatZlevel_real(file_id, trim(adjustl(scalar3D_fields_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       scalar3D_fields(:,:,z_count, field_count), &
                                    &       scalar3D_fields_info(field_count)%units, &
                                    &       scalar3D_fields_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real" )
                where (abs(scalar3D_fields(:,:,z_count, field_count)) > 1d20)
                    scalar3D_fields(:,:,z_count, field_count) = 0
                end where
            end do
        end do

        do field_count=1, num_vector2D_fields
            call getVar2D_WithAttrs_real(file_id, trim(adjustl(vector2DX_fields_info(field_count)%varname)), &
                        &       time_index, &
                        &       vector2DX_fields(:,:, field_count), &
                        &       vector2DX_fields_info(field_count)%units, &
                        &       vector2DX_fields_info(field_count)%long_name, &
                        &       "getVar2D_real vector2dx" )

            call getVar2D_WithAttrs_real(file_id, trim(adjustl(vector2DY_fields_info(field_count)%varname)), &
                        &       time_index, &
                        &       vector2DY_fields(:,:, field_count), &
                        &       vector2DY_fields_info(field_count)%units, &
                        &       vector2DY_fields_info(field_count)%long_name, &
                        &       "getVar2D_real vector2dy" )

            where (abs(vector2DX_fields(:,:, field_count)) > 1d20)
                vector2DX_fields(:,:, field_count) = 0
            end where

            where (abs(vector2DY_fields(:,:, field_count)) > 1d20)
                vector2DY_fields(:,:, field_count) = 0
            end where
        end do

        do field_count=1, num_vector3D_fields
            do z_count =1, num_zlevels
                z_index = arr_z_index(z_count)
                call getVar3DatZlevel_real(file_id, trim(adjustl(vector3DX_fields_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       vector3DX_fields(:,:,z_count, field_count), &
                                    &       vector3DX_fields_info(field_count)%units, &
                                    &       vector3DX_fields_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real vector3dx" )

                call getVar3DatZlevel_real(file_id, trim(adjustl(vector3DY_fields_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       vector3DY_fields(:,:,z_count, field_count), &
                                    &       vector3DY_fields_info(field_count)%units, &
                                    &       vector3DY_fields_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real vector3dy" )

                call getVar3DatZlevel_real(file_id, trim(adjustl(vector3DZ_fields_info(field_count)%varname)), &
                                    &       z_index, time_index, &
                                    &       vector3DZ_fields(:,:,z_count, field_count), &
                                    &       vector3DZ_fields_info(field_count)%units, &
                                    &       vector3DZ_fields_info(field_count)%long_name, &
                                    &       "getVar3DatZlevel_real vector3dz" )

                where (abs(vector3DX_fields(:,:, z_count, field_count)) > 1d20)
                    vector3DX_fields(:,:, z_count, field_count) = 0
                end where
                where (abs(vector3DY_fields(:,:, z_count, field_count)) > 1d20)
                    vector3DY_fields(:,:, z_count, field_count) = 0
                end where
                where (abs(vector3DZ_fields(:,:, z_count, field_count)) > 1d20)
                    vector3DZ_fields(:,:, z_count, field_count) = 0
                end where
            end do
        end do

        ncerr = nf90_close(file_id)
        print *, 'closed file  ', trim(adjustl(filename)), ' ... '
        
    end subroutine

    subroutine writeDirectFilteredFields(fullfilename, x_dimname, y_dimname, z_dimname, filter_dimname, time_dimname)
        character(len=*) , intent(in) :: fullfilename, x_dimname, y_dimname, z_dimname, filter_dimname, time_dimname
        integer :: file_id, xdim_id, ydim_id, zdim_id, filterdim_id, timedim_id, coords_2d(4), coords_3d(5), ncerr , &
                   timevar_id, latvar_id, lonvar_id, zvar_id, filtervar_id, field_count, numvars, var_index, counter
    
        character(len=longname_len) :: att_names(2), att_values(2)
        character(len=varname_len) :: varname

        integer(kind=int_kind), allocatable :: varids(:)

        real(kind=real_kind) , allocatable:: dummy2d(:,:,:,:), dummy3d(:,:,:,:,:)  ! x, y, z, filterlength, time

        integer, parameter :: INT_TYPE = 1, FLOAT_TYPE = 2, DOUBLE_TYPE=3


        numvars = num_scalar2D_fields + num_scalar3D_fields + 2*num_vector2D_fields + 3*num_vector3D_fields

        

        allocate(varids(numvars))
        allocate(dummy2d(nxu,nyu,1,1), dummy3d(nxu,nyu,nzu,1,1))
    
        print *, 'writing file ', trim(adjustl(fullfilename))
        
        !-------------------------------------------------------------------
        !  open netcdf file
        !-------------------------------------------------------------------
    
        ncerr = nf90_create(fullfilename, nf90_clobber, file_id)
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'nf90_open to write')
    
        !-------------------------------------------------------------------
        !  define dimensions
        !-------------------------------------------------------------------
    
        call defineDimension(file_id, nxu, xdim_id, x_dimname, 'error in defining xdim')
        call defineDimension(file_id, nyu, ydim_id, y_dimname, 'error in defining ydim')
        call defineDimension(file_id, nzu, zdim_id, z_dimname, 'error in defining zdim')
        call defineDimension(file_id, num_filterlengths, filterdim_id, filter_dimname, 'error in defining filterdim')
        !call defineDimension(file_id, 1, timedim_id, timedim_name, 'error in defining timedim')
    
        ncerr = nf90_def_dim(file_id, trim(adjustl(time_dimname)), 1, timedim_id) !NF90_UNLIMITED
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'error in defining time dimension')
    
        !-------------------------------------------------------------------
        !  define coordinates
        !-------------------------------------------------------------------
        !-- time
    
        att_names(1) = 'units'
        att_values(1) = timevar_units
    
        att_names(2) = 'calendar'
        att_values(2) = timevar_calendar
        call defineVariables(file_id, time_dimname, FLOAT_TYPE, (/timedim_id/), timevar_id, att_names, att_values )
    
        !-------------------------------------------------------------------
        !-- filterlength
        att_names(1) = 'units'
        att_values(1) = 'kilometers'
    
        att_names(2) = 'long_name'
        att_values(2) = 'coarse graining filtering isotropic xy kernel'
        call defineVariables(file_id, filter_dimname, FLOAT_TYPE, (/filterdim_id/), filtervar_id, att_names, att_values )
    
        !-------------------------------------------------------------------
        !-- vertical level
    
        att_names(1) = 'units'
        att_values(1) = ' '
    
        att_names(2) = 'long_name'
        att_values(2) = 'ROMS vertical co-ordinate(s-rho), 0 at surface -1 at bottom'
        call defineVariables(file_id, z_dimname, FLOAT_TYPE, (/zdim_id/), zvar_id, att_names, att_values )
    
    
        coords_2d(1)=xdim_id
        coords_2d(2)=ydim_id
        coords_2d(3)=filterdim_id
        coords_2d(4)=timedim_id
    
    
        coords_3d(1)=xdim_id
        coords_3d(2)=ydim_id
        coords_3d(3)=zdim_id
        coords_3d(4)=filterdim_id
        coords_3d(5)=timedim_id


        var_index = 1
        do field_count =1, num_scalar2D_fields
            varname = trim(adjustl(scalar2D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(scalar2D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(scalar2D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
            
        end do

        do field_count =1, num_scalar3D_fields
            varname = trim(adjustl(scalar3D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(scalar3D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(scalar3D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
            
        end do

        do field_count=1, num_vector2D_fields
            varname = trim(adjustl(vector2DX_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DX_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DX_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DY_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DY_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DY_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
        end do

        do field_count=1, num_vector3D_fields
            varname = trim(adjustl(vector3DX_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DX_fields_info(field_count)%units))    
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DX_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DY_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DY_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DY_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1


            varname = trim(adjustl(vector3DZ_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DZ_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DZ_fields_info(field_count)%long_name))
            
            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
            
        end do

        ncerr = nf90_enddef(file_id)
        if (ncerr /= nf90_noerr) stop 'at enddef'

        ! Put time, lengthscale and z co-ordinates value

        do counter = 1, num_zlevels
            ncerr = nf90_put_var(file_id, zvar_id,(/arr_z_index(counter)/),       &
                      start = (/counter/), &
                      count = (/1/))
            if(ncerr /= nf90_noerr) call handle_err(ncerr, 'writing z co-ordinate vals')
        end do

        do counter = 1, num_filterlengths
            ncerr = nf90_put_var(file_id, filtervar_id,(/arr_filterlengths(counter)/),       &
                      start = (/counter/), &
                      count = (/1/))
            if(ncerr /= nf90_noerr) call handle_err(ncerr, 'writing filterlength co-ordinate vals')
        end do

        ncerr = nf90_put_var(file_id, timevar_id, (/timevar_val(1)/),       &
                 start = (/1/), &
                 count = (/1/))
        if(ncerr /= nf90_noerr) call handle_err(ncerr, 'time co-ordinate vals')


        ! Start writing variables
        do counter = 1, num_filterlengths
            var_index = 1
            do field_count =1, num_scalar2D_fields
                dummy2d(:,:, 1,1) = OL_scalar2D_fields(:,:, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1    
            end do

            do field_count =1, num_scalar3D_fields
                dummy3d(:,:, :, 1,1) = OL_scalar3D_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1    
            end do

            do field_count=1, num_vector2D_fields
                dummy2d(:, :, 1,1) = OL_vector2DX_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1

                dummy2d(:, :, 1,1) = OL_vector2DY_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1
            end do

            do field_count=1, num_vector3D_fields
                dummy3d(:,:,:, 1,1) = OL_vector3DX_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1

                dummy3d(:,:,:, 1,1) = OL_vector3DY_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1
                
                dummy3d(:,:,:, 1,1) = OL_vector3DZ_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1
                
            end do
        end do

        ncerr = nf90_close(file_id)
        if (ncerr /= nf90_noerr) stop 'at close'
    
    
    end subroutine

    subroutine writeUnfiltHelmHoltzDeompFields(fullfilename, x_dimname, y_dimname, z_dimname, time_dimname)
        character(len=*) , intent(in) :: fullfilename, x_dimname, y_dimname, z_dimname, time_dimname
        integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id, coords_2d(3), coords_3d(4), ncerr , &
                   timevar_id, latvar_id, lonvar_id, zvar_id, field_count, numvars, var_index, counter
    
        character(len=longname_len) :: att_names(2), att_values(2)
        character(len=varname_len) :: varname

        integer(kind=int_kind), allocatable :: varids(:)

        real(kind=real_kind) , allocatable:: dummy2d(:,:,:), & ! x, y, time
                                             dummy3d(:,:,:,:)  ! x, y, z, time

        integer, parameter :: INT_TYPE = 1, FLOAT_TYPE = 2, DOUBLE_TYPE=3

        numvars = 6*num_vector2D_fields + &  ! phi, psi, phi_x, phi_y, psi_x, psi_y
                  6*num_vector3D_fields      ! phi, psi, phi_x, phi_y, psi_x, psi_y

        allocate(varids(numvars))
        allocate(dummy2d(nxu, nyu, 1), dummy3d(nxu, nyu, nzu,1))
    
        print *, 'writing file ', trim(adjustl(fullfilename))
        
        !-------------------------------------------------------------------
        !  open netcdf file
        !-------------------------------------------------------------------
    
        ncerr = nf90_create(fullfilename, nf90_clobber, file_id)
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'nf90_open to write')
    
        !-------------------------------------------------------------------
        !  define dimensions
        !-------------------------------------------------------------------
    
        call defineDimension(file_id, nxu, xdim_id, x_dimname, 'error in defining xdim')
        call defineDimension(file_id, nyu, ydim_id, y_dimname, 'error in defining ydim')
        call defineDimension(file_id, nzu, zdim_id, z_dimname, 'error in defining zdim')
        
        ncerr = nf90_def_dim(file_id, trim(adjustl(time_dimname)), 1, timedim_id) !NF90_UNLIMITED
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'error in defining time dimension')
    
        !-------------------------------------------------------------------
        !  define coordinates
        !-------------------------------------------------------------------
        !-- time
    
        att_names(1) = 'units'
        att_values(1) = timevar_units
    
        att_names(2) = 'calendar'
        att_values(2) = timevar_calendar
        call defineVariables(file_id, time_dimname, FLOAT_TYPE, (/timedim_id/), timevar_id, att_names, att_values )

    
        !-------------------------------------------------------------------
        !-- vertical level
    
        att_names(1) = 'units'
        att_values(1) = ' '
    
        att_names(2) = 'long_name'
        att_values(2) = 'ROMS vertical co-ordinate(s-rho), 0 at surface -1 at bottom'
        call defineVariables(file_id, z_dimname, FLOAT_TYPE, (/zdim_id/), zvar_id, att_names, att_values )
    
    
        coords_2d(1)=xdim_id
        coords_2d(2)=ydim_id
        coords_2d(3)=timedim_id
    
    
        coords_3d(1)=xdim_id
        coords_3d(2)=ydim_id
        coords_3d(3)=zdim_id
        coords_3d(4)=timedim_id


        var_index = 1
        do field_count=1, num_vector2D_fields

            varname = trim(adjustl(phi2D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(phi2D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(phi2D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(psi2D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(psi2D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(psi2D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1


            varname = trim(adjustl(vector2DX_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DX_phi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DX_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DY_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DY_phi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DY_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DX_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DX_psi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DX_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DY_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DY_psi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DY_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
        end do

        do field_count=1, num_vector3D_fields

            varname = trim(adjustl(phi3D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(phi3D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(phi3D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(psi3D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(psi3D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(psi3D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1


            varname = trim(adjustl(vector3DX_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DX_phi_fields_info(field_count)%units))    
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DX_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DY_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DY_phi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DY_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DX_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DX_psi_fields_info(field_count)%units))    
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DX_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DY_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DY_psi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DY_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
            
        end do

        ncerr = nf90_enddef(file_id)
        if (ncerr /= nf90_noerr) stop 'at enddef'

        ! Put time, lengthscale and z co-ordinates value

        do counter = 1, num_zlevels
            ncerr = nf90_put_var(file_id, zvar_id,(/arr_z_index(counter)/),       &
                      start = (/counter/), &
                      count = (/1/))
            if(ncerr /= nf90_noerr) call handle_err(ncerr, 'writing z co-ordinate vals')
        end do

        ncerr = nf90_put_var(file_id, timevar_id, (/timevar_val(1)/),       &
                 start = (/1/), &
                 count = (/1/))
        if(ncerr /= nf90_noerr) call handle_err(ncerr, 'time co-ordinate vals')

        ! Start writing variables
        var_index = 1

        do field_count=1, num_vector2D_fields
            dummy2d(:, :, 1) = phi2D_fields(:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                    start = (/1, 1, 1/), &
                    count = (/nxu, nyu, 1 /))
            var_index = var_index + 1

            dummy2d(:, :, 1) = psi2D_fields(:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                    start = (/1, 1, 1/), &
                    count = (/nxu, nyu, 1 /))
            var_index = var_index + 1

            dummy2d(:, :, 1) = vector2DX_phi_fields(:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                    start = (/1, 1, 1/), &
                    count = (/nxu, nyu, 1 /))
            var_index = var_index + 1

            dummy2d(:, :, 1) = vector2DY_phi_fields(:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                    start = (/1, 1, 1/), &
                    count = (/nxu, nyu, 1 /))
            var_index = var_index + 1

            dummy2d(:, :, 1) = vector2DX_psi_fields(:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                    start = (/1, 1, 1/), &
                    count = (/nxu, nyu, 1 /))
            var_index = var_index + 1

            dummy2d(:, :, 1) = vector2DY_psi_fields(:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                    start = (/1, 1, 1/), &
                    count = (/nxu, nyu, 1 /))
            var_index = var_index + 1
        end do

        do field_count=1, num_vector3D_fields
            dummy3d(:,:,:, 1) = phi3D_fields(:,:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                    start = (/1, 1, 1, 1/), &
                    count = (/nxu, nyu, nzu, 1 /))
            var_index = var_index + 1

            dummy3d(:,:,:, 1) = psi3D_fields(:,:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                    start = (/1, 1, 1, 1/), &
                    count = (/nxu, nyu, nzu, 1 /))
            var_index = var_index + 1

            dummy3d(:,:,:, 1) = vector3DX_phi_fields(:,:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                    start = (/1, 1, 1, 1/), &
                    count = (/nxu, nyu, nzu, 1 /))
            var_index = var_index + 1

            dummy3d(:,:,:, 1) = vector3DY_phi_fields(:,:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                    start = (/1, 1, 1, 1/), &
                    count = (/nxu, nyu, nzu, 1 /))
            var_index = var_index + 1

            dummy3d(:,:,:, 1) = vector3DX_psi_fields(:,:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                    start = (/1, 1, 1, 1/), &
                    count = (/nxu, nyu, nzu, 1 /))
            var_index = var_index + 1

            dummy3d(:,:,:, 1) = vector3DY_psi_fields(:,:, :, field_count)
            ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                    start = (/1, 1, 1, 1/), &
                    count = (/nxu, nyu, nzu, 1 /))
            var_index = var_index + 1
        end do

        ncerr = nf90_close(file_id)
        if (ncerr /= nf90_noerr) stop 'at close'

    end subroutine

    subroutine writeCoarseGrainedFields(fullfilename, x_dimname, y_dimname, z_dimname, filter_dimname, time_dimname)
        character(len=*) , intent(in) :: fullfilename, x_dimname, y_dimname, z_dimname, filter_dimname, time_dimname
        integer :: file_id, xdim_id, ydim_id, zdim_id, filterdim_id, timedim_id, coords_2d(4), coords_3d(5), ncerr , &
                   timevar_id, latvar_id, lonvar_id, zvar_id, filtervar_id, field_count, numvars, var_index, counter
    
        character(len=longname_len) :: att_names(2), att_values(2)
        character(len=varname_len) :: varname

        integer(kind=int_kind), allocatable :: varids(:)

        real(kind=real_kind) , allocatable:: dummy2d(:,:,:,:), dummy3d(:,:,:,:,:)  ! x, y, z, filterlength, time

        integer, parameter :: INT_TYPE = 1, FLOAT_TYPE = 2, DOUBLE_TYPE=3



        numvars = num_scalar2D_fields + num_scalar3D_fields + &
                  6*num_vector2D_fields + &  ! OL_psi, OL_phi, OL_phi_u, OL_phi_v , OL_psi_u, OL_psi_v 
                  7*num_vector3D_fields      ! OL_psi, OL_phi, OL_phi_u, OL_phi_v , OL_psi_u, OL_psi_v , OL_w

        

        allocate(varids(numvars))
        allocate(dummy2d(nxu,nyu,1,1), dummy3d(nxu,nyu,nzu,1,1))
    
        print *, 'writing file ', trim(adjustl(fullfilename))
        
        !-------------------------------------------------------------------
        !  open netcdf file
        !-------------------------------------------------------------------
    
        ncerr = nf90_create(fullfilename, nf90_clobber, file_id)
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'nf90_open to write')
    
        !-------------------------------------------------------------------
        !  define dimensions
        !-------------------------------------------------------------------
    
        call defineDimension(file_id, nxu, xdim_id, x_dimname, 'error in defining xdim')
        call defineDimension(file_id, nyu, ydim_id, y_dimname, 'error in defining ydim')
        call defineDimension(file_id, nzu, zdim_id, z_dimname, 'error in defining zdim')
        call defineDimension(file_id, num_filterlengths, filterdim_id, filter_dimname, 'error in defining filterdim')
        !call defineDimension(file_id, 1, timedim_id, timedim_name, 'error in defining timedim')
    
        ncerr = nf90_def_dim(file_id, trim(adjustl(time_dimname)), 1, timedim_id) !NF90_UNLIMITED
        if (ncerr /= nf90_noerr) call handle_err(ncerr, 'error in defining time dimension')
    
        !-------------------------------------------------------------------
        !  define coordinates
        !-------------------------------------------------------------------
        !-- time
    
        att_names(1) = 'units'
        att_values(1) = timevar_units
    
        att_names(2) = 'calendar'
        att_values(2) = timevar_calendar
        call defineVariables(file_id, time_dimname, FLOAT_TYPE, (/timedim_id/), timevar_id, att_names, att_values )
    
        !-------------------------------------------------------------------
        !-- filterlength
        att_names(1) = 'units'
        att_values(1) = 'kilometers'
    
        att_names(2) = 'long_name'
        att_values(2) = 'coarse graining filtering isotropic xy kernel'
        call defineVariables(file_id, filter_dimname, FLOAT_TYPE, (/filterdim_id/), filtervar_id, att_names, att_values )
    
        !-------------------------------------------------------------------
        !-- vertical level
    
        att_names(1) = 'units'
        att_values(1) = ' '
    
        att_names(2) = 'long_name'
        att_values(2) = 'ROMS vertical co-ordinate(s-rho), 0 at surface -1 at bottom'
        call defineVariables(file_id, z_dimname, FLOAT_TYPE, (/zdim_id/), zvar_id, att_names, att_values )
    
    
        coords_2d(1)=xdim_id
        coords_2d(2)=ydim_id
        coords_2d(3)=filterdim_id
        coords_2d(4)=timedim_id
    
    
        coords_3d(1)=xdim_id
        coords_3d(2)=ydim_id
        coords_3d(3)=zdim_id
        coords_3d(4)=filterdim_id
        coords_3d(5)=timedim_id


        var_index = 1
        do field_count =1, num_scalar2D_fields
            varname = trim(adjustl(scalar2D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(scalar2D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(scalar2D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
            
        end do

        do field_count =1, num_scalar3D_fields
            varname = trim(adjustl(scalar3D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(scalar3D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(scalar3D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
            
        end do

        do field_count=1, num_vector2D_fields
            varname = trim(adjustl(phi2D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(phi2D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(phi2D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(psi2D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(psi2D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(psi2D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DX_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DX_phi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DX_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DY_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DY_phi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DY_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DX_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DX_psi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DX_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector2DY_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector2DY_psi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector2DY_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, FLOAT_TYPE, coords_2d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
        end do

        do field_count=1, num_vector3D_fields
            varname = trim(adjustl(phi3D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(phi3D_fields_info(field_count)%units))    
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(phi3D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(psi3D_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(psi3D_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(psi3D_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DX_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DX_phi_fields_info(field_count)%units))    
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DX_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DY_phi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DY_phi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DY_phi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DX_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DX_psi_fields_info(field_count)%units))    
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DX_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1

            varname = trim(adjustl(vector3DY_psi_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DY_psi_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DY_psi_fields_info(field_count)%long_name))

            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1


            varname = trim(adjustl(vector3DZ_fields_info(field_count)%varname))
            att_names(1) = 'units'
            att_values(1) = trim(adjustl(vector3DZ_fields_info(field_count)%units))
            att_names(2) = 'long_name'
            att_values(2) = trim(adjustl(vector3DZ_fields_info(field_count)%long_name))
            
            call defineVariables(file_id, varname, 2, coords_3d, varids(var_index), att_names, att_values )
            var_index = var_index + 1
            
        end do

        ncerr = nf90_enddef(file_id)
        if (ncerr /= nf90_noerr) stop 'at enddef'

        ! Put time, lengthscale and z co-ordinates value

        do counter = 1, num_zlevels
            ncerr = nf90_put_var(file_id, zvar_id,(/arr_z_index(counter)/),       &
                      start = (/counter/), &
                      count = (/1/))
            if(ncerr /= nf90_noerr) call handle_err(ncerr, 'writing z co-ordinate vals')
        end do

        do counter = 1, num_filterlengths
            ncerr = nf90_put_var(file_id, filtervar_id,(/arr_filterlengths(counter)/),       &
                      start = (/counter/), &
                      count = (/1/))
            if(ncerr /= nf90_noerr) call handle_err(ncerr, 'writing filterlength co-ordinate vals')
        end do

        ncerr = nf90_put_var(file_id, timevar_id, (/timevar_val(1)/),       &
                 start = (/1/), &
                 count = (/1/))
        if(ncerr /= nf90_noerr) call handle_err(ncerr, 'time co-ordinate vals')


        ! Start writing variables
        do counter = 1, num_filterlengths
            var_index = 1
            do field_count =1, num_scalar2D_fields
                dummy2d(:,:, 1,1) = OL_scalar2D_fields(:,:, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1    
            end do

            do field_count =1, num_scalar3D_fields
                dummy3d(:,:, :, 1,1) = OL_scalar3D_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1    
            end do

            do field_count=1, num_vector2D_fields
                dummy2d(:, :, 1,1) = OL_phi2D_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1

                dummy2d(:, :, 1,1) = OL_psi2D_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1

                dummy2d(:, :, 1,1) = OL_vector2DX_phi_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1

                dummy2d(:, :, 1,1) = OL_vector2DY_phi_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1

                dummy2d(:, :, 1,1) = OL_vector2DX_psi_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1

                dummy2d(:, :, 1,1) = OL_vector2DY_psi_fields(:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy2d,       &
                        start = (/1, 1, counter, 1/), &
                        count = (/nxu, nyu, 1, 1 /))
                var_index = var_index + 1
            end do

            do field_count=1, num_vector3D_fields
                dummy3d(:,:,:, 1,1) = OL_phi3D_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1

                dummy3d(:,:,:, 1,1) = OL_psi3D_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1

                dummy3d(:,:,:, 1,1) = OL_vector3DX_phi_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1

                dummy3d(:,:,:, 1,1) = OL_vector3DY_phi_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1

                dummy3d(:,:,:, 1,1) = OL_vector3DX_psi_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1

                dummy3d(:,:,:, 1,1) = OL_vector3DY_psi_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1
                
                dummy3d(:,:,:, 1,1) = OL_vector3DZ_fields(:,:, :, field_count, counter)
                ncerr = nf90_put_var(file_id, varids(var_index), dummy3d,       &
                        start = (/1, 1, 1, counter, 1/), &
                        count = (/nxu, nyu, nzu, 1, 1 /))
                var_index = var_index + 1
                
            end do
        end do

        ncerr = nf90_close(file_id)
        if (ncerr /= nf90_noerr) stop 'at close'
    
    
    end subroutine


   


end module