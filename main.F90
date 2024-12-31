program main
    use configurationMod
    use gridModule
    use mpiwrapper
    use fields
    use filterparallel
    use input_data_info
    use read_write
    use multiGridHelmHoltz

    implicit none
    
    type(configuration):: config         ! object that defines run configuration
    integer :: file_index, time_index, z_index, facList(4)
    character (len=filename_len) :: writefilename


    facList = (/9, 5, 3, 2/)

    call startMPI()

    if (taskid .EQ. MASTER) then
        call config%construct()   ! Run the constructor for configuration object which also reads the configuration file
        call set_size(config%nx, config%ny, config%nz)
        call init_grid(config%nx, config%ny, config%nz)
        call get_grid_nc(config%InputPath, config%gridFile)

        ! setting field info 
        call set_num_scalar_fields(config%num_of_scalar_fields_to_read)
        call set_num_2Dvector_fields(config%num_of_2Dvector_fields_to_read)
        call set_num_3Dvector_fields(config%num_of_3Dvector_fields_to_read)

        call allocate_scalar_fields(nxu, nyu, nzu)
        call allocate_vector2D_fields(nxu, nyu, nzu)
        call allocate_phi_psi_fields(nxu, nyu, nzu)
        call allocate_vector3D_fields(nxu, nyu, nzu)

        call set_fieldnames(config%list_scalar_fieldsNames, &
                        &   config%list_2DvectorX_fieldsNames,  config%list_2DvectorY_fieldsNames, &
                        &   config%list_3DvectorX_fieldsNames,  config%list_3DvectorY_fieldsNames, config%list_3DvectorZ_fieldsNames)

        ! setting filterlength info
        call set_num_filterlengths(config%nfilter)
        call set_arr_filterlength(config%list_filterLength)

        ! setting input data info for reading
        call set_numfiles_numzlevels_ntimesinafile(config%num_of_files_to_read, &
                                            &      config%nz, &
                                            &      config%startTimeIndex,  &    
                                            &      config%endTimeIndex)

        timevar_name = trim(adjustl(config%timevar_name))
        vertdim_name = trim(adjustl(config%vertdim_name))

        call alloc_arr_z_index()

        call set_arr_z_index(config%list_zlevels)
                                    

    endif

    call broadCastGridInfo()
    call broadCastFieldInfo()
    call broadCastFilterInfo()
    call broadCastInputDataInfo()
    if (taskid .EQ. MASTER) print *, ' all grid info broadcasted !'

    !call 

    call setMultiGrid(facList, ULAT, ULONG, DXU, DYU, UAREA)

    call dividework()

    if (taskid .EQ. MASTER) print *, 'work division done !'

    call allocate_output_fields(nxu, nyu, nzu, num_filterlengths)

    do file_index = 1, num_files
        if (taskid .EQ. MASTER) then
            print *,"AT FILE NUMBER ",file_index,' OF ', num_files
        endif
        do time_index = start_timeindex, end_timeindex
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
            if (taskid .EQ. MASTER) then
                print *,"AT TIME INDEX ",time_index-start_timeindex+1,' OF ', end_timeindex - start_timeindex +1

                call read_fields(trim(adjustl(config%InputPath))//'/'//trim(adjustl(config%list_filenames(file_index))), &
                                    time_index)
                ! call assign_fields()
            endif

            call broadCastReadFields()


            if (taskid .EQ. MASTER) then 
                print *, ''
                print *, 'all read field info at one time instant broadcasted !'
                print *, ''
                print *, ''
                print *, ''
            end if

            call helmholtzDecompAllVecFields()

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            if (taskid == 0 ) print *, 'completed Helmholtz decomposition'

            if (taskid == 0 ) call writePolTorVelWithPsiPhi()

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call filter_allvars()

            ! OL_scalar_fields(:,:,1,:,1) = 1.0d0
            ! OL_scalar_fields(:,:,2,:,1) = 2.0d0

            ! OL_scalar_fields(:,:,1,:,2) = 10*2.0d0
            ! OL_scalar_fields(:,:,2,:,2) = 20*2.0d0

            if (taskid .EQ. MASTER) then 
                print *, ''
                print *, 'filtering all variables completed !'
                print *, ''
            end if

            call collectFilteredFields(arr_numcols_inallprocs, arr_startcolindex_inallprocs, num_filterlengths)

            if (taskid .EQ. MASTER) then 
                print *, ''
                print *, 'all filtered variables collected !'
                print *, ''
            end if
            
            if (taskid .EQ. MASTER) then
                WRITE(writefilename, "(A5,I0.3,A5,I0.3,A3)") "_file", file_index, "_time", time_index, ".nc" 
                call writeFields( trim(adjustl(config%OutputPath))//'/'//writefilename, &
                &                 'xi_rho', 'eta_rho', trim(adjustl(config%vertdim_name)), &
                &                 'Lengthscale', trim(adjustl(config%timevar_name)))
            end if

        end do !close time loop
    end do ! close

    call delMultiGrid()
    
    call deallocate_gridVars()
    call deallocate_filtervars()
    call deallocate_inputData_info()
    call delloacate_fields()
    call config%destruct()

    call MPI_Barrier(MPI_COMM_WORLD, i_err)
    if (taskid .EQ. MASTER) print *, 'Ending program'
    call stopMPI()

end program
