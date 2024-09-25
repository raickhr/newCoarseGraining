program main
    use configurationMod
    use gridModule
    use mpiwrapper
    use fields
    use filterparallel

    implicit none
    
    type(configuration):: config         ! object that defines run configuration
    

    call startMPI()

    if (taskid .EQ. MASTER) then
        config = configuration()   ! Run the constructor for configuration object which also reads the configuration file
        call set_size(config%nx, config%ny, config%nz)
        call init_grid(config%nx, config%ny, config%nz)
        call get_grid_nc(config%InputPath, config%gridFile)

        ! setting field info 
        call set_num_scalar_fields(config%num_of_scalar_fields_to_read)
        call set_num_2Dvector_fields(config%num_of_2Dvector_fields_to_read)
        call set_num_3Dvector_fields(config%num_of_3Dvector_fields_to_read)

        ! setting filterlength info
        call set_num_filterlengths(config%nfilter)
        call set_arr_filterlength(config%list_filterLength)
    endif

    call MPI_BCAST(nxu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(nyu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(nzu, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

    call MPI_Barrier(MPI_COMM_WORLD, i_err)
    call init_grid(nxu, nyu, nzu)

    call MPI_BCAST(DXU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(DYU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(ULAT, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(ULONG, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(KMU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(UAREA, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(HU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(FCORU, nxu*nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

    call MPI_Barrier(MPI_COMM_WORLD, i_err)

    call MPI_BCAST(num_scalar_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(num_2Dvector_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_BCAST(num_3Dvector_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

    call MPI_Barrier(MPI_COMM_WORLD, i_err)

    call allocate_scalar_fields(nxu, nyu, nzu)
    call allocate_vector2D_fields(nxu, nyu)
    call allocate_vector3D_fields(nxu, nyu, nzu)

    call MPI_BCAST(num_filterlengths, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    call MPI_Barrier(MPI_COMM_WORLD, i_err)
    call set_num_filterlengths(num_filterlengths)
    call MPI_BCAST(arr_filterlengths, num_filterlengths, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

    call MPI_FINALIZE(i_err)

end program
