module input_data_info
    use kinds
    use mpiwrapper
    use configurationMod
    implicit none

    integer (kind=int_kind) :: num_files, num_zlevels, start_timeindex, end_timeindex
    integer (kind=int_kind) , allocatable :: arr_z_index(:)
    real (kind=real_kind) :: timevar_val(1)
    real (kind=real_kind), allocatable :: vertDim_vals(:)
    character (len=varname_len) :: timevar_name
    character (len=varname_len) :: vertdim_name
    character (len=units_len) :: timevar_units, timevar_calendar
    contains

    subroutine set_numfiles_numzlevels_ntimesinafile(n1, n2, n3, n4, tvar_name, vdim_name)
        integer(kind=int_kind), intent(in) :: n1, n2, n3, n4
        character(len=varname_len), intent(in) :: tvar_name, vdim_name
        num_files = n1
        num_zlevels = n2
        start_timeindex = n3
        end_timeindex = n4
        timevar_name = trim(adjustl(tvar_name))
        vertdim_name = trim(adjustl(vdim_name))
    end subroutine

    subroutine set_arr_z_index(list_z_levels)
        integer(kind=int_kind) :: list_z_levels(:)
        arr_z_index = list_z_levels
    end subroutine


    subroutine deallocate_inputData_info()
        deallocate(arr_z_index)
    end subroutine

    subroutine init_inputDataInfo()
        if (taskid == MASTER) call set_numfiles_numzlevels_ntimesinafile(config%num_of_files_to_read, config%nz, &
                                                                         config%startTimeIndex, config%endTimeIndex, &
                                                                         config%timevar_name, config%vertdim_name)

        call MPI_BCAST(num_files, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(num_zlevels, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(start_timeindex, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(end_timeindex, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

        if (.not. allocated(arr_z_index)) allocate(arr_z_index(num_zlevels))
        if (.not. allocated(vertDim_vals)) allocate(vertDim_vals(num_zlevels))

        if (taskid == MASTER) then
            arr_z_index = config%list_zlevels
        endif

        call MPI_BCAST(arr_z_index, num_zlevels, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

        call MPI_Barrier(MPI_COMM_WORLD, i_err)
    end subroutine

end module
