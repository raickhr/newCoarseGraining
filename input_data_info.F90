module input_data_info
    use kinds
    use mpiwrapper
    implicit none

    integer (kind=int_kind) :: num_files, num_zlevels, start_timeindex, end_timeindex
    integer (kind=int_kind) , allocatable :: arr_z_index(:)
    real (kind=real_kind) :: timevar_val(1)
    character (len=varname_len) :: timevar_name
    character (len=units_len) :: timevar_units, timevar_calendar
    contains

    subroutine set_numfiles_numzlevels_ntimesinafile(n1, n2, n3, n4)
        integer(kind=int_kind), intent(in) :: n1, n2, n3, n4
        num_files = n1
        num_zlevels = n2
        start_timeindex = n3
        end_timeindex = n4
    end subroutine

    subroutine alloc_arr_z_index()
        if (.not. allocated(arr_z_index)) allocate(arr_z_index(num_zlevels))
    end subroutine

    subroutine set_arr_z_index(list_z_levels)
        integer(kind=int_kind) :: list_z_levels(:)
        arr_z_index = list_z_levels
    end subroutine

    subroutine broadCastInputDataInfo()
        call MPI_BCAST(num_files, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(num_zlevels, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(start_timeindex, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
        call MPI_BCAST(end_timeindex, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

        if (taskid .NE. MASTER) call alloc_arr_z_index()
        call MPI_BCAST(arr_z_index, num_zlevels, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

        call MPI_Barrier(MPI_COMM_WORLD, i_err)
    end subroutine

    subroutine deallocate_inputData_info()
        deallocate(arr_z_index)
    end subroutine

end module
