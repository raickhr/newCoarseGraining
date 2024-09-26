module filterparallel
    use kinds
    use mpiwrapper
    implicit none

    integer(kind=int_kind) :: num_filterlengths
    real(kind=real_kind), allocatable, dimension(:) :: arr_filterlengths
    contains

        subroutine set_num_filterlengths(n)
            integer(kind=int_kind), intent(in) :: n
            num_filterlengths = n
            if (.not. allocated(arr_filterlengths)) allocate(arr_filterlengths(n))
        end subroutine

        subroutine set_arr_filterlength(list_filterlength)
            real(kind=real_kind), intent(in) :: list_filterlength(:)
            arr_filterlengths = list_filterlength
        end subroutine

        subroutine broadCastFilterInfo()
            call MPI_BCAST(num_filterlengths, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
            if (taskid .NE. MASTER) call set_num_filterlengths(num_filterlengths)
            call MPI_BCAST(arr_filterlengths, num_filterlengths, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
        end subroutine

        subroutine deallocate_filtervars()
            deallocate(arr_filterlengths)
        end subroutine

end module
