module filterparallel
    use kinds
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

end module
