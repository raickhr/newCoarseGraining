module class_student
    implicit none
    type, public:: student
        real, allocatable :: marks(:)
    contains
        procedure :: allocate => allocate_marks
        procedure :: deallocate => deallocate_marks
        procedure :: setmarks => setmarks_arr
    end type student

    contains

    subroutine allocate_marks(self, n)
        class(student), intent(inout) :: self
        integer , intent(in) :: n 
        allocate(self%marks(n))
    end subroutine

    subroutine deallocate_marks(self)
        class(student), intent(inout) :: self
        deallocate(self%marks)
    end subroutine

    subroutine setmarks_arr(self, arr)
        class(student), intent(inout) :: self
        real :: arr(:)
        self%marks = arr
    end subroutine


end module class_student