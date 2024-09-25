! This module defines the data type for common data types
module kinds
    implicit none
    integer, parameter :: char_len = 150, &
                        & int_kind = kind(5), &
                        & log_kind = kind(.true.), &
                        & real_kind = selected_real_kind(6), &
                        & dbl_kind = selected_real_kind(13)
end module
