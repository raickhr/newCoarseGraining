! This module defines the data type for common data types
module kinds
    implicit none
    integer, parameter :: char_len = 150, &
                        & filename_len = 100, &
                        & pathname_len = 500, &
                        & fullfilename_len = 600, &
                        & varname_len = 50, &
                        & units_len = 100, &
                        & longname_len = 200, &
                        & dimname_len = 20, &
                        & int_kind = kind(5), &
                        & log_kind = kind(.true.), &
                        & real_kind = selected_real_kind(8), &
                        & dbl_kind = selected_real_kind(13)
end module
