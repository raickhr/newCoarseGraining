module debug
    implicit none
    logical :: verbose
#ifdef VERBOSE
    verbose = .TRUE.
#else
    verbose = .FALSE.
#endif
end module