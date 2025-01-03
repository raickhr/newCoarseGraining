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
    
    integer :: file_index, time_index, z_index
    character (len=filename_len) :: writefilename

    call startMPI()

    call init_config()

    call init_grid()

    call init_unfilt_fields()

    call init_helmholtz()

    call init_filtering()

    call setMultiGrid()

    call dividework()

    call allocate_filtered_fields()

    do file_index = 1, config%num_of_files_to_read

        if (taskid .EQ. MASTER) then
            print *,"AT FILE NUMBER ",file_index,' OF ', config%num_of_files_to_read
        endif

        do time_index = config%startTimeIndex, config%endTimeIndex
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            if (taskid .EQ. MASTER) then
                print *,"AT TIME INDEX ",time_index-config%startTimeIndex+1,' OF ', config%endTimeIndex - config%startTimeIndex +1

                call read_fields(trim(adjustl(config%InputPath))//'/'//trim(adjustl(config%list_filenames(file_index))), &
                                    time_index)

            endif

            call broadCastReadFields()

            call helmholtzDecompAllVecFields()

            if (taskid .EQ. MASTER) then
                WRITE(writefilename, "(A5,I0.3,A5,I0.3,A3)") "phi_psi_file", file_index, "_time", time_index, ".nc"
                call writeHelmHoltzDeompFields( trim(adjustl(config%OutputPath))//'/'//writefilename, &
                &   'xi_rho', 'eta_rho', trim(adjustl(config%vertdim_name)), trim(adjustl(config%timevar_name)))
            end if

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            if (taskid == 0 ) print *, 'completed Helmholtz decomposition'

            call broadCastPsiPhiFields()

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

            !call collectFilteredFields(arr_numcols_inallprocs, arr_startcolindex_inallprocs, num_filterlengths)

            call collectCoarseGrainedFields(arr_numcols_inallprocs, arr_startcolindex_inallprocs, num_filterlengths)

            if (taskid .EQ. MASTER) then 
                print *, ''
                print *, 'all filtered variables collected !'
                print *, ''
            end if
            
            if (taskid .EQ. MASTER) then
                WRITE(writefilename, "(A5,I0.3,A5,I0.3,A3)") "file", file_index, "_time", time_index, ".nc"
                call writeHelmHoltzDeompFields(fullfilename, x_dimname, y_dimname, z_dimname, time_dimname) 
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
