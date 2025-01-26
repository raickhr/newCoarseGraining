program main
#include "config.h"
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

    call init_inputDataInfo()

    call allocate_gridVars()

    if (taskid == MASTER) call get_grid_nc(config%InputPath, config%gridFile)
    
    call broadCastGridVars()

    if (taskid == 0) print *, 'GridVars bcasted'

    call init_unfilt_fields()

    if (taskid == 0) print *, 'unfilt vars initialized'

#ifdef HELMHOLTZ
    call init_helmholtz()
#endif

    if (taskid == 0) print *, 'helmholtz initialized'

    call init_filtering()

    if (taskid == 0) print *, 'filtering vars initialized'

    call dividework()

    if (taskid == 0) print *, 'work division initialized'

    call allocate_filtered_fields(nxu, nyu, nzu, num_filterlengths)

    if (taskid == 0) print *, 'filtered vars allocated'

    do file_index = 1, num_files

        if (taskid .EQ. MASTER) then
            print *,"AT FILE NUMBER ",file_index,' OF ', num_files
        endif

        do time_index = start_timeindex, end_timeindex
            if (taskid .EQ. MASTER) then
                print *, "AT TIME INDEX ",time_index-start_timeindex+1,' OF ', end_timeindex - start_timeindex +1
                
                call read_fields(trim(adjustl(config%InputPath))//'/'//trim(adjustl(config%list_filenames(file_index))), &
                                    time_index)
                
                call setAtributesOfHelmHoltzFields()

            endif

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call broadCastReadFields()

            if (taskid .EQ. MASTER) print *, ' Read fields bcasted'

#ifdef HELMHOLTZ
            call helmholtzDecompAllVecFields()
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            if (taskid == 0 ) print *, 'completed Helmholtz decomposition'

            if (taskid == 0 ) then
                WRITE(writefilename, "(A13,I0.3,A3)") "_phi_psi_time", time_index, ".nc"
                call makeFileName(config%list_filenames(file_index), writefilename)
                call writeUnfiltHelmHoltzDeompFields( trim(adjustl(config%OutputPath))//'/'//writefilename, &
                &   'xi_rho', 'eta_rho', trim(adjustl(config%vertdim_name)), trim(adjustl(config%timevar_name)))
            end if
#else
            if (taskid == 0 ) then
                WRITE(writefilename, "(A13,I0.3,A3)") "_phi_psi_time", time_index, ".nc"
                call makeFileName(config%list_filenames(file_index), writefilename)
                call readHelmHoltzDeompFields(trim(adjustl(config%OutputPath))//'/'//writefilename, time_index)
            end if
            
#endif

            call broadCastPsiPhiFields()

            if (taskid .EQ. MASTER) print *, 'phi psi broadcasted'
            
            call filter_allvars()
            
            if (taskid .EQ. MASTER) then 
                print *, ''
                print *, "TIME INDEX ",time_index-start_timeindex+1,' OF ', end_timeindex - start_timeindex +1, ' COMPLETED'
                print *, "AT FILE NUMBER ",file_index,' OF ', num_files
                print *, ''
            end if

            !call collectFilteredFields(arr_numcols_inallprocs, arr_startcolindex_inallprocs, num_filterlengths)
            call collectCoarseGrainedFields(arr_numcols_inallprocs, arr_startcolindex_inallprocs, num_filterlengths)

            if (taskid .EQ. MASTER) then 
                print *, ''
                print *, 'all filtered variables collected !'
                print *, ''
                print *, 'Now calculating the vectos from the Helmholtz Potentials !'
                call getVectorsFromFilteredHelmHoltzPotentials(num_filterlengths)

            end if
            
            if (taskid .EQ. MASTER) then
                WRITE(writefilename, "(A10,I0.3,A3)") "_filt_time", time_index, ".nc"
                call makeFileName(config%list_filenames(file_index), writefilename)

                ! call writeDirectFilteredFields( trim(adjustl(config%OutputPath))//'/'//writefilename, &
                ! &                 'xi_rho', 'eta_rho', trim(adjustl(config%vertdim_name)), &
                ! &                 'Lengthscale', trim(adjustl(config%timevar_name)))

                call writeCoarseGrainedFields( trim(adjustl(config%OutputPath))//'/'//writefilename, &
                &                 'xi_rho', 'eta_rho', trim(adjustl(config%vertdim_name)), trim(adjustl(config%wvertdim_name)), &
                &                 'Lengthscale', trim(adjustl(config%timevar_name)))
            end if

        end do !close time loop
    end do ! close

#ifdef HELMHOLTZ
    call delMultiGrid()
#endif
    
    call deallocate_gridVars()
    call deallocate_filtervars()
    call deallocate_inputData_info()
    call delloacate_fields()
    call config%destruct()

    call MPI_Barrier(MPI_COMM_WORLD, i_err)
    if (taskid .EQ. MASTER) print *, 'Ending program'
    call stopMPI()

end program
