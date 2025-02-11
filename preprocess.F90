module preprocess
    use kinds
    use constants
    use gridModule
    use mpiwrapper
    use fields
    use input_data_info
    implicit none
    
    contains

    


    subroutine laplaceSmoothLand()
        integer :: num2D_fields, avg_num2D_fields, rem_num2D_fields, &
                   counter, z_count, iti, offset, procid, ii, jj

        integer, allocatable :: start_indices(:),&
                                end_indices(:), &
                                arr_num2D_fields_in_thisproc(:)

        real(kind=real_kind), allocatable :: all2D_fields(:,:,:), &
                                             dummy2D(:,:), &
                                             dummy2D_new(:,:)
        
        num2D_fields = num_scalar2D_fields + (num_scalar3D_fields * nzu)

        if (taskid .EQ. 0) print *,' extending values at land points by laplace smoothing'

        allocate(start_indices(numtasks), &
                 end_indices(numtasks), &
                 arr_num2D_fields_in_thisproc(numtasks))

        allocate(all2D_fields(nxu, nyu, num2D_fields), &
                 dummy2D(nxu, nyu), &
                 dummy2D_new(nxu, nyu))

        offset = 1
        do counter=1, num_scalar2D_fields
            all2D_fields(:,:, offset) = scalar2D_fields(:,:, counter)
            offset = offset+1
        end do
        
        ! if (taskid .EQ. 0) print *,' 2d scalar fields'
        do counter=1, num_scalar3D_fields
            do z_count=1, nzu
                all2D_fields(:,:, offset) = scalar3D_fields(:, :, z_count, counter)
                offset = offset+1
            end do
        end do
        ! print *, 'num 2d fields', offset-1
        ! if (taskid .EQ. 0) print *,' 3d scalar fields'

        avg_num2D_fields = num2D_fields/numtasks
        rem_num2D_fields = mod(num2D_fields, numtasks)

        offset = 1

        do counter =1, numtasks
            if (counter .LE. rem_num2D_fields) then
                arr_num2D_fields_in_thisproc(counter) =  avg_num2D_fields + 1
            else
                arr_num2D_fields_in_thisproc(counter) = avg_num2D_fields
            endif
            
            if (arr_num2D_fields_in_thisproc(counter) > 0 ) then
                start_indices(counter) = offset
                end_indices(counter) = offset + arr_num2D_fields_in_thisproc(counter) -1
                offset = offset + arr_num2D_fields_in_thisproc(counter)
                
            else 
                start_indices(counter) = 0
                end_indices(counter) = 0 
            endif
        end do  
        
        call MPI_Barrier(MPI_COMM_WORLD, i_err)  
        ! Gaussian smoothing 
        if (arr_num2D_fields_in_thisproc(taskid + 1) > 0 ) then
            do counter = start_indices(taskid + 1), end_indices(taskid + 1)
                dummy2D = all2D_fields(:,:, counter)
                do iti =1, niter_laplace_smooth
                    dummy2D_new = 0.0
                    dummy2D_new = dummy2D_new + 0.25* cshift(dummy2D, shift= 1, dim=1) 
                    dummy2D_new = dummy2D_new + 0.25* cshift(dummy2D, shift=-1, dim=1)
                    dummy2D_new = dummy2D_new + 0.25* cshift(dummy2D, shift= 1, dim=2)
                    dummy2D_new = dummy2D_new + 0.25* cshift(dummy2D, shift=-1, dim=2)
                    
                    dummy2D = dummy2D_new
                    dummy2D(1,:) = dummy2D(2,:)
                    dummy2D(nxu,:) = dummy2D(nxu-1,:)
                    dummy2D(:,1) = dummy2D(:,2)
                    dummy2D(:,nyu) = dummy2D(:,nyu-1)
                    
                    where(KMU > 0)
                        dummy2D(:,:) = all2D_fields(:,:, counter)
                    endwhere
                end do
                all2D_fields(:,:, counter) = dummy2D(:,:)
            end do
        endif

        call MPI_Barrier(MPI_COMM_WORLD, i_err)
        if (taskid == 0) print *, 'smooth complete'

        ! Broadcasting smoothed fields
        do procid = 0, numtasks-1
            if (arr_num2D_fields_in_thisproc(procid + 1) > 0 ) then
                do counter = start_indices(procid + 1), end_indices(procid + 1)
                    call MPI_BCAST(all2D_fields(:,:, counter), nxu * nyu, MPI_REAL , procid, MPI_COMM_WORLD, i_err)
                    call MPI_Barrier(MPI_COMM_WORLD, i_err)
                end do
            endif
        end do

        ! reassigning smoothed fields
        offset = 1
        do counter=1, num_scalar2D_fields
            scalar2D_fields(:,:, counter) = all2D_fields(:,:, offset)
            offset = offset+1
        end do
        do counter=1, num_scalar3D_fields
            do z_count=1, nzu
                scalar3D_fields(:, :, z_count, counter) = all2D_fields(:,:, offset)
                offset = offset+1
            end do
        end do

        deallocate(start_indices, end_indices, arr_num2D_fields_in_thisproc)
        deallocate(all2D_fields, dummy2D, dummy2D_new)

        if (taskid .EQ. 0) print *,' extending values at land points by laplace smoothing DONE!'


    end subroutine
end module