module filterparallel
    use kinds
    use mpiwrapper
    use gridModule
    use fields
    implicit none

    integer(kind=int_kind) :: num_filterlengths
    real(kind=real_kind), allocatable, dimension(:) :: arr_filterlengths

    integer(kind=int_kind), allocatable, dimension(:) :: arr_numcols_inallprocs, arr_startcolindex_inallprocs
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

        subroutine dividework()
            integer :: avg_numcols, reminder, offset, counter, ncols

            allocate(arr_numcols_inallprocs(numtasks), &
                    arr_startcolindex_inallprocs(numtasks))

            avg_numcols = nyu/numtasks
            reminder = mod(nyu, numtasks)

            offset = 1

            do counter =1, numtasks
                if (counter .LE. reminder) then
                    ncols =  avg_numcols + 1
                else
                    ncols = avg_numcols
                endif
                arr_numcols_inallprocs(counter) = ncols
                arr_startcolindex_inallprocs(counter) = offset
                ! if (taskid .EQ. MASTER) then
                !     WRITE(*,'(A7, I4, A10, I4, A15, I4, A15, I4)') 'taskid',  counter-1, 'ncols', ncols, 'start', offset , 'end', offset + ncols -1
                ! endif
                offset = offset + ncols
                
            end do  
            
            call MPI_Barrier(MPI_COMM_WORLD, i_err)         
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


        subroutine filter_allvars()
            integer :: startJindex, endJindex, &
                    &  startIindex, endIindex, &
                    &  i_index, j_index

            startJindex = arr_startcolindex_inallprocs(taskid+1)
            endJindex = startJindex + arr_numcols_inallprocs(taskid+1) -1

            startIindex = 1
            endIindex = nxu

            print *, 'at proc', taskid, startJindex, endJindex

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            do j_index=startJindex, endJindex
                do i_index = startIindex, endIindex
                    ! condition for skipping the filterpoint for eg land
                    call filter_allVarsAtpoint(i_index, j_index)
                end do
            end do
        end subroutine

        subroutine filter_allVarsAtpoint(i_index, j_index)
            integer , intent(in) :: i_index, j_index
            ! local variables

            integer :: filter_counter, depth_counter, varcounter, numvars, counter, &
                   &   east_west_BoxSize, north_south_BoxSize, & 
                   &   west_cornerindex, east_cornerindex, &
                   &   south_cornerindex, north_cornerindex

            real(kind=real_kind) :: filterlengthInKM

            real(kind=real_kind), allocatable, dimension(:,:,:) :: unfiltvars ! x, y, nvars
            real(kind=real_kind), allocatable, dimension(:) :: filtered_vars ! nvars
            real(kind=real_kind), allocatable, dimension(:,:) :: kernelVal ! x, y


            numvars = num_scalar_fields * nzu + &
                      2 * num_2Dvector_fields + &
                      3 * nzu * num_3Dvector_fields   
                      
            allocate(filtered_vars(numvars))

            do filter_counter=1, num_filterlengths 

                filterlengthInKM = arr_filterlengths(filter_counter)

                call get_boxcorners_lat_lon_grid(i_index, j_index, filterlengthInKM, &
                                    &       west_cornerindex, east_cornerindex, &
                                    &       south_cornerindex, north_cornerindex, &
                                    &       east_west_BoxSize, north_south_BoxSize )

                allocate(unfiltvars(east_west_BoxSize, north_south_BoxSize, numvars))
                allocate(kernelVal(east_west_BoxSize, north_south_BoxSize) ) 
                
                call get_kernel(i_index, j_index, &       
                &       west_cornerindex, east_cornerindex, &
                &       south_cornerindex, north_cornerindex, &
                &       east_west_BoxSize, north_south_BoxSize, &
                &       filterlengthInKM, kernelVal)

                call groupUnfiltvars(unfiltvars, east_west_BoxSize, north_south_BoxSize, numvars, &
                &                    west_cornerindex, east_cornerindex, &
                &                    south_cornerindex, north_cornerindex )
                
                call get_filteredVals_allVars(east_west_BoxSize, north_south_BoxSize, numvars, unfiltvars, kernelVal, filtered_vars(:) )
                call assignFilteredVars(numvars,i_index, j_index, filter_counter, filtered_vars)

                deallocate(unfiltvars)
                deallocate(kernelVal) 

            end do
            deallocate(filtered_vars)

        end subroutine

        subroutine groupUnfiltvars(unfiltvars, nx, ny, nvars, &
            &                      west_cornerindex, east_cornerindex, &
            &                      south_cornerindex, north_cornerindex )

            integer(kind=int_kind), intent(in) :: nx, ny, nvars, &
            &                                     west_cornerindex, east_cornerindex, &
            &                                     south_cornerindex, north_cornerindex

            real(kind=real_kind), intent(out) :: unfiltvars(nx, ny, nvars)

            integer:: counter, varcounter, depth_counter

            varcounter = 1
            do counter=1, num_scalar_fields
                do depth_counter = 1, nzu
                    unfiltvars(:,:, varcounter) = &
                        &   scalar_fields(west_cornerindex:east_cornerindex,south_cornerindex:north_cornerindex, depth_counter, counter)
                    varcounter = varcounter + 1
                end do
            enddo

            do counter=1, num_2Dvector_fields
                unfiltvars(:,:, varcounter) = &
                    &   vector2DX_fields(west_cornerindex:east_cornerindex,south_cornerindex:north_cornerindex, counter)
                varcounter = varcounter + 1
                unfiltvars(:,:, varcounter) = &
                    &   vector2DY_fields(west_cornerindex:east_cornerindex,south_cornerindex:north_cornerindex, counter)
                varcounter = varcounter + 1
            end do

            do counter=1, num_3Dvector_fields
                do depth_counter = 1, nzu
                    unfiltvars(:,:, varcounter) = &
                        &   vector3DX_fields(west_cornerindex:east_cornerindex,south_cornerindex:north_cornerindex,depth_counter, counter)
                    varcounter = varcounter + 1

                    unfiltvars(:,:, varcounter) = &
                        &   vector3DY_fields(west_cornerindex:east_cornerindex,south_cornerindex:north_cornerindex,depth_counter, counter)
                    varcounter = varcounter + 1

                    unfiltvars(:,:, varcounter) = &
                        &   vector3DZ_fields(west_cornerindex:east_cornerindex,south_cornerindex:north_cornerindex,depth_counter, counter)
                    varcounter = varcounter + 1
                end do
            end do
        end subroutine

        subroutine assignFilteredVars(nvars, i_index, j_index, filter_index, filteredFields)
            integer, intent(in) :: nvars, i_index, j_index, filter_index
            real(kind=real_kind) :: filteredFields(nvars)

            integer(kind=int_kind) :: varcounter, depth_index, counter

            varcounter = 1
            do counter=1, num_scalar_fields
                do depth_index = 1, nzu
                    ! WRITE(*,'(A14, I4, A14, I4, A14, I4, A14, I4, A14, I4, A14, I4, A14, I4)')  'taskid', taskid, &
                    !                                                 'i_index', i_index, & 
                    !                                                 'j_index', j_index, & 
                    !                                                 'depth_index', depth_index, & 
                    !                                                 'counter', counter, & 
                    !                                                 'filter_index', filter_index,  &
                    !                                                 'varcounter', varcounter 

                    OL_scalar_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                end do
            enddo

            do counter=1, num_2Dvector_fields
                OL_vector2DX_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
                OL_vector2DY_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
            end do

            do counter=1, num_3Dvector_fields
                do depth_index = 1, nzu
                    OL_vector3DX_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                    OL_vector3DY_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                    OL_vector3DZ_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                end do
            end do
        end subroutine

        subroutine get_kernel(i_index, j_index, &
            &       west_cornerindex, east_cornerindex, &
            &       south_cornerindex, north_cornerindex, &
            &       east_west_BoxSize, north_south_BoxSize, &
            &       filterlengthInKM, kernelVal)

            integer(kind=int_kind), intent(in) :: i_index, j_index, &
                                        &         west_cornerindex, east_cornerindex, &
                                        &         south_cornerindex, north_cornerindex, &
                                        &         east_west_BoxSize, north_south_BoxSize

            real(kind=real_kind), intent(in) :: filterlengthInKM

            real(kind=real_kind), intent(out) :: kernelVal(east_west_BoxSize, north_south_BoxSize)

            real(kind=real_kind) :: arr_UAREA(east_west_BoxSize, north_south_BoxSize), &
                                &   arr_ULAT(east_west_BoxSize, north_south_BoxSize), &
                                &   arr_ULONG(east_west_BoxSize, north_south_BoxSize), &
                                &   great_circ_dist(east_west_BoxSize, north_south_BoxSize)

            real(kind=real_kind) :: center_long, center_lat, pointfive(east_west_BoxSize, north_south_BoxSize), &
                                    Ell_filter_by2(east_west_BoxSize, north_south_BoxSize)

            center_long = ULONG(i_index, j_index)
            center_lat = ULAT(i_index, j_index)

            arr_UAREA(:,:) = UAREA(west_cornerindex: east_cornerindex, south_cornerindex:north_cornerindex)
            arr_ULONG(:,:) = ULONG(west_cornerindex: east_cornerindex, south_cornerindex:north_cornerindex)
            arr_ULAT(:,:) = ULAT(west_cornerindex: east_cornerindex, south_cornerindex:north_cornerindex) 

            call calc_greatCircDist(center_long, center_lat, east_west_BoxSize, north_south_BoxSize, arr_ULONG, arr_ULAT, great_circ_dist)

            Ell_Filter_by2(:,:) =  filterlengthInKM/2 * 1d3 ! turn into meters
            pointfive(:,:) = 0.5

            kernelVal = pointfive-0.5*tanh((great_circ_dist-Ell_Filter_by2)/10.0) * arr_UAREA
            kernelVal = kernelVal/sum(kernelVal)
        
        end subroutine

        subroutine calc_greatCircDist(center_long, center_lat, xshape, yshape, arr_ULONG, arr_ULAT, distance)
            real(kind=real_kind), intent(in):: center_lat, center_long
            integer(kind=int_kind), intent(in) :: xshape, yshape
            real(kind=real_kind), intent(in) :: arr_ULONG( xshape, yshape), arr_ULAT( xshape, yshape)
            real(kind=real_kind), intent(out):: distance( xshape, yshape)
            
            real(kind=real_kind):: dlambda( xshape, yshape), phi1( xshape, yshape), phi2( xshape, yshape), &
                          &        dsigma( xshape, yshape), numerator( xshape, yshape), denominator( xshape, yshape)


            dlambda = arr_ULONG - center_long
            phi1 = center_lat
            phi2 = arr_ULAT

            numerator = ( cos(phi2)*sin(dlambda) )**2 + (cos(phi1)*sin(phi2) -sin(phi1)*cos(phi2)*cos(dlambda))**2
            numerator = sqrt(numerator)
            denominator = sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(dlambda)

            dsigma = atan2(numerator, denominator)
            distance = radius * dsigma

        end subroutine

        subroutine get_boxcorners_lat_lon_grid(i_index, j_index, filterlengthInKM, &
                                       &       west_cornerindex, east_cornerindex, &
                                       &       south_cornerindex, north_cornerindex, &
                                       &       east_west_BoxSize, north_south_BoxSize )

            integer(kind=int_kind), intent(in) :: i_index, j_index
            real(kind=real_kind), intent(in) :: filterlengthInKM
            integer(kind=int_kind), intent(out) ::west_cornerindex, east_cornerindex, &
                                        &       south_cornerindex, north_cornerindex, &
                                        &       east_west_BoxSize, north_south_BoxSize
                                

            integer :: half_east_west_BoxSize, half_north_south_BoxSize  

            half_east_west_BoxSize = (0.6 * filterlengthInKM* 1d3 / DXU(nxu/2, j_index))+5   !.6 and 5 for tolerance
            half_north_south_BoxSize = (0.6 * filterlengthInKM* 1d3 / DYU(nxu/2, j_index))+5   !.6 and 5 for tolerance

            west_cornerindex = i_index - half_east_west_BoxSize
            east_cornerindex = i_index + half_east_west_BoxSize
            south_cornerindex = j_index - half_north_south_BoxSize
            north_cornerindex = j_index - half_north_south_BoxSize

            if (west_cornerindex < 1) west_cornerindex = 1
            if (east_cornerindex > nxu) east_cornerindex = nxu
            if (south_cornerindex < 1) south_cornerindex = 1
            if (north_cornerindex > nyu) north_cornerindex = nyu

            east_west_BoxSize = east_cornerindex - west_cornerindex + 1
            north_south_BoxSize = north_cornerindex - south_cornerindex + 1

        end subroutine

        subroutine get_filteredVals_allVars(nx, ny, nvars, allVars, arr_kernel, filtered_vars )
            integer(kind=int_kind) :: nx, ny, nvars
            real(kind=real_kind) :: allVars(nx,ny,nvars), arr_kernel(nx,ny)
            real(kind=real_kind), intent(out):: filtered_vars(nvars)
            integer :: counter 

            do counter = 1, nvars
                filtered_vars(counter) = sum(allVars(:,:, counter) * arr_kernel(:,:))
            end do
  
        end subroutine

end module
