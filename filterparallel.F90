module filterparallel
    use kinds
    use mpiwrapper
    use gridModule
    use fields
    use configurationMod
    use netcdf_io
    implicit none

    integer(kind=int_kind) :: num_filterlengths
    real(kind=real_kind), allocatable, dimension(:) :: arr_filterlengths

    integer(kind=int_kind), allocatable, dimension(:) :: arr_numcols_inallprocs, arr_startcolindex_inallprocs
    contains

        subroutine init_filtering()
            if (taskid == MASTER) call set_num_filterlengths(config%nfilter)

            call MPI_BCAST(num_filterlengths, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            allocate(arr_filterlengths(num_filterlengths))

            if (taskid == MASTER) call set_arr_filterlength(config%list_filterLength)
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call MPI_BCAST(arr_filterlengths, num_filterlengths, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
        end subroutine

        subroutine set_num_filterlengths(n)
            integer(kind=int_kind), intent(in) :: n
            num_filterlengths = n
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
                if (taskid .EQ. MASTER) then
                    WRITE(*,'(A7, I4, A10, I4, A15, I4, A15, I4)') 'taskid',  &
                    counter-1, 'ncols', ncols, 'start', offset , 'end', offset + ncols -1
                endif
                offset = offset + ncols
                
            end do  
            
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

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            do j_index=startJindex, endJindex
                WRITE(*, '(A12, I4, A10, I4, A5, I4)') 'At taskid: ', taskid, '  column ', & 
                j_index - startJindex + 1, ' of ', endJindex - startJindex +1
                do i_index = startIindex, endIindex
                    ! condition for skipping the filterpoint for eg land
                    !call direct_filter_allVarsAtpoint(i_index, j_index)
                    !print *,'i_index, j_index, filterlengthInKM', i_index, j_index, filterlengthInKM
                    call coarseGrain_FieldsAtpoint(i_index, j_index)
                    !call MPI_Barrier(MPI_COMM_WORLD, i_err)
                    
                end do
                !print *, 'at taskid', taskid, 'filtered at all filterlengths at j_index', j_index

            end do
            
            
            
        end subroutine

        subroutine coarseGrain_FieldsAtpoint(i_index, j_index)
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


            numvars = num_scalar2D_fields + &
                      nzu * num_scalar3D_fields + &
                      2 * num_vector2D_fields   + & !psi phi ! do not filter vec_2D_phi and vec_2D_psi
                      2 * nzu * num_vector3D_fields  + &! psi phi for 2d fields only vertical is directly filtered
                      1 * (nzu+1) * num_vector3D_fields  !There is one additional layer for vertical vel
                     
                      
            

            do filter_counter=1, num_filterlengths 

                filterlengthInKM = arr_filterlengths(filter_counter)
                
                call get_boxcorners_lat_lon_grid(i_index, j_index, filterlengthInKM, &
                                    &       west_cornerindex, east_cornerindex, &
                                    &       south_cornerindex, north_cornerindex, &
                                    &       east_west_BoxSize, north_south_BoxSize )

                allocate(unfiltvars(east_west_BoxSize, north_south_BoxSize, numvars))
                allocate(filtered_vars(numvars))
                allocate(kernelVal(east_west_BoxSize, north_south_BoxSize) )
                
                call get_kernel(i_index, j_index, &       
                &       west_cornerindex, east_cornerindex, &
                &       south_cornerindex, north_cornerindex, &
                &       east_west_BoxSize, north_south_BoxSize, &
                &       filterlengthInKM, kernelVal)

                call groupUnfiltVarsForCoarseGraining(unfiltvars, east_west_BoxSize, north_south_BoxSize, numvars, &
                &                    west_cornerindex, east_cornerindex, &
                &                    south_cornerindex, north_cornerindex )
                
                call get_filteredVals_allVars(east_west_BoxSize, north_south_BoxSize, numvars, &
                                              unfiltvars, kernelVal, filtered_vars(:) )

                call assignCoarseGrainedVars(numvars,i_index, j_index, filter_counter, filtered_vars)

                deallocate(unfiltvars)
                deallocate(filtered_vars)
                deallocate(kernelVal) 

            end do
        end subroutine

        subroutine direct_filter_allVarsAtpoint(i_index, j_index)
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


            numvars = num_scalar3D_fields + &
                      nzu * num_scalar3D_fields + &
                      2 * num_vector2D_fields + &
                      3 * nzu * num_vector3D_fields   
                      
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

                call groupUnfiltvarsForDirectFiltering(unfiltvars, east_west_BoxSize, north_south_BoxSize, numvars, &
                &                    west_cornerindex, east_cornerindex, &
                &                    south_cornerindex, north_cornerindex )

                
                call get_filteredVals_allVars(east_west_BoxSize, north_south_BoxSize, numvars, &
                                              unfiltvars, kernelVal, filtered_vars(:) )
                call assignFilteredVars(numvars,i_index, j_index, filter_counter, filtered_vars)

                deallocate(unfiltvars)
                deallocate(kernelVal) 

            end do
            deallocate(filtered_vars)

        end subroutine

        subroutine groupUnfiltvarsForDirectFiltering(unfiltvars, nx, ny, nvars, &
            &                      west_cornerindex, east_cornerindex, &
            &                      south_cornerindex, north_cornerindex )

            integer(kind=int_kind), intent(in) :: nx, ny, nvars, &
            &                                     west_cornerindex, east_cornerindex, &
            &                                     south_cornerindex, north_cornerindex

            real(kind=real_kind), intent(out) :: unfiltvars(nx, ny, nvars)

            integer:: counter, varcounter, depth_counter

            varcounter = 1
            do counter=1, num_scalar2D_fields
                unfiltvars(:,:, varcounter) = &
                    &   scalar2D_fields(west_cornerindex:east_cornerindex, south_cornerindex:north_cornerindex, counter)
                varcounter = varcounter + 1
            enddo

            do counter=1, num_scalar3D_fields
                do depth_counter = 1, nzu
                    unfiltvars(:,:, varcounter) = &
                        &   scalar3D_fields(west_cornerindex:east_cornerindex, south_cornerindex:north_cornerindex, &
                                            depth_counter, counter)
                    varcounter = varcounter + 1
                end do
            enddo

            do counter=1, num_vector2D_fields
                unfiltvars(:,:, varcounter) = &
                    &   vector2DX_fields(west_cornerindex:east_cornerindex, south_cornerindex:north_cornerindex, counter)
                varcounter = varcounter + 1

                unfiltvars(:,:, varcounter) = &
                    &   vector2DY_fields(west_cornerindex:east_cornerindex, south_cornerindex:north_cornerindex, counter)
                varcounter = varcounter + 1
            end do

            do counter=1, num_vector3D_fields
                do depth_counter = 1, nzu
                    unfiltvars(:,:, varcounter) = &
                        &   vector3DX_fields(west_cornerindex:east_cornerindex, south_cornerindex:north_cornerindex, &
                                             depth_counter, counter)
                    varcounter = varcounter + 1

                    unfiltvars(:,:, varcounter) = &
                        &   vector3DY_fields(west_cornerindex:east_cornerindex, south_cornerindex:north_cornerindex, &
                                             depth_counter, counter)
                    varcounter = varcounter + 1

                    unfiltvars(:,:, varcounter) = &
                        &   vector3DZ_fields(west_cornerindex:east_cornerindex, south_cornerindex:north_cornerindex, &
                                             depth_counter, counter)
                    varcounter = varcounter + 1
                end do
            end do
        end subroutine


        subroutine groupUnfiltVarsForCoarseGraining(unfiltvars, nx, ny, nvars, &
            &                      west_cornerindex, east_cornerindex, &
            &                      south_cornerindex, north_cornerindex )

            integer(kind=int_kind), intent(in) :: nx, ny, nvars, &
            &                                     west_cornerindex, east_cornerindex, &
            &                                     south_cornerindex, north_cornerindex

            real(kind=real_kind), intent(out) :: unfiltvars(nx, ny, nvars)

            integer:: counter, varcounter, depth_counter

            varcounter = 1

            do counter=1, num_scalar2D_fields
                unfiltvars(:,:, varcounter) = &
                    &   scalar2D_fields(west_cornerindex:east_cornerindex,   &
                                        south_cornerindex:north_cornerindex, &
                                        counter)
                varcounter = varcounter + 1
            enddo

            do counter=1, num_scalar3D_fields
                do depth_counter = 1, nzu
                    unfiltvars(:,:, varcounter) = &
                        &   scalar3D_fields(west_cornerindex:east_cornerindex,   &
                                            south_cornerindex:north_cornerindex, &
                                            depth_counter, counter)
                    varcounter = varcounter + 1
                end do
            enddo

            do counter=1, num_vector2D_fields
                unfiltvars(:,:, varcounter) = &
                    &   phi2D_fields(west_cornerindex:east_cornerindex,   &
                                     south_cornerindex:north_cornerindex, &
                                     counter)
                varcounter = varcounter + 1

                unfiltvars(:,:, varcounter) = &
                    &   psi2D_fields(west_cornerindex:east_cornerindex,   &
                                     south_cornerindex:north_cornerindex, &
                                     counter)
                varcounter = varcounter + 1
            end do

            do counter=1, num_vector3D_fields
                do depth_counter = 1, nzu
                    unfiltvars(:,:, varcounter) = &
                        &   phi3D_fields(west_cornerindex:east_cornerindex,  &
                                         south_cornerindex:north_cornerindex, &
                                        depth_counter, counter)
                    varcounter = varcounter + 1

                    unfiltvars(:,:, varcounter) = &
                        &   psi3D_fields(west_cornerindex:east_cornerindex,   &
                                         south_cornerindex:north_cornerindex, &
                                         depth_counter, counter)
                    varcounter = varcounter + 1

                    unfiltvars(:,:, varcounter) = &
                        &   vector3DZ_fields(west_cornerindex:east_cornerindex, &
                                             south_cornerindex:north_cornerindex, &
                                             depth_counter, counter)
                    varcounter = varcounter + 1
                end do
                ! additional layer because they are in face
                unfiltvars(:,:, varcounter) = &
                        &   vector3DZ_fields(west_cornerindex:east_cornerindex, &
                                             south_cornerindex:north_cornerindex, &
                                             nzu+1, counter)
                varcounter = varcounter + 1
            end do
        end subroutine

        subroutine assignFilteredVars(nvars, i_index, j_index, filter_index, filteredFields)
            integer, intent(in) :: nvars, i_index, j_index, filter_index
            real(kind=real_kind) :: filteredFields(nvars)

            integer(kind=int_kind) :: varcounter, depth_index, counter

            varcounter = 1
            do counter=1, num_scalar2D_fields
                OL_scalar2D_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
            enddo

            do counter=1, num_scalar3D_fields
                do depth_index = 1, nzu
                    OL_scalar3D_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                end do
            enddo


            do counter=1, num_vector2D_fields
                OL_vector2DX_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
                OL_vector2DY_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
            end do

            do counter=1, num_vector3D_fields
                do depth_index = 1, nzu
                    OL_vector3DX_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                    OL_vector3DY_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                    OL_vector3DZ_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                end do
                ! additional layer because they are in face
                OL_vector3DZ_fields(i_index, j_index, nzu+1, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1  
            end do
        end subroutine

        subroutine assignFilteredHelmHoltzVars(nvars, i_index, j_index, filter_index, filteredFields)
            integer, intent(in) :: nvars, i_index, j_index, filter_index
            real(kind=real_kind) :: filteredFields(nvars)

            integer(kind=int_kind) :: varcounter, depth_index, counter

            varcounter = 1
            do counter=1, num_vector2D_fields
                OL_phi2D_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
                OL_psi2D_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
            end do

            do counter=1, num_vector3D_fields
                do depth_index = 1, nzu
                    OL_phi3D_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                    OL_psi3D_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                end do
            end do
        end subroutine

        subroutine assignCoarseGrainedVars(nvars, i_index, j_index, filter_index, filteredFields)
            integer, intent(in) :: nvars, i_index, j_index, filter_index
            real(kind=real_kind) :: filteredFields(nvars)

            integer(kind=int_kind) :: varcounter, depth_index, counter

            varcounter = 1

            do counter=1, num_scalar2D_fields
                OL_scalar2D_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
            enddo

            do counter=1, num_scalar3D_fields
                do depth_index = 1, nzu

                    OL_scalar3D_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1
                end do
            enddo

            do counter=1, num_vector2D_fields

                OL_phi2D_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1

                OL_psi2D_fields(i_index, j_index, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
            end do

            do counter=1, num_vector3D_fields
                do depth_index = 1, nzu

                    OL_phi3D_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1

                    OL_psi3D_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1

                    OL_vector3DZ_fields(i_index, j_index, depth_index, counter, filter_index) = filteredFields(varcounter)
                    varcounter = varcounter + 1

                end do
                ! additional layer because they are in face
                OL_vector3DZ_fields(i_index, j_index, nzu+1, counter, filter_index) = filteredFields(varcounter)
                varcounter = varcounter + 1
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

            real(kind=real_kind), allocatable, dimension(:,:) :: &
                                &   arr_UAREA, arr_ULAT, &
                                &   arr_ULONG, great_circ_dist, &
                                &   pointfive, Ell_filter_by2

            real(kind=real_kind) :: center_long, center_lat 
                                    

            center_long = ULONG(i_index, j_index)
            center_lat = ULAT(i_index, j_index)

            allocate(arr_UAREA(east_west_BoxSize, north_south_BoxSize), &
                &   arr_ULAT(east_west_BoxSize, north_south_BoxSize), &
                &   arr_ULONG(east_west_BoxSize, north_south_BoxSize), &
                &   great_circ_dist(east_west_BoxSize, north_south_BoxSize), &
                &   pointfive(east_west_BoxSize, north_south_BoxSize), &
                &   Ell_filter_by2(east_west_BoxSize, north_south_BoxSize))


            arr_UAREA(:,:) = UAREA(west_cornerindex: east_cornerindex, south_cornerindex:north_cornerindex)
            arr_ULONG(:,:) = ULONG(west_cornerindex: east_cornerindex, south_cornerindex:north_cornerindex)
            arr_ULAT(:,:) = ULAT(west_cornerindex: east_cornerindex, south_cornerindex:north_cornerindex) 

            
            call calc_greatCircDist(center_long, center_lat, east_west_BoxSize, north_south_BoxSize, &
                                arr_ULONG, arr_ULAT, great_circ_dist)

            great_circ_dist = great_circ_dist * 1d-3 ! turn into kilometers
            Ell_Filter_by2(:,:) =  filterlengthInKM/2 
            pointfive(:,:) = 0.5

            kernelVal = (pointfive-0.5*tanh((great_circ_dist-Ell_Filter_by2)/10.0)) * arr_UAREA
            where (great_circ_dist > (1.1 *Ell_filter_by2)) ! tolerance 10%
                kernelVal = 0
            end where
            kernelVal = kernelVal/sum(kernelVal)
            !if (taskid .EQ. 0) print *, sum(kernelVal)
            deallocate(arr_UAREA, arr_ULAT, &
            &   arr_ULONG, great_circ_dist, &
            &   pointfive, Ell_filter_by2)
        
        end subroutine

        subroutine calc_greatCircDist(center_long, center_lat, xshape, yshape, arr_ULONG, arr_ULAT, distance, debugflag)
            real(kind=real_kind), intent(in):: center_lat, center_long
            integer(kind=int_kind), intent(in) :: xshape, yshape
            real(kind=real_kind), intent(in) :: arr_ULONG( xshape, yshape), arr_ULAT( xshape, yshape)
            real(kind=real_kind), intent(out):: distance( xshape, yshape)
            
            real(kind=real_kind), allocatable:: dlambda(:,:), &
                                                phi1(:,:), &
                                                phi2(:,:), &
                                                dsigma(:,:), &
                                                numerator(:,:), &
                                                denominator(:,:)

            logical, optional :: debugflag

            allocate(dlambda( xshape, yshape), &
                     phi1( xshape, yshape), &
                     phi2( xshape, yshape), &
                     dsigma( xshape, yshape), &
                     numerator( xshape, yshape), &
                     denominator( xshape, yshape))

            dlambda(:,:) = arr_ULONG(:,:) - center_long
            phi1(:,:) = center_lat
            phi2(:,:) = arr_ULAT(:,:)

            numerator = ( cos(phi2)*sin(dlambda) )**2 + (cos(phi1)*sin(phi2) -sin(phi1)*cos(phi2)*cos(dlambda))**2
            numerator = sqrt(numerator)
            denominator = sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(dlambda)

            dsigma = atan2(numerator, denominator)
            distance = radius * abs(dsigma)

            deallocate(dlambda, phi1, phi2, dsigma, numerator, denominator)

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
            north_cornerindex = j_index + half_north_south_BoxSize

            if (west_cornerindex < 1) west_cornerindex = 1
            if (east_cornerindex > nxu) east_cornerindex = nxu
            if (south_cornerindex < 1) south_cornerindex = 1
            if (north_cornerindex > nyu) north_cornerindex = nyu

            east_west_BoxSize = east_cornerindex - west_cornerindex + 1
            north_south_BoxSize = north_cornerindex - south_cornerindex + 1

            ! if (i_index > 165 .AND. filterlengthInKM > 2000) then
            !     print *, 'i_index, j_index, filterlengthInKM', i_index, j_index, filterlengthInKM
            !     print *, 'east_west_BoxSize, north_south_BoxSize', east_west_BoxSize, north_south_BoxSize
            !     stop 'print ERROR in east_west_BoxSize'
            ! endif

            ! if (north_south_BoxSize > nyu) then
            !     print *, 'i_index, j_index, filterlengthInKM', i_index, j_index, filterlengthInKM
            !     print *, 'east_west_BoxSize, north_south_BoxSize', east_west_BoxSize, north_south_BoxSize
            !     stop 'print ERROR in north_south_BoxSize'
            ! endif

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
