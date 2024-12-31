module fields
    use kinds
    use gridModule
    use mpiwrapper
    implicit none

    type fieldInfo
        character (len=varname_len) :: varname
        character (len=units_len) :: units
        character (len=longname_len) :: long_name
        integer :: ndims
        character (len=dimname_len), allocatable, dimension(:) :: arr_dimnames
    end type

    integer(kind=int_kind) :: num_scalar_fields, &
                        &     num_2Dvector_fields, &
                        &     num_3Dvector_fields

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: scalar_fields ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_scalar_fields  ! x, y, z, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector2DX_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector2DY_fields  ! x, y, z, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: phi_fields  ! x, y, z, fieldid ! poloidal fields
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: psi_fields  ! x, y, z, fieldid ! toroidal fields

    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector2DX_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector2DY_fields  ! x, y, z, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DX_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DY_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DZ_fields  ! x, y, z, fieldid


    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DX_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DY_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DZ_fields  ! x, y, z, fieldid, ell

    type(fieldInfo), allocatable, dimension(:) :: scalar_field_info

    type(fieldInfo), allocatable, dimension(:) :: vector2DX_field_info
    type(fieldInfo), allocatable, dimension(:) :: vector2DY_field_info

    type(fieldInfo), allocatable, dimension(:) :: phi_field_info
    type(fieldInfo), allocatable, dimension(:) :: psi_field_info

    type(fieldInfo), allocatable, dimension(:) :: vector3DX_field_info
    type(fieldInfo), allocatable, dimension(:) :: vector3DY_field_info
    type(fieldInfo), allocatable, dimension(:) :: vector3DZ_field_info

    contains

        subroutine set_num_scalar_fields(n)
            integer(kind=int_kind), intent(in) :: n
            num_scalar_fields = n
        end subroutine

        subroutine set_num_2dvector_fields(n)
            integer(kind=int_kind), intent(in) :: n
            num_2Dvector_fields = n
        end subroutine

        subroutine set_num_3dvector_fields(n)
            integer(kind=int_kind), intent(in) :: n
            num_3Dvector_fields = n
        end subroutine

        subroutine allocate_scalar_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            if (.not. allocated(scalar_fields) .AND. num_scalar_fields > 0) then 
                allocate(scalar_fields(nx, ny, nz, num_scalar_fields))
            endif

            if (.not. allocated(scalar_field_info) .AND. num_scalar_fields > 0) then 
                allocate(scalar_field_info(num_scalar_fields))
            endif
                
        end subroutine

        subroutine allocate_vector2D_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            if (.not. allocated(vector2DX_fields) .AND. num_2Dvector_fields > 0) then 
                allocate(vector2DX_fields(nx, ny, nz, num_2Dvector_fields), &
                    &    vector2DY_fields(nx, ny, nz, num_2Dvector_fields))
            endif

            if (.not. allocated(vector2DX_field_info) .AND. num_2Dvector_fields > 0) then 
                allocate(vector2DX_field_info(num_2Dvector_fields), &
                    &    vector2DY_field_info(num_2Dvector_fields))
            endif
        end subroutine

        subroutine allocate_phi_psi_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            if (.not. allocated(phi_fields) .AND. num_2Dvector_fields > 0) then 
                allocate(phi_fields(nx, ny, nz, num_2Dvector_fields), &
                    &    psi_fields(nx, ny, nz, num_2Dvector_fields))
            endif

            if (.not. allocated(psi_field_info) .AND. num_2Dvector_fields > 0) then 
                allocate(phi_field_info(num_2Dvector_fields), &
                    &    psi_field_info(num_2Dvector_fields))
            endif
        end subroutine

        subroutine allocate_vector3D_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            if (.not.allocated(vector3DX_fields) .AND. num_3Dvector_fields > 0) then
                allocate(vector3DX_fields(nx, ny, nz, num_3Dvector_fields), &
                    &    vector3DY_fields(nx, ny, nz, num_3Dvector_fields), &
                    &    vector3DZ_fields(nx, ny, nz, num_3Dvector_fields)) 
            end if
            if (.not. allocated(vector3DX_field_info) .AND. num_3Dvector_fields > 0) then 
                allocate(vector3DX_field_info(num_3Dvector_fields), &
                    &    vector3DY_field_info(num_3Dvector_fields), &
                    &    vector3DZ_field_info(num_3Dvector_fields)) 
            end if
        end subroutine

        subroutine allocate_output_fields(nx, ny, nz, nell)
            integer(kind=int_kind), intent(in) :: nx, ny, nz, nell

            if (.not. allocated(OL_scalar_fields) .AND. num_scalar_fields > 0) then 
                allocate(OL_scalar_fields(nx, ny, nz, num_scalar_fields, nell))
            endif

            if (.not. allocated(OL_vector2DX_fields) .AND. num_2Dvector_fields > 0) then 
                allocate(OL_vector2DX_fields(nx, ny, nz, num_2Dvector_fields, nell), &
                    &    OL_vector2DY_fields(nx, ny, nz, num_2Dvector_fields, nell))
            endif

            if (.not. allocated(OL_vector3DX_fields) .AND. num_3Dvector_fields > 0) then 
                allocate(OL_vector3DX_fields(nx, ny, nz, num_3Dvector_fields, nell), &
                    &    OL_vector3DY_fields(nx, ny, nz, num_3Dvector_fields, nell), &
                    &    OL_vector3DZ_fields(nx, ny, nz, num_3Dvector_fields, nell))
            endif

        end subroutine

        subroutine set_fieldnames(list_scalarfieldnames, &
            &                     list_2dvectorxfieldnames,  list_2dvectoryfieldnames, &
            &                     list_3dvectorxfieldnames, list_3dvectoryfieldnames, list_3dvectorzfieldnames)
    
            character(len=*), intent(in), dimension(:) :: list_scalarfieldnames, &
            &                                             list_2dvectorxfieldnames,  list_2dvectoryfieldnames, &
            &                                             list_3dvectorxfieldnames, list_3dvectoryfieldnames, list_3dvectorzfieldnames
    
            integer(kind=int_kind) :: counter
    
    
            do counter = 1, num_scalar_fields
                scalar_field_info(counter)%varname = trim(adjustl(list_scalarfieldnames(counter)))
            end do
    
            do counter = 1, num_2Dvector_fields
                vector2DX_field_info(counter)%varname = trim(adjustl(list_2dvectorxfieldnames(counter)))
                vector2DY_field_info(counter)%varname = trim(adjustl(list_2dvectoryfieldnames(counter)))
            end do
    
            do counter = 1, num_3Dvector_fields
                vector3DX_field_info(counter)%varname = trim(adjustl(list_3dvectorxfieldnames(counter)))
    
                vector3DY_field_info(counter)%varname = trim(adjustl(list_3dvectoryfieldnames(counter)))

                vector3DZ_field_info(counter)%varname = trim(adjustl(list_3dvectorzfieldnames(counter)))
            end do    
            
        end subroutine 

        subroutine broadCastFieldInfo()
            call MPI_BCAST(num_scalar_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(num_2Dvector_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(num_3Dvector_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
    
            if (taskid .NE. MASTER) then
                call allocate_scalar_fields(nxu, nyu, nzu)
                call allocate_vector2D_fields(nxu, nyu, nzu)
                call allocate_phi_psi_fields(nxu, nyu, nzu)
                call allocate_vector3D_fields(nxu, nyu, nzu)
            end if
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
        end subroutine

        subroutine broadCastReadFields()
            integer :: counter
            do counter = 1, num_scalar_fields
                call MPI_BCAST(scalar_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do
    
            do counter = 1, num_2Dvector_fields
                call MPI_BCAST(vector2DX_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector2DY_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do
    
            do counter = 1, num_3Dvector_fields
                call MPI_BCAST(vector3DX_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector3DY_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector3DZ_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do 
            
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
        end subroutine

        subroutine collectFilteredFields(arr_numcols_inallprocs, arr_startcolindex_inallprocs, num_filterlengths)
            integer, intent(in) :: arr_numcols_inallprocs(numtasks), &
                                   arr_startcolindex_inallprocs(numtasks), &
                                   num_filterlengths
            integer :: offset, counter, zcounter, js, je, send_size, recv_size(numtasks), displacements(numtasks), filter_index

            real (kind=real_kind), allocatable :: send_buffer(:,:), recv_buffer(:,:)

            js = arr_startcolindex_inallprocs(taskid + 1)
            je = js + arr_numcols_inallprocs(taskid + 1) -1

            recv_size = arr_numcols_inallprocs * nxu
            send_size = nxu * arr_numcols_inallprocs(taskid + 1)

            displacements(:) = 0
            offset = 0
            do counter = 2, numtasks
                offset = offset + recv_size(counter-1)
                displacements(counter) = offset
            end do
            
            allocate(send_buffer(nxu, arr_numcols_inallprocs(taskid + 1)))
            allocate(recv_buffer(nxu, nyu))
            

            call MPI_Barrier(MPI_COMM_WORLD, i_err)


            do filter_index=1, num_filterlengths

                do counter = 1, num_scalar_fields
                    do zcounter = 1, nzu
                        send_buffer(:,:) = OL_scalar_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_scalar_fields(:,  :  ,zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                        call MPI_Barrier(MPI_COMM_WORLD, i_err)
                    end do
                    
                end do

                
                do counter = 1, num_2Dvector_fields
                    do zcounter = 1, nzu
                        send_buffer(:,:) = OL_vector2DX_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                                &           OL_vector2DX_fields(:,:, zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                                &           MASTER, MPI_COMM_WORLD, i_err)
                        
                        send_buffer(:,:) = OL_vector2DY_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                                &           OL_vector2DY_fields(:,:, zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                                &           MASTER, MPI_COMM_WORLD, i_err)
                    end do
                end do

                call MPI_Barrier(MPI_COMM_WORLD, i_err)
        
                do counter = 1, num_3Dvector_fields
                    do zcounter = 1, nzu
                        send_buffer(:,:) = OL_vector3DX_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_vector3DX_fields(:,:,zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)

                        send_buffer(:,:) = OL_vector3DY_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_vector3DY_fields(:,:,zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)

                        send_buffer(:,:) = OL_vector3DZ_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_vector3DZ_fields(:,:,zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                    end do
                    call MPI_Barrier(MPI_COMM_WORLD, i_err)
                end do 

            end do

            deallocate(send_buffer)
            
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

        end subroutine

        subroutine delloacate_fields()
            if (allocated(scalar_fields)) deallocate(scalar_fields)

            if (allocated(vector2DX_fields)) deallocate(vector2DX_fields)
            if (allocated(vector2DY_fields)) deallocate(vector2DY_fields)

            if (allocated(psi_fields)) deallocate(psi_fields)
            if (allocated(phi_fields)) deallocate(phi_fields)

            if (allocated(vector3DX_fields)) deallocate(vector3DX_fields)
            if (allocated(vector3DY_fields)) deallocate(vector3DY_fields)
            if (allocated(vector3DZ_fields)) deallocate(vector3DZ_fields)

            if (allocated(scalar_field_info)) deallocate(scalar_field_info)

            if (allocated(vector2DX_field_info)) deallocate(vector2DX_field_info) 
            if (allocated(vector2DY_field_info)) deallocate(vector2DY_field_info)

            if (allocated(psi_field_info)) deallocate(psi_field_info) 
            if (allocated(phi_field_info)) deallocate(phi_field_info)

            if (allocated(vector3DX_field_info)) deallocate(vector3DX_field_info) 
            if (allocated(vector3DY_field_info)) deallocate(vector3DY_field_info) 
            if (allocated(vector3DZ_field_info)) deallocate(vector3DZ_field_info)

            if (allocated(OL_scalar_fields)) deallocate(OL_scalar_fields)
            if (allocated(OL_vector2DX_fields)) deallocate(OL_vector2DX_fields)
            if (allocated(OL_vector2DY_fields)) deallocate(OL_vector2DY_fields)
            if (allocated(OL_vector3DX_fields)) deallocate(OL_vector3DX_fields)
            if (allocated(OL_vector3DY_fields)) deallocate(OL_vector3DY_fields)
            if (allocated(OL_vector3DZ_fields)) deallocate(OL_vector3DZ_fields)

        end subroutine

end module
