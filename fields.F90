module fields
    use kinds
    use gridModule
    use configurationMod
    use mpiwrapper
    implicit none

    type fieldInfo
        character (len=varname_len) :: varname
        character (len=units_len) :: units
        character (len=longname_len) :: long_name
        integer :: ndims
        character (len=dimname_len), allocatable, dimension(:) :: arr_dimnames
    end type

    integer(kind=int_kind) :: num_scalar2D_fields, &
                        &     num_scalar3D_fields, &
                        &     num_vector2D_fields, &
                        &     num_vector3D_fields

    real(kind=real_kind), allocatable, dimension(:,:,:):: scalar2D_fields ! x, y,  fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_scalar2D_fields  ! x, y,  fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: scalar3D_fields ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_scalar3D_fields  ! x, y, z, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DX_fields  ! x, y, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DY_fields  ! x, y, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_vector2DX_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_vector2DY_fields  ! x, y, z, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:):: phi2D_fields  ! x, y, fieldid ! poloidal fields
    real(kind=real_kind), allocatable, dimension(:,:,:):: psi2D_fields  ! x, y, fieldid ! toroidal fields

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_phi2D_fields  ! x, y, fieldid, ell ! poloidal fields
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_psi2D_fields  ! x, y, fieldid, ell ! toroidal fields

    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DX_phi_fields  ! x, y,  fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DY_phi_fields  ! x, y,  fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_vector2DX_phi_fields  ! x, y, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_vector2DY_phi_fields  ! x, y, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DX_psi_fields  ! x, y, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DY_psi_fields  ! x, y, fieldid
    
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_vector2DX_psi_fields  ! x, y, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: OL_vector2DY_psi_fields  ! x, y, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DX_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DY_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DZ_fields  ! x, y, z, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: phi3D_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: psi3D_fields  ! x, y, z, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_phi3D_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_psi3D_fields  ! x, y, z, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DX_phi_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DY_phi_fields  ! x, y, z, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DX_phi_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DY_phi_fields  ! x, y, z, fieldid, ell

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DX_psi_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DY_psi_fields  ! x, y, z, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DX_psi_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DY_psi_fields  ! x, y, z, fieldid, ell
    
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DX_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DY_fields  ! x, y, z, fieldid, ell
    real(kind=real_kind), allocatable, dimension(:,:,:,:,:):: OL_vector3DZ_fields  ! x, y, z, fieldid, ell

    type(fieldInfo), allocatable, dimension(:) :: scalar2D_fields_info
    type(fieldInfo), allocatable, dimension(:) :: scalar3D_fields_info

    type(fieldInfo), allocatable, dimension(:) :: vector2DX_fields_info
    type(fieldInfo), allocatable, dimension(:) :: vector2DY_fields_info

    type(fieldInfo), allocatable, dimension(:) :: phi2D_fields_info
    type(fieldInfo), allocatable, dimension(:) :: psi2D_fields_info

    type(fieldInfo), allocatable, dimension(:) :: vector2DX_phi_fields_info
    type(fieldInfo), allocatable, dimension(:) :: vector2DY_phi_fields_info

    type(fieldInfo), allocatable, dimension(:) :: vector2DX_psi_fields_info
    type(fieldInfo), allocatable, dimension(:) :: vector2DY_psi_fields_info

    type(fieldInfo), allocatable, dimension(:) :: vector3DX_fields_info
    type(fieldInfo), allocatable, dimension(:) :: vector3DY_fields_info
    type(fieldInfo), allocatable, dimension(:) :: vector3DZ_fields_info

    type(fieldInfo), allocatable, dimension(:) :: phi3D_fields_info
    type(fieldInfo), allocatable, dimension(:) :: psi3D_fields_info

    type(fieldInfo), allocatable, dimension(:) :: vector3DX_phi_fields_info
    type(fieldInfo), allocatable, dimension(:) :: vector3DY_phi_fields_info

    type(fieldInfo), allocatable, dimension(:) :: vector3DX_psi_fields_info
    type(fieldInfo), allocatable, dimension(:) :: vector3DY_psi_fields_info

    contains

        subroutine init_unfilt_fields()
            if (taskid == MASTER) then
                call set_num_scalar2D_fields(config%num_of_scalar2D_fields_to_read)
                call set_num_scalar3D_fields(config%num_of_scalar3D_fields_to_read)
                call set_num_vector2D_fields(config%num_of_vector2D_fields_to_read)
                call set_num_vector3D_fields(config%num_of_vector3D_fields_to_read)

                call set_fieldnames(config%list_scalar2D_fieldsNames, config%list_scalar3D_fieldsNames, &
                                &   config%list_vector2DX_fieldsNames,  config%list_vector2DY_fieldsNames, &
                                &   config%list_vector3DX_fieldsNames,  config%list_vector3DY_fieldsNames, config%list_vector3DZ_fieldsNames)
            endif

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call MPI_BCAST(num_scalar2D_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(num_scalar3D_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(num_vector2D_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(num_vector3D_fields, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
    
            call allocate_scalar2D_fields(nxu, nyu)
            call allocate_scalar3D_fields(nxu, nyu, nzu)

            call allocate_vector2D_fields(nxu, nyu)
            call allocate_phipsi2D_fields(nxu, nyu)

            call allocate_vector3D_fields(nxu, nyu, nzu)
            call allocate_phipsi3D_fields(nxu, nyu, nzu)
            
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

        end subroutine

        subroutine set_num_scalar2D_fields(n)
            integer(kind=int_kind), intent(in) :: n
            num_scalar2D_fields = n
        end subroutine

        subroutine set_num_scalar3D_fields(n)
            integer(kind=int_kind), intent(in) :: n
            num_scalar3D_fields = n
        end subroutine

        subroutine set_num_vector2d_fields(n)
            integer(kind=int_kind), intent(in) :: n
            num_vector2D_fields = n
        end subroutine

        subroutine set_num_vector3d_fields(n)
            integer(kind=int_kind), intent(in) :: n
            num_vector3D_fields = n
        end subroutine

        subroutine allocate_scalar2D_fields(nx, ny)
            integer(kind=int_kind), intent(in) :: nx, ny
            if (.not. allocated(scalar2D_fields) .AND. num_scalar2D_fields > 0) then 
                allocate(scalar2D_fields(nx, ny, num_scalar2D_fields))
            endif

            if (.not. allocated(scalar2D_fields_info) .AND. num_scalar2D_fields > 0) then 
                allocate(scalar2D_fields_info(num_scalar2D_fields))
            endif
                
        end subroutine

        subroutine allocate_scalar3D_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            if (.not. allocated(scalar3D_fields) .AND. num_scalar3D_fields > 0) then 
                allocate(scalar3D_fields(nx, ny, nz, num_scalar3D_fields))
            endif

            if (.not. allocated(scalar3D_fields_info) .AND. num_scalar3D_fields > 0) then 
                allocate(scalar3D_fields_info(num_scalar3D_fields))
            endif
                
        end subroutine

        subroutine allocate_vector2D_fields(nx, ny)
            integer(kind=int_kind), intent(in) :: nx, ny
            if (.not. allocated(vector2DX_fields) .AND. num_vector2D_fields > 0) then 
                allocate(vector2DX_fields(nx, ny, num_vector2D_fields), &
                    &    vector2DY_fields(nx, ny, num_vector2D_fields))
            endif

            if (.not. allocated(vector2DX_fields_info) .AND. num_vector2D_fields > 0) then 
                allocate(vector2DX_fields_info(num_vector2D_fields), &
                    &    vector2DY_fields_info(num_vector2D_fields))
            endif
        end subroutine

        subroutine allocate_vector3D_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            if (.not. allocated(vector3DX_fields) .AND. num_vector3D_fields > 0) then 
                allocate(vector3DX_fields(nx, ny, nz, num_vector3D_fields), &
                    &    vector3DY_fields(nx, ny, nz, num_vector3D_fields), &
                    &    vector3DZ_fields(nx, ny, nz, num_vector3D_fields))
            endif

            if (.not. allocated(vector3DX_fields_info) .AND. num_vector3D_fields > 0) then 
                allocate(vector3DX_fields_info(num_vector3D_fields), &
                    &    vector3DY_fields_info(num_vector3D_fields), &
                    &    vector3DZ_fields_info(num_vector3D_fields))
            endif
            
        end subroutine

        subroutine allocate_phipsi2D_fields(nx, ny)
            integer(kind=int_kind), intent(in) :: nx, ny
            if (.not. allocated(phi2D_fields) .AND. num_vector2D_fields > 0) then 
                allocate(phi2D_fields(nx, ny, num_vector2D_fields), &
                    &    psi2D_fields(nx, ny, num_vector2D_fields))
            endif

            if (.not. allocated(psi2D_fields_info) .AND. num_vector2D_fields > 0) then 
                allocate(phi2D_fields_info(num_vector2D_fields), &
                    &    psi2D_fields_info(num_vector2D_fields))
            endif

            if (.not. allocated(vector2DX_psi_fields) .AND. num_vector2D_fields > 0) then 
                allocate(vector2DX_psi_fields(nx, ny, num_vector2D_fields), &
                    &    vector2DY_psi_fields(nx, ny, num_vector2D_fields))
            endif

            if (.not. allocated(vector2DX_phi_fields) .AND. num_vector2D_fields > 0) then 
                allocate(vector2DX_phi_fields(nx, ny, num_vector2D_fields), &
                    &    vector2DY_phi_fields(nx, ny, num_vector2D_fields))
            endif
        end subroutine


        subroutine allocate_phipsi3D_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            if (.not. allocated(phi3D_fields) .AND. num_vector3D_fields > 0) then 
                allocate(phi3D_fields(nx, ny, nz, num_vector3D_fields), &
                    &    psi3D_fields(nx, ny, nz, num_vector3D_fields))
            endif

            if (.not. allocated(psi3D_fields_info) .AND. num_vector3D_fields > 0) then 
                allocate(phi3D_fields_info(num_vector3D_fields), &
                    &    psi3D_fields_info(num_vector3D_fields))
            endif

            if (.not. allocated(vector3DX_psi_fields) .AND. num_vector3D_fields > 0) then 
                allocate(vector3DX_psi_fields(nx, ny, nz, num_vector3D_fields), &
                    &    vector3DY_psi_fields(nx, ny, nz, num_vector3D_fields))
            endif

            if (.not. allocated(vector3DX_phi_fields) .AND. num_vector3D_fields > 0) then 
                allocate(vector3DX_phi_fields(nx, ny, nz, num_vector3D_fields), &
                    &    vector3DY_phi_fields(nx, ny, nz, num_vector3D_fields))
            endif
        end subroutine


        subroutine allocate_filtered_fields()
            integer(kind=int_kind) :: nx, ny, nz, nell

            nx = nxu 
            ny = nyu
            nz = nzu
            nell = num_filterlengths

            if (.not. allocated(OL_scalar2D_fields) .AND. num_scalar2D_fields > 0) then 
                allocate(OL_scalar2D_fields(nx, ny, num_scalar2D_fields, nell))
            endif

            if (.not. allocated(OL_scalar3D_fields) .AND. num_scalar3D_fields > 0) then 
                allocate(OL_scalar3D_fields(nx, ny, nz, num_scalar3D_fields, nell))
            endif

            if (.not. allocated(OL_vector2DX_fields) .AND. num_vector2D_fields > 0) then 
                allocate(OL_vector2DX_fields(nx, ny, num_vector2D_fields, nell), &
                    &    OL_vector2DY_fields(nx, ny, num_vector2D_fields, nell))

                allocate(OL_phi2D_fields(nx, ny, num_vector2D_fields, nell), &
                    &    OL_psi2D_fields(nx, ny, num_vector2D_fields, nell))

                allocate(OL_vector2DX_phi_fields(nx, ny, num_vector2D_fields, nell), &
                    &    OL_vector2DY_phi_fields(nx, ny, num_vector2D_fields, nell))

                allocate(OL_vector2DX_psi_fields(nx, ny, num_vector2D_fields, nell), &
                    &    OL_vector2DY_psi_fields(nx, ny, num_vector2D_fields, nell))
            endif

            if (.not. allocated(OL_vector3DX_fields) .AND. num_vector3D_fields > 0) then 
                allocate(OL_vector3DX_fields(nx, ny, nz, num_vector3D_fields, nell), &
                    &    OL_vector3DY_fields(nx, ny, nz, num_vector3D_fields, nell), &
                    &    OL_vector3DZ_fields(nx, ny, nz, num_vector3D_fields, nell))

                allocate(OL_phi3D_fields(nx, ny, nz, num_vector3D_fields, nell), &
                    &    OL_psi3D_fields(nx, ny, nz, num_vector3D_fields, nell))

                allocate(OL_vector3DX_phi_fields(nx, ny, nz, num_vector3D_fields, nell), &
                    &    OL_vector3DY_phi_fields(nx, ny, nz, num_vector3D_fields, nell))

                allocate(OL_vector3DX_psi_fields(nx, ny, nz, num_vector3D_fields, nell), &
                    &    OL_vector3DY_psi_fields(nx, ny, nz, num_vector3D_fields, nell))
            endif

        end subroutine

        subroutine set_fieldnames(list_scalar2d_fieldnames, list_scalar3d_fieldnames, &
            &                     list_vector2dx_fieldnames,  list_vector2dy_fieldnames, &
            &                     list_vector3dx_fieldnames, list_vector3dy_fieldnames, list_vector3dzfieldnames)
    
            character(len=*), intent(in), dimension(:) :: list_scalar2d_fieldnames, list_scalar3d_fieldnames, &
            &                                             list_vector2dx_fieldnames,  list_vector2dy_fieldnames, &
            &                                             list_vector3dx_fieldnames, list_vector3dy_fieldnames, list_vector3dzfieldnames
    
            integer(kind=int_kind) :: counter
    
    
            do counter = 1, num_scalar2D_fields
                scalar2D_fields_info(counter)%varname = trim(adjustl(list_scalar2d_fieldnames(counter)))
            end do

            do counter = 1, num_scalar3D_fields
                scalar3D_fields_info(counter)%varname = trim(adjustl(list_scalar3d_fieldnames(counter)))
            end do
    
            do counter = 1, num_vector2D_fields
                vector2DX_fields_info(counter)%varname = trim(adjustl(list_vector2dx_fieldnames(counter)))
                vector2DY_fields_info(counter)%varname = trim(adjustl(list_vector2dy_fieldnames(counter)))
            end do
    
            do counter = 1, num_vector3D_fields
                vector3DX_fields_info(counter)%varname = trim(adjustl(list_vector3dx_fieldnames(counter)))
                vector3DY_fields_info(counter)%varname = trim(adjustl(list_vector3dy_fieldnames(counter)))
                vector3DZ_fields_info(counter)%varname = trim(adjustl(list_vector3dzfieldnames(counter)))
            end do    
            
        end subroutine 

        subroutine broadCastReadFields()
            integer :: counter
            do counter = 1, num_scalar2D_fields
                call MPI_BCAST(scalar2D_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do

            do counter = 1, num_scalar3D_fields
                call MPI_BCAST(scalar3D_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do
    
            do counter = 1, num_vector2D_fields
                call MPI_BCAST(vector2DX_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector2DY_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do
    
            do counter = 1, num_vector3D_fields
                call MPI_BCAST(vector3DX_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector3DY_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector3DZ_fields(:,:,:, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do 
            
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
        end subroutine

        subroutine broadCastPsiPhiFields()
            integer :: counter
            do counter = 1, num_vector2D_fields
                call MPI_BCAST(phi2D_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(psi2D_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

                call MPI_BCAST(vector2DX_phi_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector2DX_psi_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

                call MPI_BCAST(vector2DY_phi_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector2DY_psi_fields(:,:, counter), nxu * nyu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            end do

            do counter = 1, num_vector3D_fields
                call MPI_BCAST(phi3D_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(psi3D_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

                call MPI_BCAST(vector3DX_phi_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector3DX_psi_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)

                call MPI_BCAST(vector3DY_phi_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
                call MPI_BCAST(vector3DY_psi_fields(:,:, :, counter), nxu * nyu * nzu, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
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

                do counter = 1, num_scalar2D_fields
                    send_buffer(:,:) = OL_scalar2D_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                        &           OL_scalar2D_fields(:,  :  , counter, filter_index), recv_size, displacements, MPI_REAL, &
                        &           MASTER, MPI_COMM_WORLD, i_err)
                    call MPI_Barrier(MPI_COMM_WORLD, i_err)    
                end do

                do counter = 1, num_scalar3D_fields
                    do zcounter = 1, nzu
                        send_buffer(:,:) = OL_scalar3D_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_scalar3D_fields(:,  :  ,zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                        call MPI_Barrier(MPI_COMM_WORLD, i_err)
                    end do
                    
                end do

                
                do counter = 1, num_vector2D_fields
                    send_buffer(:,:) = OL_vector2DX_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_vector2DX_fields(:,:, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                    
                    send_buffer(:,:) = OL_vector2DY_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_vector2DY_fields(:,:, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                end do

                call MPI_Barrier(MPI_COMM_WORLD, i_err)
        
                do counter = 1, num_vector3D_fields
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

        subroutine collectCoarseGrainedFields(arr_numcols_inallprocs, arr_startcolindex_inallprocs, num_filterlengths)
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

                do counter = 1, num_scalar2D_fields
                    send_buffer(:,:) = OL_scalar2D_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                        &           OL_scalar2D_fields(:,  :  , counter, filter_index), recv_size, displacements, MPI_REAL, &
                        &           MASTER, MPI_COMM_WORLD, i_err)
                    call MPI_Barrier(MPI_COMM_WORLD, i_err)
                end do

                do counter = 1, num_scalar3D_fields
                    do zcounter = 1, nzu
                        send_buffer(:,:) = OL_scalar3D_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_scalar3D_fields(:,  :  ,zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                        call MPI_Barrier(MPI_COMM_WORLD, i_err)
                    end do
                    
                end do
                
                do counter = 1, num_vector2D_fields
                    send_buffer(:,:) = OL_vector2DX_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_vector2DX_fields(:,:, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                    
                    send_buffer(:,:) = OL_vector2DY_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_vector2DY_fields(:,:, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)

                    send_buffer(:,:) = OL_phi2D_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_phi2D_fields(:,:, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                    
                    send_buffer(:,:) = OL_psi2D_fields(:,js:je, counter, filter_index)
                    call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                            &           OL_psi2D_fields(:,:, counter, filter_index), recv_size, displacements, MPI_REAL, &
                            &           MASTER, MPI_COMM_WORLD, i_err)
                end do

                call MPI_Barrier(MPI_COMM_WORLD, i_err)
        
                do counter = 1, num_vector3D_fields
                    do zcounter = 1, nzu
                        send_buffer(:,:) = OL_phi3D_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                                &           OL_phi3D_fields(:,:, zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                                &           MASTER, MPI_COMM_WORLD, i_err)
                        
                        send_buffer(:,:) = OL_psi3D_fields(:,js:je, zcounter, counter, filter_index)
                        call MPI_GATHERV(send_buffer, send_size, MPI_REAL , &
                                &           OL_psi3D_fields(:,:, zcounter, counter, filter_index), recv_size, displacements, MPI_REAL, &
                                &           MASTER, MPI_COMM_WORLD, i_err)

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
            !! SCALAR FIELDS !!
            if (allocated(scalar2D_fields)) deallocate(scalar2D_fields)
            if (allocated(scalar3D_fields)) deallocate(scalar3D_fields)

            if (allocated(scalar2D_fields_info)) deallocate(scalar2D_fields_info)
            if (allocated(scalar3D_fields_info)) deallocate(scalar3D_fields_info)

            !! 2D VECTOR FIELDS !!

            if (allocated(vector2DX_fields)) deallocate(vector2DX_fields)
            if (allocated(vector2DY_fields)) deallocate(vector2DY_fields)

            if (allocated(vector2DX_fields_info)) deallocate(vector2DX_fields_info) 
            if (allocated(vector2DY_fields_info)) deallocate(vector2DY_fields_info)

            if (allocated(vector2DX_phi_fields)) deallocate(vector2DX_phi_fields)
            if (allocated(vector2DY_phi_fields)) deallocate(vector2DY_phi_fields)

            if (allocated(vector2DX_phi_fields_info)) deallocate(vector2DX_phi_fields_info) 
            if (allocated(vector2DY_phi_fields_info)) deallocate(vector2DY_phi_fields_info)

            if (allocated(vector2DX_psi_fields)) deallocate(vector2DX_psi_fields)
            if (allocated(vector2DY_psi_fields)) deallocate(vector2DY_psi_fields)

            if (allocated(vector2DX_psi_fields_info)) deallocate(vector2DX_psi_fields_info) 
            if (allocated(vector2DY_psi_fields_info)) deallocate(vector2DY_psi_fields_info)

            if (allocated(psi2D_fields)) deallocate(psi2D_fields)
            if (allocated(phi2D_fields)) deallocate(phi2D_fields)

            if (allocated(psi2D_fields_info)) deallocate(psi2D_fields_info) 
            if (allocated(phi2D_fields_info)) deallocate(phi2D_fields_info)

            !! 2D VECTOR FIELDS !!

            if (allocated(vector3DX_fields)) deallocate(vector3DX_fields)
            if (allocated(vector3DY_fields)) deallocate(vector3DY_fields)
            if (allocated(vector3DZ_fields)) deallocate(vector3DZ_fields)

            if (allocated(vector3DX_fields_info)) deallocate(vector3DX_fields_info) 
            if (allocated(vector3DY_fields_info)) deallocate(vector3DY_fields_info) 
            if (allocated(vector3DZ_fields_info)) deallocate(vector3DZ_fields_info)

            if (allocated(vector3DX_phi_fields)) deallocate(vector3DX_phi_fields)
            if (allocated(vector3DY_phi_fields)) deallocate(vector3DY_phi_fields)

            if (allocated(vector3DX_phi_fields_info)) deallocate(vector3DX_phi_fields_info) 
            if (allocated(vector3DY_phi_fields_info)) deallocate(vector3DY_phi_fields_info)

            if (allocated(vector3DX_psi_fields)) deallocate(vector3DX_psi_fields)
            if (allocated(vector3DY_psi_fields)) deallocate(vector3DY_psi_fields)

            if (allocated(vector3DX_psi_fields_info)) deallocate(vector3DX_psi_fields_info) 
            if (allocated(vector3DY_psi_fields_info)) deallocate(vector3DY_psi_fields_info)

            if (allocated(psi3D_fields)) deallocate(psi3D_fields)
            if (allocated(phi3D_fields)) deallocate(phi3D_fields)

            if (allocated(psi3D_fields_info)) deallocate(psi3D_fields_info) 
            if (allocated(phi3D_fields_info)) deallocate(phi3D_fields_info)

            !! FILTERED FIELDS !!
            !! FILTERED SCALAR FIELDS !!

            if (allocated(OL_scalar2D_fields)) deallocate(OL_scalar2D_fields)
            if (allocated(OL_scalar3D_fields)) deallocate(OL_scalar3D_fields)

            !! FILTERED 2D VECTOR FIELDS !!

            if (allocated(OL_vector2DX_fields)) deallocate(OL_vector2DX_fields)
            if (allocated(OL_vector2DY_fields)) deallocate(OL_vector2DY_fields)

            if (allocated(OL_phi2D_fields)) deallocate(OL_phi2D_fields)
            if (allocated(OL_psi2D_fields)) deallocate(OL_psi2D_fields)

            if (allocated(OL_vector2DX_phi_fields)) deallocate(OL_vector2DX_phi_fields)
            if (allocated(OL_vector2DY_phi_fields)) deallocate(OL_vector2DY_phi_fields)

            if (allocated(OL_vector2DX_psi_fields)) deallocate(OL_vector2DX_psi_fields)
            if (allocated(OL_vector2DY_psi_fields)) deallocate(OL_vector2DY_psi_fields)

            !! FILTERED 3D VECTOR FIELDS !!

            if (allocated(OL_vector3DX_fields)) deallocate(OL_vector3DX_fields)
            if (allocated(OL_vector3DY_fields)) deallocate(OL_vector3DY_fields)
            if (allocated(OL_vector3DZ_fields)) deallocate(OL_vector3DZ_fields)

            if (allocated(OL_vector3DX_phi_fields)) deallocate(OL_vector3DX_phi_fields)
            if (allocated(OL_vector3DY_phi_fields)) deallocate(OL_vector3DY_phi_fields)

            if (allocated(OL_vector3DX_psi_fields)) deallocate(OL_vector3DX_psi_fields)
            if (allocated(OL_vector3DY_psi_fields)) deallocate(OL_vector3DY_psi_fields)

            if (allocated(OL_phi3D_fields)) deallocate(OL_phi3D_fields)
            if (allocated(OL_psi3D_fields)) deallocate(OL_psi3D_fields)

        end subroutine

end module
