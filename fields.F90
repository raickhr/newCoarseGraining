module fields
    use kinds
    use gridModule
    use mpiwrapper
    implicit none

    type fieldInfo
        character (len=varname_len) :: varname
        character (len=units_len) :: units
        character (len=longname_len) :: long_name
    end type

    integer(kind=int_kind) :: num_scalar_fields, &
                        &     num_2Dvector_fields, &
                        &     num_3Dvector_fields

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: scalar_fields  ! x, y, z, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DX_fields  ! x, y, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DY_fields  ! x, y, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DX_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DY_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:,:):: vector3DZ_fields  ! x, y, z, fieldid

    type(fieldInfo), allocatable, dimension(:) :: scalar_field_info

    type(fieldInfo), allocatable, dimension(:) :: vector2DX_field_info
    type(fieldInfo), allocatable, dimension(:) :: vector2DY_field_info

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
            if (.not. allocated(scalar_fields) .AND. num_scalar_fields > 0) allocate(scalar_fields(nx, ny, nz, num_scalar_fields))
            if (.not. allocated(scalar_field_info) .AND. num_scalar_fields > 0) allocate(scalar_field_info(num_scalar_fields))
        end subroutine

        subroutine allocate_vector2D_fields(nx, ny)
            integer(kind=int_kind), intent(in) :: nx, ny
            if (.not. allocated(vector2DX_fields) .AND. num_2Dvector_fields > 0) then 
            allocate(vector2DX_fields(nx, ny, num_2Dvector_fields), &
                &    vector2DY_fields(nx, ny, num_2Dvector_fields))
            endif
            if (.not. allocated(vector2DX_field_info) .AND. num_2Dvector_fields > 0) then 
            allocate(vector2DX_field_info(num_2Dvector_fields), &
                &    vector2DY_field_info(num_2Dvector_fields))
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
                call allocate_vector2D_fields(nxu, nyu)
                call allocate_vector3D_fields(nxu, nyu, nzu)
            end if
            call MPI_Barrier(MPI_COMM_WORLD, i_err)
        end subroutine

        subroutine delloacate_fields()
            if (allocated(scalar_fields)) deallocate(scalar_fields)
            if (allocated(vector2DX_fields)) deallocate(vector2DX_fields)
            if (allocated(vector2DY_fields)) deallocate(vector2DY_fields)
            if (allocated(vector3DX_fields)) deallocate(vector3DX_fields)
            if (allocated(vector3DY_fields)) deallocate(vector3DY_fields)
            if (allocated(vector3DZ_fields)) deallocate(vector3DZ_fields)
            if (allocated(scalar_field_info)) deallocate(scalar_field_info)
            if (allocated(vector2DX_field_info)) deallocate(vector2DX_field_info) 
            if (allocated(vector2DY_field_info)) deallocate(vector2DY_field_info) 
            if (allocated(vector3DX_field_info)) deallocate(vector3DX_field_info) 
            if (allocated(vector3DY_field_info)) deallocate(vector3DY_field_info) 
            if (allocated(vector3DZ_field_info)) deallocate(vector3DZ_field_info)
        end subroutine

end module
