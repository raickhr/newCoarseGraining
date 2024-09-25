module fields
    use kinds
    use gridModule
    implicit none

    type fieldInfo
        character (len=50) :: fileName
        character (len=100) :: units
        character (len=200) :: long_name      
    end type

    integer(kind=int_kind) :: num_scalar_fields, &
                        &     num_2Dvector_fields, &
                        &     num_3Dvector_fields

    real(kind=real_kind), allocatable, dimension(:,:,:,:):: scalar_fields  ! x, y, z, fieldid

    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DX_fields  ! x, y, z, fieldid
    real(kind=real_kind), allocatable, dimension(:,:,:):: vector2DY_fields  ! x, y, z, fieldid

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
            allocate(scalar_fields(nx, ny, nz, num_scalar_fields))
            allocate(scalar_field_info(num_scalar_fields))
        end subroutine

        subroutine allocate_vector2D_fields(nx, ny)
            integer(kind=int_kind), intent(in) :: nx, ny
            allocate(vector2DX_fields(nx, ny, num_2Dvector_fields), &
                &    vector2DY_fields(nx, ny, num_2Dvector_fields))
            allocate(vector2DX_field_info(num_2Dvector_fields), &
                &    vector2DY_field_info(num_2Dvector_fields))
        end subroutine

        subroutine allocate_vector3D_fields(nx, ny, nz)
            integer(kind=int_kind), intent(in) :: nx, ny, nz
            allocate(vector3DX_fields(nx, ny, nz, num_3Dvector_fields), &
                &    vector3DY_fields(nx, ny, nz, num_3Dvector_fields), &
                &    vector3DZ_fields(nx, ny, nz, num_3Dvector_fields)) 
            allocate(vector3DX_field_info(num_3Dvector_fields), &
                &    vector3DY_field_info(num_3Dvector_fields), &
                &    vector3DZ_field_info(num_3Dvector_fields)) 
        end subroutine

end module
