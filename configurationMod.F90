module configurationMod
    use kinds
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  THIS MODULE IS FOR SETTING THE CONFIGURATION TO RUN THE CODE
    !         A CONFIGURATION FILE IS PROVIDED AT RUNTIME
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !public
        type, public:: configuration
            !Defining the configuration elements as a Data type
            character (len=pathname_len) :: InputPath  
            character (len=filename_len) :: gridFile  !Name of the grid files

            integer :: num_of_files_to_read   !Number of input flies excluding gridfile
            character (len = filename_len), allocatable :: list_filenames(:)

            integer :: num_of_scalar_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable :: list_scalar_fieldsNames(:)

            integer :: num_of_2Dvector_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable :: list_2DvectorX_fieldsNames(:)
            character (len = varname_len), allocatable :: list_2DvectorY_fieldsNames(:)
            
            integer :: num_of_3Dvector_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable :: list_3DvectorX_fieldsNames(:)
            character (len = varname_len), allocatable :: list_3DvectorY_fieldsNames(:)
            character (len = varname_len), allocatable :: list_3DvectorZ_fieldsNames(:)

            character (len = varname_len) :: timevar_name
            integer :: startTimeIndex   !Start Time index for each file
            integer :: endTimeIndex   !End Time index for each file
            integer :: nx, ny, nz, nt        ! size of the array in each file
            integer, allocatable :: list_zlevels(:)  ! list of zlevels to read

            integer :: nfilter               ! number of filterlengths
            real, allocatable :: list_filterLength(:)   ! array of the filterlength

            character (len=pathname_len):: OutputPath 
        contains
            procedure :: destruct => deallocate_vars     ! Method to delloacate variables
            procedure :: construct => new_configuration

        end type

        ! !This is for the constructor function for the datatype configuration
        ! interface configuration
        !     module procedure new_configuration
        ! end interface

        contains

        subroutine new_configuration(self)
            class(configuration), intent(inout) :: self

            !Defining the configuration elements as a Data type
            character (len=pathname_len) :: InputPath  
            character (len=filename_len) :: gridFile  !Name of the grid files

            integer :: num_of_files_to_read   !Number of input flies excluding gridfile
            character (len = filename_len), allocatable , dimension(:):: list_filenames

            integer :: num_of_scalar_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable , dimension(:):: list_scalar_fieldsNames

            integer :: num_of_2Dvector_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable , dimension(:):: list_2DvectorX_fieldsNames
            character (len = varname_len), allocatable , dimension(:):: list_2DvectorY_fieldsNames

            integer :: num_of_3Dvector_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable , dimension(:):: list_3DvectorX_fieldsNames
            character (len = varname_len), allocatable , dimension(:):: list_3DvectorY_fieldsNames
            character (len = varname_len), allocatable , dimension(:):: list_3DvectorZ_fieldsNames


            character (len = varname_len) :: timevar_name

            integer :: startTimeIndex   !Start Time index for each file
            integer :: endTimeIndex   !End Time index for each file
            integer :: nx, ny, nz, nt        ! size of the array in each file
            integer, allocatable , dimension(:):: list_zlevels

            integer :: nfilter               ! number of filterlengths
            real, allocatable , dimension(:):: list_filterLength   ! array of the filterlength
            
            character (len=pathname_len):: OutputPath   


            

            integer :: counter
            
            namelist /input/ &
                    & InputPath,&
                    & gridFile,&
                    & num_of_files_to_read, &
                    & num_of_scalar_fields_to_read, &
                    & num_of_2Dvector_fields_to_read, &
                    & num_of_3Dvector_fields_to_read, &
                    & timevar_name, &
                    & startTimeIndex, &
                    & endTimeIndex, &
                    & nx, ny, nz, nt, &
                    & nfilter, &
                    & OutputPath

            namelist /Lists/ &
                    & list_filenames, &
                    & list_scalar_fieldsNames, &
                    & list_2DvectorX_fieldsNames, &
                    & list_2DvectorY_fieldsNames, &
                    & list_3DvectorX_fieldsNames, &
                    & list_3DvectorY_fieldsNames, &
                    & list_3DvectorZ_fieldsNames, &
                    & list_zlevels, &
                    & list_filterLength


            print *, 'reading configuration...'
            read(*,input)

            WRITE(*,'(A50, A)') 'InputPath :', trim(adjustl(InputPath))
            WRITE(*,'(A50, A)') 'gridFile :', trim(adjustl(gridFile))
            WRITE(*,'(A50, I4)') 'num_of_files_to_read :', num_of_files_to_read
            WRITE(*,'(A50, I4)') 'num_of_scalar_fields_to_read :', num_of_scalar_fields_to_read
            WRITE(*,'(A50, I4)') 'num_of_2Dvector_fields_to_read :', num_of_2Dvector_fields_to_read
            WRITE(*,'(A50, I4)') 'num_of_3Dvector_fields_to_read :', num_of_3Dvector_fields_to_read
            WRITE(*,'(A50, A30)') 'timevar_name :', timevar_name
            WRITE(*,'(A50, I4)') 'startTimeIndex :', startTimeIndex
            WRITE(*,'(A50, I4)') 'endTimeIndex :', endTimeIndex
            WRITE(*,'(A50, I4, I4, I4, I4)') 'nx, ny, nz, nt :', nx, ny, nz, nt
            WRITE(*,'(A50, I4)') 'nfilter :', nfilter
            WRITE(*,'(A50, A)') 'OutputPath :', trim(adjustl(OutputPath))

            print *, ' '
            
            self%InputPath  = trim(adjustl(InputPath))
            self%OutputPath   = trim(adjustl(OutputPath))
            self%gridFile  = trim(adjustl(gridFile))
            self%timevar_name = trim(adjustl(timevar_name))
            self%startTimeIndex   = startTimeIndex
            self%endTimeIndex   = endTimeIndex
            self%nx = nx
            self%ny = ny
            self%nz = nz
            self%nt = nt    

            self%num_of_files_to_read   = num_of_files_to_read
            self%num_of_scalar_fields_to_read = num_of_scalar_fields_to_read
            self%num_of_2Dvector_fields_to_read = num_of_2Dvector_fields_to_read
            self%num_of_3Dvector_fields_to_read = num_of_3Dvector_fields_to_read
            self%nfilter = nfilter

            allocate(list_filenames(num_of_files_to_read))
            allocate(list_scalar_fieldsNames(num_of_scalar_fields_to_read))
            allocate(list_2DvectorX_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(list_2DvectorY_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(list_3DvectorX_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(list_3DvectorY_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(list_3DvectorZ_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(list_zlevels(nz))
            allocate(list_filterLength(nfilter))

            allocate(self%list_filenames(num_of_files_to_read))
            allocate(self%list_scalar_fieldsNames(num_of_scalar_fields_to_read))
            allocate(self%list_2DvectorX_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(self%list_2DvectorY_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(self%list_3DvectorX_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(self%list_3DvectorY_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(self%list_3DvectorZ_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(self%list_zlevels(nz))
            allocate(self%list_filterLength(nfilter))

            read(*, Lists)
            
            print *, ' Following are the filenames that will be read'
            do counter=1, num_of_files_to_read
                WRITE(*,'(A50)') trim(adjustl(list_filenames(counter)))
                self%list_filenames(counter) = trim(adjustl(list_filenames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the scalar fieldnames that will be read'
            do counter=1, num_of_scalar_fields_to_read
                WRITE(*,'(A25)') trim(adjustl(list_scalar_fieldsNames(counter)))
                self%list_scalar_fieldsNames(counter) = trim(adjustl(list_scalar_fieldsNames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the 2D vector fieldnames that will be read'
            do counter=1, num_of_2Dvector_fields_to_read
                WRITE(*,'(A25, A25)') trim(adjustl(list_2DvectorX_fieldsNames(counter))), trim(adjustl(list_2DvectorY_fieldsNames(counter)))
                self%list_2DvectorX_fieldsNames(counter) = trim(adjustl(list_2DvectorX_fieldsNames(counter)))
                self%list_2DvectorY_fieldsNames(counter) = trim(adjustl(list_2DvectorY_fieldsNames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the 3D vector fieldnames that will be read'
            do counter=1, num_of_3Dvector_fields_to_read
                WRITE(*,'(A25,A25, A25)') trim(adjustl(list_3DvectorX_fieldsNames(counter))), &
                                    &     trim(adjustl(list_3DvectorY_fieldsNames(counter))), &
                                    &     trim(adjustl(list_3DvectorZ_fieldsNames(counter)))
                self%list_3DvectorX_fieldsNames(counter) = trim(adjustl(list_3DvectorX_fieldsNames(counter)))
                self%list_3DvectorY_fieldsNames(counter) = trim(adjustl(list_3DvectorY_fieldsNames(counter)))
                self%list_3DvectorZ_fieldsNames(counter) = trim(adjustl(list_3DvectorZ_fieldsNames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the vertical level indices that will be read'
            do counter=1, nz
                WRITE(*,'(I50)') list_zlevels(counter)
                self%list_zlevels(counter) = list_zlevels(counter)
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the filterlengths ' 
            do counter=1, nfilter
                WRITE(*,'(F50.3)') list_filterLength(counter)
                self%list_filterLength(counter) = list_filterLength(counter)
            end do
            print *, ' '
            print *, ' '

	        print *, "Configuration read and set SUCCESS "

            deallocate(list_filenames)
            deallocate(list_scalar_fieldsNames)
            deallocate(list_2DvectorX_fieldsNames)
            deallocate(list_2DvectorY_fieldsNames)
            deallocate(list_3DvectorX_fieldsNames)
            deallocate(list_3DvectorY_fieldsNames)
            deallocate(list_3DvectorZ_fieldsNames)
            deallocate(list_zlevels)
            deallocate(list_filterLength)

        end subroutine new_configuration

        subroutine deallocate_vars(self)
            class(configuration), intent(inout) :: self
            if (allocated(self%list_filenames)) deallocate(self%list_filenames)
            if (allocated(self%list_scalar_fieldsNames)) deallocate(self%list_scalar_fieldsNames)
            if (allocated(self%list_2DvectorX_fieldsNames)) deallocate(self%list_2DvectorX_fieldsNames)
            if (allocated(self%list_2DvectorY_fieldsNames)) deallocate(self%list_2DvectorY_fieldsNames)
            if (allocated(self%list_3DvectorX_fieldsNames)) deallocate(self%list_3DvectorX_fieldsNames)
            if (allocated(self%list_3DvectorY_fieldsNames)) deallocate(self%list_3DvectorY_fieldsNames)
            if (allocated(self%list_3DvectorZ_fieldsNames)) deallocate(self%list_3DvectorZ_fieldsNames)
            if (allocated(self%list_zlevels)) deallocate(self%list_zlevels)
            if (allocated(self%list_filterLength)) deallocate(self%list_filterLength)

        end subroutine

end module
