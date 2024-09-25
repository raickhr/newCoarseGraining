module configurationMod

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !  THIS MODULE IS FOR SETTING THE CONFIGURATION TO RUN THE CODE
    !         A CONFIGURATION FILE IS PROVIDED AT RUNTIME
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !public
        type configuration
            !Defining the configuration elements as a Data type
            character (len=250) :: InputPath  
            character (len=150) :: gridFile  !Name of the grid files

            integer :: num_of_files_to_read   !Number of input flies excluding gridfile
            character (len = 150), allocatable :: list_filenames(:)

            integer :: num_of_scalar_fields_to_read !Number of fields to filter
            character (len = 150), allocatable :: list_scalar_fieldsNames(:)

            integer :: num_of_2Dvector_fields_to_read !Number of fields to filter
            character (len = 150), allocatable :: list_2DvectorX_fieldsNames(:)
            character (len = 150), allocatable :: list_2DvectorY_fieldsNames(:)
            
            integer :: num_of_3Dvector_fields_to_read !Number of fields to filter
            character (len = 150), allocatable :: list_3DvectorX_fieldsNames(:)
            character (len = 150), allocatable :: list_3DvectorY_fieldsNames(:)
            character (len = 150), allocatable :: list_3DvectorZ_fieldsNames(:)

            integer :: startTimeIndex   !Start Time index for each file
            integer :: endTimeIndex   !End Time index for each file
            integer :: nx, ny, nz, nt        ! size of the array in each file

            integer :: nfilter               ! number of filterlengths
            real, allocatable :: list_filterLength(:)   ! array of the filterlength

            character (len=250):: OutputPath   

        end type

        !This is for the constructor function for the datatype configuration
        interface configuration
            module procedure new_configuration
        end interface

        contains

        function new_configuration()
            !Defining the configuration elements as a Data type
            character (len=250) :: InputPath  
            character (len=150) :: gridFile  !Name of the grid files

            integer :: num_of_files_to_read   !Number of input flies excluding gridfile
            character (len = 150), allocatable :: list_filenames(:)

            integer :: num_of_scalar_fields_to_read !Number of fields to filter
            character (len = 150), allocatable :: list_scalar_fieldsNames(:)

            integer :: num_of_2Dvector_fields_to_read !Number of fields to filter
            character (len = 150), allocatable :: list_2DvectorX_fieldsNames(:)
            character (len = 150), allocatable :: list_2DvectorY_fieldsNames(:)

            integer :: num_of_3Dvector_fields_to_read !Number of fields to filter
            character (len = 150), allocatable :: list_3DvectorX_fieldsNames(:)
            character (len = 150), allocatable :: list_3DvectorY_fieldsNames(:)
            character (len = 150), allocatable :: list_3DvectorZ_fieldsNames(:)

            integer :: startTimeIndex   !Start Time index for each file
            integer :: endTimeIndex   !End Time index for each file
            integer :: nx, ny, nz, nt        ! size of the array in each file

            integer :: nfilter               ! number of filterlengths
            real, allocatable :: list_filterLength(:)   ! array of the filterlength
            
            character (len=250):: OutputPath   


            type(configuration) :: new_configuration

            integer :: counter
            
            namelist /input/ &
                    & InputPath,&
                    & gridFile,&
                    & num_of_files_to_read, &
                    & num_of_scalar_fields_to_read, &
                    & num_of_2Dvector_fields_to_read, &
                    & num_of_3Dvector_fields_to_read, &
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
                    & list_filterLength


            print *, 'reading configuration...'
            read(*,input)

            WRITE(*,'(A50, A)') 'InputPath :', InputPath
            WRITE(*,'(A50, A)') 'gridFile :', gridFile
            WRITE(*,'(A50, I4)') 'num_of_files_to_read :', num_of_files_to_read
            WRITE(*,'(A50, I4)') 'num_of_scalar_fields_to_read :', num_of_scalar_fields_to_read
            WRITE(*,'(A50, I4)') 'num_of_2Dvector_fields_to_read :', num_of_2Dvector_fields_to_read
            WRITE(*,'(A50, I4)') 'num_of_3Dvector_fields_to_read :', num_of_3Dvector_fields_to_read
            WRITE(*,'(A50, I4)') 'startTimeIndex :', startTimeIndex
            WRITE(*,'(A50, I4)') 'endTimeIndex :', endTimeIndex
            WRITE(*,'(A50, I4, I4, I4, I4)') 'nx, ny, nz, nt :', nx, ny, nz, nt
            WRITE(*,'(A50, I4)') 'nfilter :', nfilter
            WRITE(*,'(A50, A)') 'OutputPath :', OutputPath

            print *, ' '
            
            new_configuration%InputPath  = trim(adjustl(InputPath))
            new_configuration%OutputPath   = trim(adjustl(OutputPath))
            new_configuration%gridFile  = trim(adjustl(gridFile))
            new_configuration%startTimeIndex   = startTimeIndex
            new_configuration%endTimeIndex   = endTimeIndex
            new_configuration%nx = nx
            new_configuration%ny = ny
            new_configuration%nz = nz
            new_configuration%nt = nt    

            new_configuration%num_of_files_to_read   = num_of_files_to_read
            new_configuration%num_of_scalar_fields_to_read = num_of_scalar_fields_to_read
            new_configuration%num_of_2Dvector_fields_to_read = num_of_2Dvector_fields_to_read
            new_configuration%num_of_3Dvector_fields_to_read = num_of_3Dvector_fields_to_read
            new_configuration%nfilter = nfilter

            allocate(list_filenames(num_of_files_to_read))
            allocate(list_scalar_fieldsNames(num_of_scalar_fields_to_read))
            allocate(list_2DvectorX_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(list_2DvectorY_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(list_3DvectorX_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(list_3DvectorY_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(list_3DvectorZ_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(list_filterLength(nfilter))

            allocate(new_configuration%list_filenames(num_of_files_to_read))
            allocate(new_configuration%list_scalar_fieldsNames(num_of_scalar_fields_to_read))
            allocate(new_configuration%list_2DvectorX_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(new_configuration%list_2DvectorY_fieldsNames(num_of_2Dvector_fields_to_read))
            allocate(new_configuration%list_3DvectorX_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(new_configuration%list_3DvectorY_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(new_configuration%list_3DvectorZ_fieldsNames(num_of_3Dvector_fields_to_read))
            allocate(new_configuration%list_filterLength(nfilter))

            read(*, Lists)
            
            print *, ' Following are the filenames that will be read'
            do counter=1, num_of_files_to_read
                WRITE(*,'(A50)') trim(adjustl(list_filenames(counter)))
                new_configuration%list_filenames(counter) = trim(adjustl(list_filenames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the scalar fieldnames that will be read'
            do counter=1, num_of_scalar_fields_to_read
                WRITE(*,'(A50)') trim(adjustl(list_scalar_fieldsNames(counter)))
                new_configuration%list_scalar_fieldsNames(counter) = trim(adjustl(list_scalar_fieldsNames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the 2D vector fieldnames that will be read'
            do counter=1, num_of_2Dvector_fields_to_read
                WRITE(*,'(A50, A50)') trim(adjustl(list_2DvectorX_fieldsNames(counter))), trim(adjustl(list_2DvectorY_fieldsNames(counter)))
                new_configuration%list_2DvectorX_fieldsNames(counter) = trim(adjustl(list_2DvectorX_fieldsNames(counter)))
                new_configuration%list_2DvectorY_fieldsNames(counter) = trim(adjustl(list_2DvectorY_fieldsNames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the 3D vector fieldnames that will be read'
            do counter=1, num_of_3Dvector_fields_to_read
                WRITE(*,'(A50,A50, A50)') trim(adjustl(list_3DvectorX_fieldsNames(counter))), &
                                    &     trim(adjustl(list_3DvectorY_fieldsNames(counter))), &
                                    &     trim(adjustl(list_3DvectorZ_fieldsNames(counter)))
                new_configuration%list_3DvectorX_fieldsNames(counter) = trim(adjustl(list_3DvectorX_fieldsNames(counter)))
                new_configuration%list_3DvectorY_fieldsNames(counter) = trim(adjustl(list_3DvectorY_fieldsNames(counter)))
                new_configuration%list_3DvectorZ_fieldsNames(counter) = trim(adjustl(list_3DvectorZ_fieldsNames(counter)))
            end do
            print *, ' '
            print *, ' '

            print *, ' Following are the filterlengths '
            do counter=1, nfilter
                WRITE(*,'(F50.3)') list_filterLength(counter)
                new_configuration%list_filterLength(counter) = list_filterLength(counter)
            end do
            print *, ' '
            print *, ' '

	        print *, "Configuration read and set SUCCESS "

        end function new_configuration

end module
