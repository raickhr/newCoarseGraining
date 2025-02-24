module configurationMod
    use kinds
    use mpiwrapper
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

            integer :: num_of_scalar2D_fields_to_read !Number of 2D scalar fields to filter
            character (len = varname_len), allocatable :: list_scalar2D_fieldsNames(:)

            integer :: num_of_scalar3D_fields_to_read !Number of 3d scalar fields to filter
            character (len = varname_len), allocatable :: list_scalar3D_fieldsNames(:)

            integer :: num_of_vector2D_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable :: list_vector2DX_fieldsNames(:)
            character (len = varname_len), allocatable :: list_vector2DY_fieldsNames(:)
            
            integer :: num_of_vector3D_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable :: list_vector3DX_fieldsNames(:)
            character (len = varname_len), allocatable :: list_vector3DY_fieldsNames(:)
            character (len = varname_len), allocatable :: list_vector3DZ_fieldsNames(:)

            character (len = varname_len) :: timevar_name
            character (len = varname_len) :: vertdim_name
            character (len = varname_len) :: wvertdim_name
            integer :: startTimeIndex   !Start Time index for each file
            integer :: endTimeIndex   !End Time index for each file
            integer :: nx, ny, nz, nt        ! size of the array in each file
            integer, allocatable :: list_zlevels(:) ! list of zlevels to read
            integer, allocatable :: list_wzlevels(:) ! list of zlevels to read on top and bottom cell

            integer :: ncoarse_levels
            integer, allocatable :: list_coarse_factor_levels(:)

            integer :: max_iterations   ! for helmholtz decomp
            real(kind=real_kind) :: abs_tol, rel_tol, div_tol

            integer :: nfilter               ! number of filterlengths
            real, allocatable :: list_filterLength(:)   ! array of the filterlength

            integer:: num_iti_laplace_smooth  ! number of iteration for laplace smoothing for extending values at land values

            character (len=pathname_len):: OutputPath 
        contains
            procedure :: destruct => deallocate_vars     ! Method to delloacate variables
            procedure :: construct => new_configuration

        end type

        type(configuration):: config! object that defines run configuration

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

            integer :: num_of_scalar2D_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable , dimension(:):: list_scalar2D_fieldsNames

            integer :: num_of_scalar3D_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable , dimension(:):: list_scalar3D_fieldsNames

            integer :: num_of_vector2D_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable , dimension(:):: list_vector2DX_fieldsNames
            character (len = varname_len), allocatable , dimension(:):: list_vector2DY_fieldsNames

            integer :: num_of_vector3D_fields_to_read !Number of fields to filter
            character (len = varname_len), allocatable , dimension(:):: list_vector3DX_fieldsNames
            character (len = varname_len), allocatable , dimension(:):: list_vector3DY_fieldsNames
            character (len = varname_len), allocatable , dimension(:):: list_vector3DZ_fieldsNames


            character (len = varname_len) :: timevar_name
            character (len = varname_len) :: vertdim_name
            character (len = varname_len) :: wvertdim_name

            integer :: startTimeIndex   !Start Time index for each file
            integer :: endTimeIndex   !End Time index for each file
            integer :: nx, ny, nz, nt        ! size of the array in each file

            integer, allocatable :: list_zlevels(:) ! list of zlevels to read
            integer, allocatable :: list_wzlevels(:) ! list of zlevels to read on top and bottom cell

            integer :: ncoarse_levels ! number of coarsening levels for Helmholtz decomposition
            integer, allocatable :: list_coarse_factor_levels(:)

            integer :: max_iterations ! maximum iterations for Helmholtz decomposition
            real(kind=real_kind) :: abs_tol, rel_tol, div_tol ! absolute tolerance, relative tolerance and divergence tolerance

            integer :: nfilter               ! number of filterlengths
            real, allocatable , dimension(:):: list_filterLength   ! array of the filterlength

            integer:: num_iti_laplace_smooth  ! number of iteration for laplace smoothing for extending values at land values
            
            character (len=pathname_len):: OutputPath   


            

            integer :: counter
            
            namelist /input/ &
                    & InputPath,&
                    & gridFile,&
                    & num_of_files_to_read, &
                    & num_of_scalar2D_fields_to_read, &
                    & num_of_scalar3D_fields_to_read, &
                    & num_of_vector2D_fields_to_read, &
                    & num_of_vector3D_fields_to_read, &
                    & timevar_name, &
                    & vertdim_name, &
                    & wvertdim_name, &
                    & startTimeIndex, &
                    & endTimeIndex, &
                    & nx, ny, nz, nt, &
                    & ncoarse_levels, &
                    & max_iterations, &
                    & abs_tol, rel_tol, div_tol, &
                    & nfilter, &
                    & num_iti_laplace_smooth, &
                    & OutputPath

            namelist /fileList/ &
                    & list_filenames

            namelist /scalar2DList/ &
                    & list_scalar2D_fieldsNames

            namelist /scalar3DList/ &
                    & list_scalar3D_fieldsNames

            namelist /vector2DList/ &
                    & list_vector2DX_fieldsNames, &
                    & list_vector2DY_fieldsNames

            namelist /vector3DList/ &
                    & list_vector3DX_fieldsNames, &
                    & list_vector3DY_fieldsNames, &
                    & list_vector3DZ_fieldsNames

            namelist /zlevelLists/ &
                    & list_zlevels, list_wzlevels

            namelist /coarseFactorList/ &
                    & list_coarse_factor_levels
            namelist /filterLenList/ &
                    & list_filterLength


            print *, 'reading configuration...'
            read(*,input)

            WRITE(*,'(A50, A)') 'InputPath :', trim(adjustl(InputPath))
            WRITE(*,'(A50, A)') 'gridFile :', trim(adjustl(gridFile))
            WRITE(*,'(A50, I4)') 'num_of_files_to_read :', num_of_files_to_read
            WRITE(*,'(A50, I4)') 'num_of_scalar2D_fields_to_read :', num_of_scalar2D_fields_to_read
            WRITE(*,'(A50, I4)') 'num_of_scalar3D_fields_to_read :', num_of_scalar3D_fields_to_read
            WRITE(*,'(A50, I4)') 'num_of_vector2D_fields_to_read :', num_of_vector2D_fields_to_read
            WRITE(*,'(A50, I4)') 'num_of_vector3D_fields_to_read :', num_of_vector3D_fields_to_read
            WRITE(*,'(A50, A30)') 'timevar_name :', timevar_name
            WRITE(*,'(A50, A30)') 'vertdim_name :', vertdim_name
            WRITE(*,'(A50, A30)') 'wvertdim_name :', wvertdim_name
            WRITE(*,'(A50, I4)') 'startTimeIndex :', startTimeIndex
            WRITE(*,'(A50, I4)') 'endTimeIndex :', endTimeIndex
            WRITE(*,'(A50, I4, I4, I4, I4)') 'nx, ny, nz, nt :', nx, ny, nz, nt
            WRITE(*,'(A50, I4)') 'ncoarse_levels :', ncoarse_levels
            WRITE(*,'(A50, I4)') 'max_iterations :', max_iterations
            WRITE(*,'(A50, E9.2, E9.2, E9.2)') 'abs_tol, rel_tol, div_tol :', abs_tol, rel_tol, div_tol
            WRITE(*,'(A50, I4)') 'ncoarse_levels :', ncoarse_levels
            WRITE(*,'(A50, I4)') 'nfilter :', nfilter
            WRITE(*,'(A50, I4)') 'iterations laplace smooth :', num_iti_laplace_smooth
            WRITE(*,'(A50, A)') 'OutputPath :', trim(adjustl(OutputPath))

            print *, ' '
            
            self%InputPath  = trim(adjustl(InputPath))
            self%OutputPath   = trim(adjustl(OutputPath))
            self%gridFile  = trim(adjustl(gridFile))
            self%timevar_name = trim(adjustl(timevar_name))
            self%vertdim_name = trim(adjustl(vertdim_name))
            self%wvertdim_name = trim(adjustl(wvertdim_name))
            self%startTimeIndex   = startTimeIndex
            self%endTimeIndex   = endTimeIndex
            self%nx = nx
            self%ny = ny
            self%nz = nz
            self%nt = nt    

            self%num_of_files_to_read   = num_of_files_to_read
            self%num_of_scalar2D_fields_to_read = num_of_scalar2D_fields_to_read
            self%num_of_scalar3D_fields_to_read = num_of_scalar3D_fields_to_read
            self%num_of_vector2D_fields_to_read = num_of_vector2D_fields_to_read
            self%num_of_vector3D_fields_to_read = num_of_vector3D_fields_to_read
            self%ncoarse_levels = ncoarse_levels
            self%max_iterations = max_iterations
            self%abs_tol = abs_tol
            self%rel_tol = rel_tol
            self%div_tol = div_tol
            self%nfilter = nfilter
            self%num_iti_laplace_smooth = num_iti_laplace_smooth

            if ( num_of_files_to_read > 0) then 
                allocate(list_filenames(num_of_files_to_read))
                allocate(self%list_filenames(num_of_files_to_read))
                read(*,fileList)

                print *, ' Following are the filenames that will be read'
                do counter=1, num_of_files_to_read
                    WRITE(*,'(A50)') trim(adjustl(list_filenames(counter)))
                    self%list_filenames(counter) = trim(adjustl(list_filenames(counter)))
                end do
                print *, ' '
                print *, ' '

                deallocate(list_filenames)
            endif

            if ( num_of_scalar2D_fields_to_read > 0) then 
                
                allocate(list_scalar2D_fieldsNames(num_of_scalar2D_fields_to_read))
                allocate(self%list_scalar2D_fieldsNames(num_of_scalar2D_fields_to_read))

                read(*,scalar2DList)

                print *, ' Following are the 2D scalar fieldnames that will be read'
                do counter=1, num_of_scalar2D_fields_to_read
                    WRITE(*,'(A25)') trim(adjustl(list_scalar2D_fieldsNames(counter)))
                    self%list_scalar2D_fieldsNames(counter) = trim(adjustl(list_scalar2D_fieldsNames(counter)))
                end do
                print *, ' '
                print *, ' '
                deallocate(list_scalar2D_fieldsNames)

            endif

            if ( num_of_scalar3D_fields_to_read > 0) then 
               
                allocate(list_scalar3D_fieldsNames(num_of_scalar3D_fields_to_read))
                allocate(self%list_scalar3D_fieldsNames(num_of_scalar3D_fields_to_read))

                read(*,scalar3DList)

                print *, ' Following are the 3D scalar fieldnames that will be read'
                do counter=1, num_of_scalar3D_fields_to_read
                    WRITE(*,'(A25)') trim(adjustl(list_scalar3D_fieldsNames(counter)))
                    self%list_scalar3D_fieldsNames(counter) = trim(adjustl(list_scalar3D_fieldsNames(counter)))
                end do
                print *, ' '
                print *, ' '

                deallocate(list_scalar3D_fieldsNames)
            endif


            if ( num_of_vector2D_fields_to_read > 0) then
                
                allocate(list_vector2DX_fieldsNames(num_of_vector2D_fields_to_read))
                allocate(list_vector2DY_fieldsNames(num_of_vector2D_fields_to_read))

                allocate(self%list_vector2DX_fieldsNames(num_of_vector2D_fields_to_read))
                allocate(self%list_vector2DY_fieldsNames(num_of_vector2D_fields_to_read))

                read(*,vector2DList)

                print *, ' Following are the 2D vector fieldnames that will be read'
                do counter=1, num_of_vector2D_fields_to_read
                    WRITE(*,'(A25, A25)') trim(adjustl(list_vector2DX_fieldsNames(counter))), &
                                        trim(adjustl(list_vector2DY_fieldsNames(counter)))
                    self%list_vector2DX_fieldsNames(counter) = trim(adjustl(list_vector2DX_fieldsNames(counter)))
                    self%list_vector2DY_fieldsNames(counter) = trim(adjustl(list_vector2DY_fieldsNames(counter)))
                end do
                print *, ' '
                print *, ' '

                deallocate(list_vector2DX_fieldsNames) 
                deallocate(list_vector2DY_fieldsNames) 

            endif

            if ( num_of_vector3D_fields_to_read > 0) then
                
                allocate(list_vector3DX_fieldsNames(num_of_vector3D_fields_to_read))
                allocate(list_vector3DY_fieldsNames(num_of_vector3D_fields_to_read))
                allocate(list_vector3DZ_fieldsNames(num_of_vector3D_fields_to_read))

                allocate(self%list_vector3DX_fieldsNames(num_of_vector3D_fields_to_read))
                allocate(self%list_vector3DY_fieldsNames(num_of_vector3D_fields_to_read))
                allocate(self%list_vector3DZ_fieldsNames(num_of_vector3D_fields_to_read))

                read(*,vector3DList)

                print *, ' Following are the 3D vector fieldnames that will be read'
                do counter=1, num_of_vector3D_fields_to_read
                    WRITE(*,'(A25,A25, A25)') trim(adjustl(list_vector3DX_fieldsNames(counter))), &
                                        &     trim(adjustl(list_vector3DY_fieldsNames(counter))), &
                                        &     trim(adjustl(list_vector3DZ_fieldsNames(counter)))
                    self%list_vector3DX_fieldsNames(counter) = trim(adjustl(list_vector3DX_fieldsNames(counter)))
                    self%list_vector3DY_fieldsNames(counter) = trim(adjustl(list_vector3DY_fieldsNames(counter)))
                    self%list_vector3DZ_fieldsNames(counter) = trim(adjustl(list_vector3DZ_fieldsNames(counter)))
                end do
                print *, ' '
                print *, ' '

                deallocate(list_vector3DX_fieldsNames)
                deallocate(list_vector3DY_fieldsNames)
                deallocate(list_vector3DZ_fieldsNames)

            endif

            if ( nz > 0) then
                
                allocate(list_zlevels(nz), list_wzlevels(nz+1))
                allocate(self%list_zlevels(nz), self%list_wzlevels(nz+1))
                read(*,zlevelLists)

                print *, ' Following are the vertical level indices that will be read'
                do counter=1, nz
                    WRITE(*,'(I50)') list_zlevels(counter)
                    self%list_zlevels(counter) = list_zlevels(counter)
                end do
                print *, ' '
                print *, ' '

                print *, ' Following are the vertical level indices for top/bottom faces that will be read'
                do counter=1, nz+1
                    WRITE(*,'(I50)') list_wzlevels(counter)
                    self%list_wzlevels(counter) = list_wzlevels(counter)
                end do
                print *, ' '
                print *, ' '

                deallocate(list_zlevels, list_wzlevels)

            endif

            if ( ncoarse_levels > 0) then 
                allocate(list_coarse_factor_levels(ncoarse_levels))
                allocate(self%list_coarse_factor_levels(ncoarse_levels))
                read(*,coarseFactorList)
                print *, ' Following are the coarsening levels for Helmholtz decomposition ' 
                do counter=1, ncoarse_levels
                    WRITE(*,'(I50)') list_coarse_factor_levels(counter)
                    self%list_coarse_factor_levels(counter) = list_coarse_factor_levels(counter)
                end do
                print *, ' '
                print *, ' '
                deallocate(list_coarse_factor_levels)
            endif

            if ( nfilter > 0) then 
                allocate(list_filterLength(nfilter))
                allocate(self%list_filterLength(nfilter))
                read(*, filterLenList)
                print *, ' Following are the filterlengths ' 
                do counter=1, nfilter
                    WRITE(*,'(F50.3)') list_filterLength(counter)
                    self%list_filterLength(counter) = list_filterLength(counter)
                end do
                print *, ' '
                print *, ' '
                deallocate(list_filterLength)
            endif

            print *, "Configuration read and set SUCCESS"

        end subroutine new_configuration

        subroutine deallocate_vars(self)
            class(configuration), intent(inout) :: self
            if (allocated(self%list_filenames)) deallocate(self%list_filenames)
            if (allocated(self%list_scalar2D_fieldsNames)) deallocate(self%list_scalar2D_fieldsNames) 
            if (allocated(self%list_scalar3D_fieldsNames)) deallocate(self%list_scalar3D_fieldsNames)
            if (allocated(self%list_vector2DX_fieldsNames)) deallocate(self%list_vector2DX_fieldsNames)
            if (allocated(self%list_vector2DY_fieldsNames)) deallocate(self%list_vector2DY_fieldsNames)
            if (allocated(self%list_vector3DX_fieldsNames)) deallocate(self%list_vector3DX_fieldsNames)
            if (allocated(self%list_vector3DY_fieldsNames)) deallocate(self%list_vector3DY_fieldsNames)
            if (allocated(self%list_vector3DZ_fieldsNames)) deallocate(self%list_vector3DZ_fieldsNames)
            if (allocated(self%list_zlevels)) deallocate(self%list_zlevels)
            if (allocated(self%list_wzlevels)) deallocate(self%list_wzlevels)
            if (allocated(self%list_filterLength)) deallocate(self%list_filterLength)

        end subroutine

        subroutine init_config()
            ! Run the constructor for configuration object which also reads the configuration file
            if (taskid .EQ. MASTER) call config%construct()   
        end subroutine

end module
