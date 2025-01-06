module multiGridHelmHoltz
    use kinds   
    use configurationMod
    use coarsening
    use interpolation
    use mpiwrapper
    use operators
    use helmHoltzDecomp
    use fields
    use gridModule
    use read_write
    implicit none

    type :: grid
        real(kind=real_kind), allocatable, dimension(:,:) :: centerDx, centerDy, & ! dx and dy array
                                                             cellArea, &  ! cell area
                                                             lat, lon, & ! coordinates
                                                             uvel, vvel, &  ! vectors for Helmholtz decomp
                                                             psi, phi ! toroidal and poloidal scalars
        integer :: nx, ny, coarseLevel
    contains
            procedure :: setGrid => set_grid     ! method to set grid value at grid level
            procedure :: delGrid => del_grid
            procedure :: setUvelVvelFromOrig => set_uvel_vvelFromOrig
            procedure :: setUvelVvelByVal => set_uvel_vvel_byVal
            procedure :: resetUvelVvel => reset_uvel_vvel
    end type

    integer :: ncoarse_factors, maximum_iterations
    real(kind=real_kind) :: absolute_tolerance, relative_tolerance, divergence_tolerance
    integer, allocatable :: factorList(:)

    class(grid), allocatable :: multiGrid(:)
    
    contains

        subroutine init_helmholtz()
            call initPETSC()

            if (taskid == MASTER) then
                call set_ncoarse_factors(config%ncoarse_levels)
                call set_tolerances(config%abs_tol, config%rel_tol, config%div_tol)
                call set_maximum_iterations(config%max_iterations)
            endif

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call MPI_BCAST(ncoarse_factors, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(absolute_tolerance, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(relative_tolerance, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(divergence_tolerance, 1, MPI_REAL , MASTER, MPI_COMM_WORLD, i_err)
            call MPI_BCAST(maximum_iterations, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)
    
            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call MPI_BCAST(ncoarse_factors, 1, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

            allocate(factorList(ncoarse_factors))
            
            if (taskid == MASTER) call set_coarse_factorlist(config%list_coarse_factor_levels)

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            call MPI_BCAST(factorList, ncoarse_factors, MPI_INTEGER , MASTER, MPI_COMM_WORLD, i_err)

            call MPI_Barrier(MPI_COMM_WORLD, i_err)

            allocate(multiGrid(ncoarse_factors))
            allocate(multiGridMats(ncoarse_factors))

            call setMultiGrid()

        end subroutine

        subroutine set_ncoarse_factors(n)
            integer(kind=int_kind), intent(in) :: n
            ncoarse_factors = n
        end subroutine

        subroutine set_coarse_factorlist(list_coarse_factor_levels)
            integer, intent(in) :: list_coarse_factor_levels(:)
            factorList = list_coarse_factor_levels
        end subroutine

        subroutine set_maximum_iterations(n)
            integer(kind=int_kind), intent(in) :: n
            maximum_iterations = n
        end subroutine

        subroutine set_tolerances(abs_tol, rel_tol, div_tol)
            real(kind=real_kind), intent(in) :: abs_tol, rel_tol, div_tol
            absolute_tolerance = abs_tol
            relative_tolerance = rel_tol
            divergence_tolerance = div_tol
        end subroutine

        subroutine set_grid(self, coarseLevel, org_lat, org_lon, org_centerDx, org_centerDy, org_cellArea)
            class(grid), intent(inout) :: self
            integer, intent(in) :: coarseLevel
            real, intent(in) :: org_lat(:,:), org_lon(:,:), org_centerDx(:,:), org_centerDy(:,:), org_cellArea(:,:)
            
            integer :: nx, ny, shapeArr(2)

            shapeArr = shape(org_lat)
            nx = shapeArr(1)
            ny = shapeArr(2)

            if (coarseLevel > 1) then 
                call coarsenLatLon(nx, ny, coarseLevel, org_lat, org_lon, self%lat, self%lon)
                call coarsenDXDY(nx, ny, coarseLevel, org_centerDx, org_centerDy, self%centerDx, self%centerDy, downCenterUp = 0)
                call coarsenAREA(nx, ny, coarseLevel, org_cellArea, self%cellArea)
            else if (coarseLevel == 1) then
                allocate(self%lat(nx, ny), self%lon(nx, ny), self%cellArea(nx, ny), &
                         self%centerDx(nx, ny), self%centerDy(nx, ny))
                self%lat = org_lat
                self%lon = org_lon
                self%centerDx = org_centerDx
                self%centerDy = org_centerDy
                self%cellArea = org_cellArea
            end if

            shapeArr = shape(self%lat)
            self%nx = shapeArr(1)
            self%ny = shapeArr(2)
            self%coarseLevel = coarseLevel

        end subroutine

        subroutine del_grid(self)
            class(grid), intent(inout) :: self
            deallocate(self%lat, self%lon)
            deallocate(self%centerDx, self%centerDy)
            deallocate(self%cellArea)
            deallocate(self%uvel, self%vvel, self%psi, self%phi) 
        end subroutine

        subroutine set_uvel_vvelFromOrig(self, factor, org_uvel, org_vvel, org_cellArea)
            class(grid), intent(inout) :: self
            real, intent(in) :: org_uvel(:,:), org_vvel(:,:), org_cellArea(:,:)
            
            integer :: factor, nx, ny, shapeArr(2)

            if (taskid == 0) then
                shapeArr = shape(org_uvel)
                nx = shapeArr(1)
                ny = shapeArr(2)

                if (factor > 1) then
                    call coarsenField(nx, ny, factor, org_uvel, org_cellArea, self%uvel)
                    call coarsenField(nx, ny, factor, org_vvel, org_cellArea, self%vvel)
                    
                else if (factor == 1 ) then 
                    if (.not. allocated(self%uvel)) allocate(self%uvel(nx,ny), self%vvel(nx, ny))
                    self%uvel = org_uvel
                    self%vvel = org_vvel
                end if

                shapeArr = shape(self%uvel)
                nx = shapeArr(1)
                ny = shapeArr(2)
                if (nx .NE. self%nx .or. ny .NE. self%ny) stop &
                'Error in coarsening UVEl for HelmHoltz Decomp. First Set Multi Grid'


                shapeArr = shape(self%vvel)
                nx = shapeArr(1)
                ny = shapeArr(2)
                if (nx .NE. self%nx .or. ny .NE. self%ny) stop &
                'Error in coarsening VVEl for HelmHoltz Decomp. First Set Multi Grid'
            endif
        end subroutine

        subroutine set_uvel_vvel_byVal(self, uvel, vvel)
            class(grid), intent(inout) :: self
            real, intent(in) :: uvel(:,:), vvel(:,:)

            integer :: factor, nx, ny, shapeArr(2)

            if (taskid == 0 ) then
                shapeArr = shape(uvel)
                nx = shapeArr(1)
                ny = shapeArr(2)

                self%uvel = uvel
                self%vvel = vvel

                shapeArr = shape(self%uvel)
                nx = shapeArr(1)
                ny = shapeArr(2)
                if (nx .NE. self%nx .or. ny .NE. self%ny) stop &
                'Error in coarsening UVEl for HelmHoltz Decomp. First Set Multi Grid Or Check Multigrid'


                shapeArr = shape(self%vvel)
                nx = shapeArr(1)
                ny = shapeArr(2)
                if (nx .NE. self%nx .or. ny .NE. self%ny) stop &
                'Error in coarsening VVEl for HelmHoltz Decomp. First Set Multi Grid Or Check Multigrid'
            endif

        end subroutine

        subroutine reset_uvel_vvel(self)
            class(grid) :: self
            if (taskid == 0) deallocate(self%uvel, self%vvel)
        end subroutine

        subroutine setMultiGrid()
            integer :: i
            real(kind=real_kind), allocatable, dimension(:,:) :: dummy
            
            do i = 1, ncoarse_factors
                call multiGrid(i)%setGrid(factorList(i), ULAT, ULONG, DXU, DYU, UAREA)
                call MPI_Barrier(MPI_COMM_WORLD, i_err)

                call multiGridMats(i)%setMat(multiGrid(i)%nx, multiGrid(i)%ny, &
                                            multiGrid(i)%centerDx, multiGrid(i)%centerDy)
                call MPI_Barrier(MPI_COMM_WORLD, i_err)

                allocate(dummy(multiGrid(i)%nx, multiGrid(i)%ny))
                dummy = 0.0

                call multiGridMats(i)%setRHS(multiGrid(i)%nx, multiGrid(i)%ny, &
                                            multiGrid(i)%centerDx, multiGrid(i)%centerDy, dummy, dummy)
                call MPI_Barrier(MPI_COMM_WORLD, i_err)

                call multiGridMats(i)%setLHS(multiGrid(i)%nx, multiGrid(i)%ny, dummy, dummy)
                deallocate(dummy)
            end do
        end subroutine

        subroutine delMultiGrid()
            deallocate(multiGrid, factorList)
            call finalizePETSC()
        end subroutine

        subroutine solveByMultiGrid(nx, ny, uvel, vvel, cellarea, phi, psi, &
                                    max_Iter, rel_Tol, abs_Tol, div_Tol)
            integer , intent(in) :: nx, ny
            real(kind=real_kind), intent(in) :: uvel(nx, ny), vvel(nx, ny), cellArea(nx, ny)
            real(kind=real_kind), intent(out) :: phi(nx, ny), psi(nx, ny)

            integer, intent(in) :: max_Iter
            real, intent(in) :: rel_Tol, abs_Tol, div_Tol
            integer :: i, factor, nfactors

            real(kind=real_kind), allocatable, dimension(:,:) :: crs_phi, crs_psi,  &
                                                                 wrk_phi, wrk_psi, &
                                                                 phi_seed, psi_seed, &
                                                                 res_uvel, res_vvel, &
                                                                 uvel_pol, vvel_pol, &
                                                                 uvel_tor, vvel_tor



            nfactors = size(factorList)

            do i = 1, nfactors
                factor = factorList(i)
                if (rank == 0) print *, 'Setting UVEL VVEL at coarse scale', factor
                call multiGrid(i)%setUvelVvelFromOrig(factor, uvel, vvel, cellArea)
            end do

            call MPI_Barrier(MPI_COMM_WORLD, i_err)
            if (taskid == 0) print *, 'Starting Multigrid'

            do i = 1, nfactors 
                factor = factorList(i)
                if (taskid == 0 ) then
                    print *, ''
                    print *, ''
                    print *, 'helmholtz at coarse level', factor  
                    allocate(wrk_phi(multiGrid(i)%nx, multiGrid(i)%ny), wrk_psi(multiGrid(i)%nx, multiGrid(i)%ny), &
                             phi_seed(multiGrid(i)%nx, multiGrid(i)%ny), psi_seed(multiGrid(i)%nx, multiGrid(i)%ny), &
                             res_uvel(multiGrid(i)%nx, multiGrid(i)%ny), res_vvel(multiGrid(i)%nx, multiGrid(i)%ny), &
                             uvel_pol(multiGrid(i)%nx, multiGrid(i)%ny), vvel_pol(multiGrid(i)%nx, multiGrid(i)%ny), &
                             uvel_tor(multiGrid(i)%nx, multiGrid(i)%ny), vvel_tor(multiGrid(i)%nx, multiGrid(i)%ny))

                    if (i == 1) then
                        wrk_phi = 0
                        wrk_psi = 0
                    else
                        call blinearInterpolationLatLon(multiGrid(i-1)%lat(1,:), multiGrid(i-1)%lon(:,1), &   ! previous grid
                                                    multiGrid(i)%lat(1,:), multiGrid(i)%lon(:,1), &       ! current grid
                                                    crs_psi, wrk_psi)

                        call blinearInterpolationLatLon(multiGrid(i-1)%lat(1,:), multiGrid(i-1)%lon(:,1), &   ! previous grid
                                                        multiGrid(i)%lat(1,:), multiGrid(i)%lon(:,1), &       ! current grid
                                                        crs_phi, wrk_phi)
                        
                        deallocate(crs_phi, crs_psi)
                    endif
                endif

                call MPI_Barrier(MPI_COMM_WORLD, i_err)
                
                call multiGridMats(i)%resetRHS(multiGrid(i)%nx, multiGrid(i)%ny, &
                                               multiGrid(i)%centerDx, multiGrid(i)%centerDy, &
                                               multiGrid(i)%uvel, multiGrid(i)%vvel)

                call multiGridMats(i)%resetLHS(multiGrid(i)%nx, multiGrid(i)%ny, wrk_phi, wrk_psi)

                call multiGridMats(i)%solve(maxIter = max_Iter, relTol = rel_Tol, absTol = abs_Tol, divTol = div_Tol)

                if (taskid == 0) allocate(crs_phi(multiGrid(i)%nx, multiGrid(i)%ny), &
                                          crs_psi(multiGrid(i)%nx, multiGrid(i)%ny))

                call multiGridMats(i)%getSol(multiGrid(i)%nx, multiGrid(i)%ny, crs_phi, crs_psi)

                if (taskid == 0 ) then                    
                    phi_seed = crs_phi
                    psi_seed = crs_psi
                    
                    call getPolTorVelFD(crs_psi, crs_phi, multiGrid(i)%centerDx, multiGrid(i)%centerDy, &
                                        uvel_pol, uvel_tor, vvel_pol, vvel_tor)

                    res_uvel = multiGrid(i)%uvel - uvel_pol - uvel_tor 
                    res_vvel = multiGrid(i)%vvel - vvel_pol - vvel_tor
                    wrk_phi = 0.0 
                    wrk_psi = 0.0

                    print *, 'residual calculation complete'

                endif
                call MPI_Barrier(MPI_COMM_WORLD, i_err)

                call multiGridMats(i)%resetRHS(multiGrid(i)%nx, multiGrid(i)%ny, &
                                               multiGrid(i)%centerDx, multiGrid(i)%centerDy,&
                                               res_uvel, res_vvel)

                call multiGridMats(i)%resetLHS(multiGrid(i)%nx, multiGrid(i)%ny, wrk_phi, wrk_psi)

                if (rank == 0 ) print *, 'LHS RHS from residual set'

                call MPI_Barrier(MPI_COMM_WORLD, i_err)

                call multiGridMats(i)%solve(maxIter = max_Iter, relTol = rel_Tol, &
                                            absTol = abs_Tol, divTol = div_Tol)
                call multiGridMats(i)%getSol(multiGrid(i)%nx, multiGrid(i)%ny, crs_phi, crs_psi)
                
                if (rank == 0) then 
                    crs_phi = crs_phi + phi_seed
                    crs_psi = crs_psi + psi_seed
                    if (i == nfactors) then
                        phi = crs_phi
                        psi = crs_psi 
                        print *, 'collected solutions in original grid'
                    endif
                end if

                call multiGridMats(i)%resetLHS(multiGrid(i)%nx, multiGrid(i)%ny, crs_phi, crs_psi)

                call MPI_Barrier(MPI_COMM_WORLD, i_err)

                if (rank == 0 ) then 
                    deallocate( wrk_phi, wrk_psi, &
                                phi_seed, psi_seed, &
                                res_uvel, res_vvel, &
                                uvel_pol, vvel_pol, &
                                uvel_tor, vvel_tor)

                    if ( i == nfactors) deallocate(crs_phi, crs_psi)
                    print *, 'Loop completed'
                    print *, ''
                    print *, ''
                endif

                call MPI_Barrier(MPI_COMM_WORLD, i_err)


            end do
            

        end subroutine


        subroutine helmholtzDecompAllVecFields()
            integer :: counter, z_counter, nx, ny, nz
            integer :: max_Iter
            real :: rel_Tol, abs_Tol, div_Tol

            max_Iter = maximum_iterations
            rel_Tol = relative_tolerance
            abs_Tol = absolute_tolerance
            div_Tol = divergence_tolerance

            nx = nxu
            ny = nyu
            nz = nzu

            do counter=1, num_vector2D_fields
                if (taskid == 0 ) then
                    print *, 'Starting Helmholtz Decomposition for 2D vector field number', counter
                end if

                call solveByMultiGrid(nx, ny, &
                                    vector2DX_fields(:, :, counter), &
                                    vector2DY_fields(:, :, counter), &
                                    UAREA, &
                                    phi2D_fields(:,:, counter), &
                                    psi2D_fields(:,:, counter), &
                                    max_Iter, rel_Tol, abs_Tol, div_Tol)

                if (taskid == 0) then
                    call getPolTorVelFD(psi2D_fields(:,:, counter), phi2D_fields(:,:, counter), DXU, DYU,  &
                                        vector2DX_phi_fields(:, :, counter), vector2DX_psi_fields(:, :, counter), &
                                        vector2DY_phi_fields(:, :, counter), vector2DY_psi_fields(:, :, counter))
                    call write2dVar('psi2D.nc', 'psi2D', psi2D_fields(:,:, counter))
                    call write2dVar('phi2D.nc', 'phi2D', phi2D_fields(:,:, counter))
                endif
            end do


            do counter=1, num_vector3D_fields
                do z_counter =1, nz
                    if (taskid == 0 ) then
                        print *, 'Starting Helmholtz Decomposition for 3D vector field number', counter, &
                                 ' at z count ', z_counter
                    end if

                    call solveByMultiGrid(nx, ny, &
                                        vector3DX_fields(:, :, z_counter, counter), &
                                        vector3DY_fields(:, :, z_counter, counter), &
                                        UAREA, &
                                        phi3D_fields(:,:,z_counter, counter), &
                                        psi3D_fields(:,:,z_counter, counter), &
                                        max_Iter, rel_Tol, abs_Tol, div_Tol)

                    if (taskid == 0) then
                        call getPolTorVelFD(psi3D_fields(:,:,z_counter, counter), &
                                            phi3D_fields(:,:,z_counter, counter), DXU, DYU, &
                                            vector3DX_phi_fields(:, :, z_counter, counter), &
                                            vector3DX_psi_fields(:, :, z_counter, counter), &
                                            vector3DY_phi_fields(:, :, z_counter, counter), &
                                            vector3DY_psi_fields(:, :, z_counter, counter))
                    endif
                end do
            end do

        end subroutine

        subroutine getVectorsFromFilteredHelmHoltzPotentials(num_filterlengths)
            integer, intent(in) :: num_filterlengths
            integer :: ell_counter, counter, z_counter, nx, ny, nz
            integer :: max_Iter
            real :: rel_Tol, abs_Tol, div_Tol

            max_Iter = maximum_iterations
            rel_Tol = relative_tolerance
            abs_Tol = absolute_tolerance
            div_Tol = divergence_tolerance

            nx = nxu
            ny = nyu
            nz = nzu

            do ell_counter =1, num_filterlengths
                do counter=1, num_vector2D_fields
                    if (taskid == 0) then
                        call getPolTorVelFD(OL_psi2D_fields(:,:, counter, ell_counter),&
                                            OL_phi2D_fields(:,:, counter, ell_counter),&
                                            DXU, DYU,  &
                                            OL_vector2DX_phi_fields(:, :, counter, ell_counter), &
                                            OL_vector2DX_psi_fields(:, :, counter, ell_counter), &
                                            OL_vector2DY_phi_fields(:, :, counter, ell_counter), &
                                            OL_vector2DY_psi_fields(:, :, counter, ell_counter))
                    endif
                end do


                do counter=1, num_vector3D_fields
                    do z_counter =1, nz
                        if (taskid == 0) then
                            call getPolTorVelFD(OL_psi3D_fields(:,:,z_counter, counter, ell_counter), &
                                                OL_phi3D_fields(:,:,z_counter, counter, ell_counter), &
                                                DXU, DYU, &
                                                OL_vector3DX_phi_fields(:, :, z_counter, counter, ell_counter), &
                                                OL_vector3DX_psi_fields(:, :, z_counter, counter, ell_counter), &
                                                OL_vector3DY_phi_fields(:, :, z_counter, counter, ell_counter), &
                                                OL_vector3DY_psi_fields(:, :, z_counter, counter, ell_counter))
                        endif
                    end do
                end do
            end do
        end subroutine

        

    
    

end module
