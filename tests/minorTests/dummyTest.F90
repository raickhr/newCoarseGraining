program mpi_corrected_solver
    use petsc
    implicit none
  
    ! Declarations
    Mat :: A                ! System matrix
    Vec :: x, b             ! Solution and RHS vectors
    KSP :: ksp              ! Linear solver
    DM :: da                ! Distributed array
    PetscErrorCode :: ierr
    PetscMPIInt :: rank, size
    PetscInt :: mx, my, xs, ys, xm, ym, i, j, idx_phi, idx_psi, idx_div, idx_curl
    PetscReal :: hx, hy, val, rhs_val
  
    ! Domain and grid parameters
    mx = 10                ! Number of grid points in x
    my = 10                ! Number of grid points in y
    hx = 1.0 / (mx - 1)    ! Grid spacing in x
    hy = 1.0 / (my - 1)    ! Grid spacing in y
  
    ! Initialize PETSc and MPI
    call PetscInitialize(PETSC_NULL_CHARACTER, PETSC_NULL_INTEGER, PETSC_NULL_CHARACTER, 'MPI-enabled solver')
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
  
    ! Create distributed array (DMDA)
    call DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, &
         mx, my, PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, da, ierr)
    call DMSetUp(da, ierr)
  
    ! Get local grid ranges
    call DMDAGetCorners(da, xs, ys, PETSC_NULL_INTEGER, xm, ym, PETSC_NULL_INTEGER, ierr)
  
    ! Create global matrix and vectors
    call DMCreateMatrix(da, A, ierr)
    call DMCreateGlobalVector(da, b, ierr)
    call VecDuplicate(b, x, ierr)
  
    ! Assemble the system matrix A
    do j = ys, ys + ym - 1
       do i = xs, xs + xm - 1
          idx_phi = i + j * mx                ! Index for phi
          idx_psi = idx_phi + mx * my         ! Index for psi
          idx_div = idx_phi + 2 * mx * my     ! Index for -div V
          idx_curl = idx_phi + 3 * mx * my    ! Index for -curl V
  
          ! First row: phi equation
          call MatSetValue(A, idx_phi, idx_phi, -4.0 / (hx * hx) - 4.0 / (hy * hy), INSERT_VALUES, ierr)
          if (i > 0) call    MatSetValue(A, idx_phi, idx_phi-1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (i < mx-1) call MatSetValue(A, idx_phi, idx_phi+1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (j > 0) call    MatSetValue(A, idx_phi, idx_phi-mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
          if (j < my-1) call MatSetValue(A, idx_phi, idx_phi+mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
  
          ! Second row: psi equation
          call MatSetValue(A, idx_psi, idx_psi, -4.0 / (hx * hx) - 4.0 / (hy * hy), INSERT_VALUES, ierr)
          if (i > 0) call    MatSetValue(A, idx_psi, idx_psi-1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (i < mx-1) call MatSetValue(A, idx_psi, idx_psi+1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (j > 0) call    MatSetValue(A, idx_psi, idx_psi-mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
          if (j < my-1) call MatSetValue(A, idx_psi, idx_psi+mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
  
          ! Third row: -div V equation (Laplacian for phi)
          call MatSetValue(A, idx_div, idx_phi, -4.0 / (hx * hx) - 4.0 / (hy * hy), INSERT_VALUES, ierr)
          if (i > 0) call    MatSetValue(A, idx_div, idx_phi-1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (i < mx-1) call MatSetValue(A, idx_div, idx_phi+1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (j > 0) call    MatSetValue(A, idx_div, idx_phi-mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
          if (j < my-1) call MatSetValue(A, idx_div, idx_phi+mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
  
          ! Fourth row: -curl V equation (Laplacian for psi)
          call MatSetValue(A, idx_curl, idx_psi, -4.0 / (hx * hx) - 4.0 / (hy * hy), INSERT_VALUES, ierr)
          if (i > 0) call MatSetValue(A, idx_curl, idx_psi-1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (i < mx-1) call MatSetValue(A, idx_curl, idx_psi+1, 1.0 / (hx * hx), INSERT_VALUES, ierr)
          if (j > 0) call MatSetValue(A, idx_curl, idx_psi-mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
          if (j < my-1) call MatSetValue(A, idx_curl, idx_psi+mx, 1.0 / (hy * hy), INSERT_VALUES, ierr)
       end do
    end do
  
    ! Finalize matrix assembly
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
  
    ! Assemble the RHS vector b
    do j = ys, ys + ym - 1
       do i = xs, xs + xm - 1
          idx_phi = i + j * mx
          idx_psi = idx_phi + mx * my
          idx_div = idx_phi + 2 * mx * my
          idx_curl = idx_phi + 3 * mx * my
  
          ! Example RHS values: Replace with actual data
          call VecSetValue(b, idx_phi, 1.0, INSERT_VALUES, ierr)     ! RHS for phi
          call VecSetValue(b, idx_psi, 1.0, INSERT_VALUES, ierr)     ! RHS for psi
          call VecSetValue(b, idx_div, 0.0, INSERT_VALUES, ierr)     ! RHS for -div V
          call VecSetValue(b, idx_curl, 0.0, INSERT_VALUES, ierr)    ! RHS for -curl V
       end do
    end do
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)
  
    ! Create the KSP solver
    call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    call KSPSetFromOptions(ksp, ierr)
  
    ! Solve the linear system Ax = b
    call KSPSolve(ksp, b, x, ierr)
  
    ! Print solution (rank 0 only)
    if (rank == 0) then
       call PetscPrintf(PETSC_COMM_SELF, 'Solution vector:\n', ierr)
       call VecView(x, PETSC_VIEWER_STDOUT_SELF, ierr)
    end if
  
    ! Clean up
    call MatDestroy(A, ierr)
    call VecDestroy(b, ierr)
    call VecDestroy(x, ierr)
    call KSPDestroy
  