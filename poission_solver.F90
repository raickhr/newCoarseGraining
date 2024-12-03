program poisson_solver
    use mkl
    implicit none

    integer, parameter :: n = 100  ! Size of the grid (n x n)
    integer :: i, j, n_eqns, ierr, iparm(64)
    double precision, dimension(n*n) :: A, b, phi
    integer, dimension(n*n) :: ia, ja
    integer :: mtype, maxfct, mnum, phase

    ! Initialize grid
    n_eqns = n * n

    ! Setup the sparse matrix (A), right-hand side vector (b)
    ! For example, assuming a simple 5-point stencil discretization
    call setup_poisson(A, ia, ja, b, n)

    ! Initialize PARDISO solver parameters
    mtype = 11  ! Real symmetric positive definite matrix
    maxfct = 1
    mnum = 1
    phase = 1  ! Initialize PARDISO
    iparm = 0
    iparm(1) = 1  ! Fill-in reordering method
    iparm(2) = 2  ! Pivoting control

    ! Call PARDISO to solve the system A * phi = b
    call pardiso(iparm, maxfct, mnum, mtype, phase, n_eqns, A, ia, ja, 0, 0, iparm, 0, phi, ierr)

    if (ierr .ne. 0) then
        print *, 'Error in PARDISO solver: ', ierr
        stop
    end if

    ! Print the solution
    print *, 'Solution phi:'
    do i = 1, n_eqns
        print *, phi(i)
    end do

end program poisson_solver

subroutine setup_poisson(A, ia, ja, b, n)
    implicit none
    integer, parameter :: nmax = 10000
    integer :: n
    double precision, dimension(nmax) :: A
    integer, dimension(nmax) :: ia, ja
    double precision, dimension(n) :: b
    integer :: i, j, k, nnz

    nnz = 0
    do i = 1, n
        do j = 1, n
            k = (i - 1) * n + j
            ! Set up the matrix entries for the 5-point stencil
            if (i > 1) then
                nnz = nnz + 1
                A(nnz) = -1.0d0
                ia(nnz) = k
                ja(nnz) = k - n
            end if
            if (i < n) then
                nnz = nnz + 1
                A(nnz) = -1.0d0
                ia(nnz) = k
                ja(nnz) = k + n
            end if
            if (j > 1) then
                nnz = nnz + 1
                A(nnz) = -1.0d0
                ia(nnz) = k
                ja(nnz) = k - 1
            end if
            if (j < n) then
                nnz = nnz + 1
                A(nnz) = -1.0d0
                ia(nnz) = k
                ja(nnz) = k + 1
            end if
            A(nnz+1) = 4.0d0
            ia(nnz+1) = k
            ja(nnz+1) = k

            b(k) = 1.0d0  ! Example right-hand side
        end do
    end do
end subroutine setup_poisson
