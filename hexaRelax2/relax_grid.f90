MODULE grid

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE grid_setup

    INTEGER :: namelen, num_hex_cells
    CHARACTER (LEN = MPI_MAX_PROCESSOR_NAME) :: procname

    CALL MPI_DIMS_CREATE(mpisize, mpidir, nproc, ierr) !nproc = dims

    IF (periodic .EQ. 1) THEN
      periods(1) = .TRUE.
      periods(2) = .TRUE.
    ELSE
      periods(1) = .FALSE.
      periods(2) = .FALSE.
    END IF

    PRINT*, MPI_COMM_WORLD, mpidir, nproc, periods, .TRUE., comm, ierr
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, mpidir, nproc, periods, .TRUE., comm, ierr)

    CALL MPI_COMM_RANK(comm, rank, ierr)

    CALL MPI_CART_COORDS(comm, rank, mpidir, coords, ierr)

    CALL MPI_GET_PROCESSOR_NAME(procname, namelen, ierr)
    ! PRINT*, rank, procname(1:namelen), coords(1), coords(2)

    CALL MPI_BARRIER(comm, ierr)

    CALL MPI_CART_SHIFT(comm, 0, 1, left, right, ierr)
    CALL MPI_CART_SHIFT(comm, 1, 1, down, up, ierr)
    ! PRINT*, 'Left/right', left, rank, right
    ! PRINT*, 'Up/down', down, rank, up

    ! Get rank of process at (0,0)
    CALL MPI_CART_RANK(comm, (/0, 0/), rankstart, ierr)

    IF (rank .EQ. rankstart) THEN
      OPEN(UNIT = 42, &
           FILE = parameters_file, &
           FORM = 'FORMATTED', &
           STATUS = 'OLD')
      READ(42, *) num_hex_cells
      READ(42, *) nxglobal, nyglobal, nzglobal
      READ(42, *) length_cm
      READ(42, *) time_s
      CLOSE(42)
    END IF
    CALL MPI_BCAST(num_hex_cells, 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(nxglobal, 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(nyglobal , 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(nzglobal, 1, MPI_INTEGER, rankstart, comm, ierr)
    CALL MPI_BCAST(length_cm, 1, MPI_real, rankstart, comm, ierr)
    CALL MPI_BCAST(time_s, 1, MPI_real, rankstart, comm, ierr)
                                                      !
    CALL MPI_BARRIER(comm,ierr)

    nz = nzglobal            ! no mpi in z direction
    nx = nxglobal / nproc(1) ! no of cells in x associated with each process
    ny = nyglobal / nproc(2) ! no of cells in y associated with each process

    IF ( ((nx * nproc(1)) .NE. nxglobal) .or.  &
         ((ny * nproc(2)) .NE. nyglobal) ) THEN
        IF ( (nx * nproc(1)) .NE. nxglobal) THEN
           PRINT*,'Unable to subdivide equally in x. Fix grid'
        ENDIF
        IF ( (ny * nproc(2)) .NE. nyglobal) THEN
           PRINT*,'Unable to subdivide equally in y. Fix grid'
        ENDIF
        CALL MPI_FINALIZE(ierr)
        STOP
    ENDIF

  END SUBROUTINE grid_setup

  SUBROUTINE setup_param

    ! This subroutine reads the setup file to check if periodic / open
    ! boundary conditions are used.

    ! Dummy strings to read files
    CHARACTER (LEN = 10) :: s1
    CHARACTER (LEN = 6 ) :: s2
    INTEGER :: i

    OPEN(UNIT   = 42, &
         FILE   = setup_file, &
         FORM   = 'FORMATTED', &
         STATUS = 'OLD', &
         ACTION = 'READ')

    ! Skip first 6 lines of file
    DO i = 1, 6
      READ(42, *)
    END DO

    ! Read periodic and open values
    READ(42, *) s1
    READ(42, *) s2

    ! Get integer values from the strings
    READ(s1(10:10), *) periodic
    READ(s2(6:6)  , *) open

    CLOSE(42)

  END SUBROUTINE setup_param

  SUBROUTINE arrayaloc

    ALLOCATE(aax(nx,   ny+1, nz+1))
    ALLOCATE(aay(nx+1, ny,   nz+1))
    ALLOCATE(aaz(nx+1, ny+1, nz  ))

  END SUBROUTINE arrayaloc

  SUBROUTINE arraydealoc

    DEALLOCATE(aax, aay, aaz)

  END SUBROUTINE arraydealoc


END MODULE grid
