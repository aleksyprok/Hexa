MODULE grid

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE grid_setup

      CALL MPI_DIMS_CREATE(mpisize, mpidir, dims, ierr)
      nproc = dims
      PRINT*, mpisize, mpidir, dims, ierr, nproc

      IF (periodic .EQ. 1) THEN
        periods(1) = .TRUE.
        periods(2) = .TRUE.
        PRINT*, 'PERIODIC'
      ELSE
        periods(1) = .FALSE.
        periods(2) = .FALSE.
        PRINT*, 'CLOSED'
      END IF



    ! sets of cartesian geometry with boundaries and places communicator in index comm.

      ! call MPI_CART_CREATE(MPI_COMM_WORLD,mpidir,nproc,periods,.TRUE.,comm,ierr)
      !
      ! call MPI_COMM_RANK(comm,rank,ierr)
      !
      ! CALL MPI_CART_COORDS(comm,rank,mpidir,coords,ierr)
      !
      ! call MPI_Get_processor_name(procname,namelen,ierr)

  END SUBROUTINE grid_setup

  SUBROUTINE setup_param

    ! This subroutine reads the setup file to check if periodic / open
    ! boundary conditions are used.

    ! Dummy strings (used to read file)
    CHARACTER (LEN = 10) :: s1
    CHARACTER (LEN = 6 ) :: s2
    INTEGER :: i

    IF (rank .EQ. 0) THEN

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

      CALL MPI_BCAST(periodic, 1, MPI_INTEGER, 0, comm,ierr)
      CALL MPI_BCAST(open    , 1, MPI_INTEGER, 0, comm,ierr)

    END IF

  END SUBROUTINE setup_param


END MODULE grid
