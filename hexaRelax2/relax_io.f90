MODULE io

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE readdata(filename)

    CHARACTER (LEN = *), INTENT(IN) :: filename
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: aax_global, aay_global, aaz_global
    INTEGER, DIMENSION(mpidir) :: dumcord
    INTEGER :: i, j, k, opt

    ! Get aa_global and calculate aa for processor at (0, 0) cartesian coordinate
    IF (rank .EQ. rankstart) THEN

      ALLOCATE(aax_global(nxglobal,     nyglobal + 1, nzglobal + 1))
      ALLOCATE(aay_global(nxglobal + 1, nyglobal,     nzglobal + 1))
      ALLOCATE(aaz_global(nxglobal + 1, nyglobal + 1, nzglobal    ))

      PRINT*, 'Reading 3D model from ' // filename

      OPEN(UNIT = 42, &
           FILE = filename, &
           FORM = 'UNFORMATTED', &
           STATUS = 'OLD')

      READ(42) opt

      IF (.NOT. (opt .EQ. 1)) STOP 'Invalid option'

      READ(42) (((aax_global(i, j, k), i = 1, nxglobal    ), j = 1, nyglobal + 1), k = 1, nzglobal + 1)
      READ(42) (((aay_global(i, j, k), i = 1, nxglobal + 1), j = 1, nyglobal    ), k = 1, nzglobal + 1)
      READ(42) (((aaz_global(i, j, k), i = 1, nxglobal + 1), j = 1, nyglobal + 1), k = 1, nzglobal    )

      CLOSE(42)

      aax = aax_global(1:nx  , 1:ny+1, 1:nz+1)
      aay = aay_global(1:nx+1, 1:ny,   1:nz+1)
      aaz = aaz_global(1:nx+1, 1:ny+1, 1:nz  )

    END IF

    ! Get aa for the remaining processors
    IF(rank .EQ. rankstart) THEN

      DO j = 0, nproc(2) - 1
        DO i = 0, nproc(1) - 1
          IF ((i .EQ. 0) .AND. (j .EQ. 0)) CYCLE ! Already calculated aa for processor at (0, 0)
          dumcord(1) = i
          dumcord(2) = j
          CALL MPI_CART_RANK(comm, dumcord, nextrank, ierr)
          CALL MPI_SEND(aax_global((i * nx) + 1 : (i + 1) * nx, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   1 : nz + 1), &
                        nx * (ny + 1) * (nz + 1), &
                        MPI_REAL, nextrank, tag, comm, ierr)
          CALL MPI_SEND(aay_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny, &
                                   1 : nz + 1), &
                        (nx + 1) * ny * (nz + 1), &
                        MPI_REAL, nextrank, tag, comm, ierr)
          CALL MPI_SEND(aaz_global((i * nx) + 1 : (i + 1) * nx + 1, &
                                   (j * ny) + 1 : (j + 1) * ny + 1, &
                                   1 : nz), &
                        (nx + 1) * (ny + 1) * nz, &
                        MPI_REAL, nextrank, tag, comm, ierr)
        END DO
      END DO

    ELSE

      CALL MPI_RECV(aax, nx * (ny + 1) * (nz + 1), MPI_REAL, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(aay, (nx + 1) * ny * (nz + 1), MPI_REAL, rankstart, tag, comm, stat, ierr)
      CALL MPI_RECV(aaz, (nx + 1) * (ny + 1) * nz, MPI_REAL, rankstart, tag, comm, stat, ierr)

    END IF

  END SUBROUTINE readdata

  SUBROUTINE writedata(n)

    INTEGER, INTENT(IN) :: n
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: bbx_global, bby_global, bbz_global

    IF (rank .EQ. rankstart) THEN

       ALLOCATE(bbx_global(1:nxglobal+1, 0:nyglobal+1, 0:nzglobal+1))
       ALLOCATE(bby_global(0:nxglobal+1, 1:nyglobal+1, 0:nzglobal+1))
       ALLOCATE(bbz_global(0:nxglobal+1, 0:nyglobal+1, 1:nzglobal+1))

       aax_global(1:nx,1:ny+1,1:nz+1) = aax(1:nx,1:ny+1,1:nz+1)
       aay_global(1:nx+1,1:ny,1:nz+1) = aay(1:nx+1,1:ny,1:nz+1)
       aaz_global(1:nx+1,1:ny+1,1:nz)   = aaz(1:nx+1,1:ny+1,1:nz)
    ENDIF

  END SUBROUTINE writedata

END MODULE io
