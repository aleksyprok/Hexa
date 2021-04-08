MODULE cal

  USE var_global
  USE io

  IMPLICIT none

CONTAINS

  SUBROUTINE calc_boundary_field

    REAL, DIMENSION(1:nx+1, 0:ny+1) :: bbx1
    REAL, DIMENSION(0:nx+1, 1:ny+1) :: bby1

    CALL readdata(evolution_field_file)

    ! Calculate bbx0

    bbx1(:, 1:ny) =  &
        (aaz(:, 2:ny+1, 1) - aaz(:, 1:ny, 1)) / dely  &
      - (aay(:, :     , 2) - aay(:, :   , 1)) / delz

    ! Do vertical transfer
    call MPI_SENDRECV(bbx1(:, ny), nx + 1, MPI_REAL, up,   tag, &
                      bbx1(:, 0 ), nx + 1, MPI_REAL, down, tag, &
                      comm, stat, ierr)
    call MPI_SENDRECV(bbx1(:, 1   ), nx + 1, MPI_REAL, down, tag, &
                      bbx1(:, ny+1), nx + 1, MPI_REAL, up,   tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (down .EQ. MPI_PROC_NULL) bbx1(:, 0   ) = bbx1(:, 1 )
    IF (up   .EQ. MPI_PROC_NULL) bbx1(:, ny+1) = bbx1(:, ny)

    bbx0 = bbx1

    ! Calculate bby0

    bby1(1:nx, :) =  &
        (aax(:,      :, 2) - aax(:,    :, 1)) / delz  &
      - (aaz(2:nx+1, :, 1) - aaz(1:nx, :, 1)) / delx
    ! Do horizontal transfer
    call MPI_SENDRECV(bby1(nx, :), ny + 1, MPI_REAL, right, tag, &
                      bby1(0,  :), ny + 1, MPI_REAL, left,  tag, &
                      comm, stat, ierr)
    call MPI_SENDRECV(bby1(1   , :), ny + 1, MPI_REAL, left,  tag, &
                      bby1(nx+1, :), ny + 1, MPI_REAL, right, tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (left  .EQ. MPI_PROC_NULL) bby1(0,    :) = bby1(1,  :)
    IF (right .EQ. MPI_PROC_NULL) bby1(nx+1, :) = bby1(nx, :)

    bby0 = bby1

  END SUBROUTINE calc_boundary_field

  SUBROUTINE calc_initial_field

    CALL readdata(potential_field_file)

    ! Caclulate bbx

    bbx(:, 1:ny, 1:nz) =  &
        (aaz(:, 2:ny+1, :     ) - aaz(:, 1:ny, :   )) / dely  &
      - (aay(:, :     , 2:nz+1) - aay(:, :   , 1:nz)) / delz

    ! Do vertical transfer
    call MPI_SENDRECV(bbx(:, ny, 1:nz), (nx + 1) * nz, MPI_REAL, up,   tag, &
                      bbx(:, 0 , 1:nz), (nx + 1) * nz, MPI_REAL, down, tag, &
                      comm, stat, ierr)
    call MPI_SENDRECV(bbx(:, 1   , 1:nz), (nx + 1) * nz, MPI_REAL, down, tag, &
                      bbx(:, ny+1, 1:nz), (nx + 1) * nz, MPI_REAL, up,   tag, &
                      comm, stat, ierr)

  END SUBROUTINE calc_initial_field

  SUBROUTINE boundary_conditions

    IF (down .EQ. MPI_PROC_NULL) THEN
      bbx(:, 0, :) = bbx(:, 1, :)
    END IF

    IF (up .EQ. MPI_PROC_NULL) THEN
      bbx(:, ny+1, :) = bbx(:, ny, :)
    END IF

  END SUBROUTINE boundary_conditions

END MODULE cal
