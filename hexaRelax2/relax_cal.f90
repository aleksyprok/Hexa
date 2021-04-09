
MODULE cal

  USE var_global
  USE io

  IMPLICIT none

CONTAINS


  SUBROUTINE calc_boundary_field

    REAL, DIMENSION(1:nx+1, 0:ny+1) :: bbx1
    REAL, DIMENSION(0:nx+1, 1:ny+1) :: bby1
    REAL, DIMENSION(0:nx+1, 0:ny+1) :: bbz1

    CALL readdata(evolution_field_file)

    ! Calculate bbx1

    bbx1(:, 1:ny) =  &
        (aaz(:, 2:ny+1, 1) - aaz(:, 1:ny, 1)) / dely  &
      - (aay(:, :     , 2) - aay(:, :   , 1)) / delz

    ! Vertical transfer
    CALL MPI_SENDRECV(bbx1(:, ny), nx + 1, MPI_REAL, up,   tag, &
                      bbx1(:, 0 ), nx + 1, MPI_REAL, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbx1(:, 1   ), nx + 1, MPI_REAL, down, tag, &
                      bbx1(:, ny+1), nx + 1, MPI_REAL, up,   tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (down .EQ. MPI_PROC_NULL) bbx1(:, 0   ) = bbx1(:, 1 )
    IF (up   .EQ. MPI_PROC_NULL) bbx1(:, ny+1) = bbx1(:, ny)

    ! Calculate bby1

    bby1(1:nx, :) =  &
        (aax(:,      :, 2) - aax(:,    :, 1)) / delz  &
      - (aaz(2:nx+1, :, 1) - aaz(1:nx, :, 1)) / delx
    ! Horizontal transfer
    CALL MPI_SENDRECV(bby1(nx, :), ny + 1, MPI_REAL, right, tag, &
                      bby1(0,  :), ny + 1, MPI_REAL, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bby1(1   , :), ny + 1, MPI_REAL, left,  tag, &
                      bby1(nx+1, :), ny + 1, MPI_REAL, right, tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (left  .EQ. MPI_PROC_NULL) bby1(0,    :) = bby1(1,  :)
    IF (right .EQ. MPI_PROC_NULL) bby1(nx+1, :) = bby1(nx, :)


    ! Calculate bbz1

    bbz1(1:nx, 1:ny) =  &
        (aay(2:nx+1,      :, 1) - aay(1:nx,    :, 1)) / delx  &
      - (aax(     :, 2:ny+1, 1) - aax(   :, 1:ny, 1)) / dely
    ! Horizontal transfer
    CALL MPI_SENDRECV(bbz1(nx, :), ny + 2, MPI_REAL, right, tag, &
                      bbz1(0,  :), ny + 2, MPI_REAL, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz1(1   , :), ny + 2, MPI_REAL, left,  tag, &
                      bbz1(nx+1, :), ny + 2, MPI_REAL, right, tag, &
                      comm, stat, ierr)
    ! Vertical transfer
    CALL MPI_SENDRECV(bbz1(:, ny), nx + 2, MPI_REAL, up,   tag, &
                      bbz1(:, 0 ), nx + 2, MPI_REAL, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz1(:, 1   ), nx + 2, MPI_REAL, down, tag, &
                      bbz1(:, ny+1), nx + 2, MPI_REAL, up,   tag, &
                      comm, stat, ierr)
    ! Apply boundary conditions
    IF (left  .EQ. MPI_PROC_NULL) bbz1(0,    :) = bbz1(1,  :)
    IF (right .EQ. MPI_PROC_NULL) bbz1(nx+1, :) = bbz1(nx, :)
    IF (down  .EQ. MPI_PROC_NULL) bbz1(:, 0   ) = bbz1(:, 1 )
    IF (up    .EQ. MPI_PROC_NULL) bbz1(:, ny+1) = bbz1(:, ny)

    ! Calculate bbx0, bby0
    bbx0 = bbx1 - delz / delx * (bbz1(1:nx+1, :) - bbz1(0:nx, :))
    bby0 = bby1 - delz / dely * (bbz1(:, 1:ny+1) - bbz1(:, 0:ny))
    ! bbx0 = bbx1
    ! bby0 = bby1

  END SUBROUTINE calc_boundary_field

  SUBROUTINE calc_initial_field

    CALL readdata(potential_field_file)

    bbx(:, 1:ny, 1:nz) =  &
        (aaz(:, 2:ny+1, :     ) - aaz(:, 1:ny, :   )) / dely  &
      - (aay(:, :     , 2:nz+1) - aay(:, :   , 1:nz)) / delz

    bby(1:nx, :, 1:nz) =  &
        (aax(:     , :, 2:nz+1) - aax(:   , :, 1:nz)) / delz  &
      - (aaz(2:nx+1, :, :     ) - aaz(1:nx, :, :   )) / delx

    bbz(1:nx, 1:ny, :) =  &
        (aay(2:nx+1, :     , :) - aay(1:nx, :   , :)) / delx  &
      - (aax(:     , 2:ny+1, :) - aax(:   , 1:ny, :)) / dely

    CALL horizontal_transfer
    CALL vertical_transfer
    CALL boundary_conditions

  END SUBROUTINE calc_initial_field

  SUBROUTINE boundary_conditions

    IF (left .EQ. MPI_PROC_NULL) THEN
      bby(0, :, :) = bby(1, :, :)
      bbz(0, :, :) = bbz(1, :, :)
    END IF

    IF (right .EQ. MPI_PROC_NULL) THEN
      bby(nx+1, :, :) = bby(nx, :, :)
      bbz(nx+1, :, :) = bbz(nx, :, :)
    END IF

    IF (down .EQ. MPI_PROC_NULL) THEN
      bbx(:, 0, :) = bbx(:, 1, :)
      bbz(:, 0, :) = bbz(:, 1, :)
    END IF

    IF (up .EQ. MPI_PROC_NULL) THEN
      bbx(:, ny+1, :) = bbx(:, ny, :)
      bbz(:, ny+1, :) = bbz(:, ny, :)
    END IF

    bbx(:,:,0) = bbx0
    bby(:,:,0) = bby0

    bbx(:,:,nz+1) = bbx(:,:,nz)
    bby(:,:,nz+1) = bby(:,:,nz)

  END SUBROUTINE boundary_conditions

  SUBROUTINE horizontal_transfer

    CALL MPI_SENDRECV(bby(nx, :, :), (ny + 1) * (nz + 2), MPI_REAL, right, tag, &
                      bby(0,  :, :), (ny + 1) * (nz + 2), MPI_REAL, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bby(1   , :, :), (ny + 1) * (nz + 2), MPI_REAL, left,  tag, &
                      bby(nx+1, :, :), (ny + 1) * (nz + 2), MPI_REAL, right, tag, &
                      comm, stat, ierr)

    CALL MPI_SENDRECV(bbz(nx, :, :), (ny + 2) * (nz + 1), MPI_REAL, right, tag, &
                      bbz(0,  :, :), (ny + 2) * (nz + 1), MPI_REAL, left,  tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz(1   , :, :), (ny + 2) * (nz + 1), MPI_REAL, left,  tag, &
                      bbz(nx+1, :, :), (ny + 2) * (nz + 1), MPI_REAL, right, tag, &
                      comm, stat, ierr)

  END SUBROUTINE horizontal_transfer

  SUBROUTINE vertical_transfer

    CALL MPI_SENDRECV(bbx(:, ny, :), (nx + 1) * (nz + 2), MPI_REAL, up,   tag, &
                      bbx(:, 0 , :), (nx + 1) * (nz + 2), MPI_REAL, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbx(:, 1   , :), (nx + 1) * (nz + 2), MPI_REAL, down, tag, &
                      bbx(:, ny+1, :), (nx + 1) * (nz + 2), MPI_REAL, up,   tag, &
                      comm, stat, ierr)

    CALL MPI_SENDRECV(bbz(:, ny, :), (nx + 2) * (nz + 1), MPI_REAL, up,   tag, &
                      bbz(:, 0 , :), (nx + 2) * (nz + 1), MPI_REAL, down, tag, &
                      comm, stat, ierr)
    CALL MPI_SENDRECV(bbz(:, 1   , :), (nx + 2) * (nz + 1), MPI_REAL, down, tag, &
                      bbz(:, ny+1, :), (nx + 2) * (nz + 1), MPI_REAL, up,   tag, &
                      comm, stat, ierr)

  END SUBROUTINE vertical_transfer

END MODULE cal
