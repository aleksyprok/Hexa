MODULE cal

  USE var_global
  USE io

  IMPLICIT none

CONTAINS

  SUBROUTINE calc_boundary_field

    REAL, DIMENSION(1:nx+1, 0:ny+1) :: bbx1
    REAL, DIMENSION(0:nx+1, 1:ny+1) :: bby1

    CALL readdata(evolution_field_file)

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

  END SUBROUTINE calc_boundary_field

  SUBROUTINE calc_initial_field

    ! ! Caclulate bbx
    !
    ! ! Caclulate curl of A
    ! bbx(:, 1:ny, 1:nz) =  &
    !     (aaz(:, 2:ny+1, 1:nz  ) - aaz(:, 1:ny, 1:nz)) / dely  &
    !   - (aay(:, 1:ny  , 2:nz+1) - aay(:, 1:ny, 1:nz)) / delz
    !
    !
    ! ! Impose boundary conditions
    ! IF (down .EQ.MPI_PROC_NULL ) THEN
    !   bbx(:, 0, 2:nz+1)=bbx(:,   2,2:nz+1)
    ! ENDIF
    ! IF (up .EQ. MPI_PROC_NULL) THEN
    !   bbx(1:nx+1,ny+2,2:nz+1)=bbx(1:nx+1,ny+1,2:nz+1)
    ! ENDIF
    !
    ! bbx(1:nx+1,1:ny+2,   1)=bbx(1:nx+1,1:ny+2,   2)  &
    !      -delz/delx*(bbz(2:nx+2,1:ny+2,1)-bbz(1:nx+1,1:ny+2,1))

  END SUBROUTINE calc_initial_field

END MODULE cal
