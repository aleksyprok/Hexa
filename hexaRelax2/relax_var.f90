MODULE var_global

  IMPLICIT NONE

  INCLUDE "mpif.h"

  REAL, PARAMETER :: pi = 3.141592654
  REAL, PARAMETER :: frc_coef = 3000. ! Frictional coefficient (km^2/s)

  ! Grid parameters:
  REAL :: delx, dely, delz

  ! Number of cells (not corners):
  INTEGER :: nxglobal, nyglobal, nzglobal
  INTEGER :: nx, ny, nz

  REAL, DIMENSION(:, :, :), ALLOCATABLE :: bbx, bby, bbz, bx, by, bz, bb
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: ccx, ccy, ccz, cx, cy, cz
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: eex, eey, eez, ex, ey, ez
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: vx, vy, vz

  ! variables associated with MPI
  INTEGER :: mpisize, ierr, rank, comm
  INTEGER :: left, right, up, down, nextrank, rankstart, rankend

END MODULE var_global
