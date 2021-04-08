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

  INTEGER :: periodic, open

  ! File and directory names
  CHARACTER (LEN = *), PARAMETER :: potential_field_file = 'run1/poten_00003p'
  CHARACTER (LEN = *), PARAMETER :: evolution_field_file = 'run1/run1_00003p'
  CHARACTER (LEN = *), PARAMETER :: parameters_file      = 'run1/param1'       ! Used to get nx, ny, nz
  CHARACTER (LEN = *), PARAMETER :: setup_file           = 'run1/run1_setup'   ! Used to check for periodic/open boundaries
  CHARACTER (LEN = *), PARAMETER :: output_dir           = 'run1'

  REAL, DIMENSION(:, :, :), ALLOCATABLE :: aax, aay, aaz
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: bbx, bby, bbz, bx, by, bz, bb
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: ccx, ccy, ccz, cx, cy, cz
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: eex, eey, eez, ex, ey, ez
  REAL, DIMENSION(:, :, :), ALLOCATABLE :: vx, vy, vz

  REAL, DIMENSION(:, :), ALLOCATABLE :: bbx0, bby0

  ! Variables associated with MPI
  INTEGER :: mpisize, ierr, rank, comm
  INTEGER :: left, right, up, down, nextrank, rankstart, rankend
  INTEGER, PARAMETER :: mpidir = 2, tag = 100
  INTEGER, DIMENSION(mpidir) :: nproc, coords
  integer, DIMENSION(MPI_STATUS_SIZE) :: stat
  LOGICAL, DIMENSION(mpidir) :: periods

  ! Variables associated with length and time
  REAL :: length_cm !length of one Hexa length unit in cm
  REAL :: time_s !time of one hexa time unit in seconds

END MODULE var_global
