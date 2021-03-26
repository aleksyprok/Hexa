PROGRAM relax

  USE var_global
  USE grid

  IMPLICIT NONE

  CALL MPI_INIT(ierr) ! initiate MPI

  CALL MPI_COMM_SIZE(comm, mpisize, ierr) ! get # of processes

  CALL setup_param ! Check if periodic or open boundary conditions are used
  CALL grid_setup  ! Sets up computational grid and MPI (and reads in param1)

  PRINT*, rank, nx, ny, nz

  CALL MPI_BARRIER(comm, ierr)

  CALL MPI_FINALIZE(ierr)

END PROGRAM relax
