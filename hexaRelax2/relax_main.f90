PROGRAM relax

  USE var_global
  USE grid
  USE io
  USE cal

  IMPLICIT NONE

  CALL MPI_INIT(ierr) ! Initiate MPI

  CALL MPI_COMM_SIZE(comm, mpisize, ierr) ! Get # of processes

  CALL setup_param ! Check if periodic or open boundary conditions are used
  CALL grid_setup  ! Sets up computational grid and MPI (and reads in param1)

  CALL arrayaloc ! Allocates arrays

  CALL calc_boundary_field

  CALL MPI_BARRIER(comm, ierr)

  CALL MPI_FINALIZE(ierr)

END PROGRAM relax
