PROGRAM relax

  USE var_global
  USE grid

  IMPLICIT NONE

  CALL MPI_INIT(ierr) ! initiate MPI

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr) ! get # of processes

  CALL grid_setup

  CALL MPI_BARRIER(comm, ierr)

  CALL MPI_Finalize(ierr)

END PROGRAM relax
