PROGRAM relax

  USE var_global

  IMPLICIT none

  CALL MPI_INIT(ierr) ! initiate MPI

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, ierr) ! get # of processes

  PRINT*, mpisize

  CALL MPI_BARRIER(comm, ierr)

  CALL MPI_Finalize(ierr)

END PROGRAM relax
