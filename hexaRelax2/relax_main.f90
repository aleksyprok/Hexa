PROGRAM relax

  USE var_global
  USE grid
  USE io

  IMPLICIT NONE

  CALL MPI_INIT(ierr) ! Initiate MPI

  CALL MPI_COMM_SIZE(comm, mpisize, ierr) ! Get # of processes

  CALL setup_param ! Check if periodic or open boundary conditions are used
  CALL grid_setup  ! Sets up computational grid and MPI (and reads in param1)

  CALL arrayaloc ! Allocates arrays

  CALL readdata(evolution_field_file)
  CALL readdata(potential_field_file) ! Read in initial A

  CALL MPI_BARRIER(comm, ierr)

  CALL MPI_FINALIZE(ierr)

END PROGRAM relax
