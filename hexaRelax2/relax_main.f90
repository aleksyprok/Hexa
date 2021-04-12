PROGRAM relax

  USE var_global
  USE grid
  USE io
  USE cal

  IMPLICIT NONE

  timestep_s = 60. ! Timestep in seconds.

  CALL MPI_INIT(ierr) ! Initiate MPI

  CALL MPI_COMM_SIZE(comm, mpisize, ierr) ! Get # of processes

  CALL setup_param ! Check if periodic or open boundary conditions are used
  CALL grid_setup  ! Sets up computational grid and MPI (and reads in param1)

  basedt = timestep_s / time_s
  dt = MINVAL( (/ basedt, 0.1 /) ) ! Set dt to be basedt initially or 0.1 of the magnetogram cadence (whichever is smaller)

  CALL arrayaloc ! Allocates arrays

  CALL calc_boundary_field
  CALL calc_initial_field

  CALL writedata(0)

  CALL relax_routine

  CALL MPI_BARRIER(comm, ierr)

  CALL MPI_FINALIZE(ierr)

END PROGRAM relax
