MODULE grid

  USE var_global

  IMPLICIT NONE

CONTAINS

  SUBROUTINE grid_setup

    !
    ! Determine sub array division of global arrays.
    ! USE MPI routine to split up  grid - dims gives number of processes in
    ! each di rection.

      call MPI_DIMS_CREATE(mpisize,mpidir,dims,ierr)
      nproc = dims
      print*,mpisize,mpidir,dims,ierr,nproc

    ! USE MPI routine to define cartesian grid.

      ! if (periodic .eq. 1) then
      !   periods(1)=.true.
      !   periods(2)=.true.
      !   print *, 'PERIODIC'
      ! else
      !   periods(1)=.false.
      !   periods(2)=.false.
      !   print *, 'CLOSED'
      ! endif


    ! sets of cartesian geometry with boundaries and places communicator in index comm.

      ! call MPI_CART_CREATE(MPI_COMM_WORLD,mpidir,nproc,periods,.TRUE.,comm,ierr)
      !
      ! call MPI_COMM_RANK(comm,rank,ierr)
      !
      ! CALL MPI_CART_COORDS(comm,rank,mpidir,coords,ierr)
      !
      ! call MPI_Get_processor_name(procname,namelen,ierr)

  END SUBROUTINE grid_setup

  ! SUBROUTINE setup_param
  !
  !     IF (rank .EQ. 0) THEN
  !
  !       print*, 'Read model parameters: ', setup_file
  !       OPEN(UNIT = 3, &
  !            FILE = setup_file, &
  !            FORM = 'FORMATTED', &
  !            STATUS = 'OLD')
  !       vsetup=get_value(3,'vsetup')    ! setup file version
  !       nmajor=get_value(3,'nmajor')
  !       nstrt  =get_value(3,'nstrt')
  !       nend  =get_value(3,'nend')
  !       etaia  =get_value(3,'etaia')
  !       eta4a  =get_value(3,'eta4a')
  !       periodic =get_value(3,'periodic')
  !       open=get_value(3,'open')
  !
  !       close(3)
  !
  !       print *,'nmajor=',nmajor
  !       print *,'nstrt=',nstrt
  !       print *,'nend =',nend
  !       print *,'etaia =',etaia
  !       print *,'eta4a =',eta4a
  !       print *,'periodic=',periodic
  !       print *, 'open =',open
  !
  !    endif
  !
  !   CALL MPI_BCAST(vsetup,1,MPI_INTEGER,0,comm,ierr)
  !   CALL MPI_BCAST(nmajor,1,MPI_INTEGER,0,comm,ierr)
  !   CALL MPI_BCAST(nminor,1,MPI_INTEGER,0,comm,ierr)
  !
  !   CALL MPI_BCAST(nstrt,1,MPI_INTEGER,0,comm,ierr)
  !   CALL MPI_BCAST(nend,1,MPI_INTEGER,0,comm,ierr)
  !   CALL MPI_BCAST(etaia,1,MPI_REAL,0,comm,ierr)
  !   CALL MPI_BCAST(eta4a,1,MPI_REAL,0,comm,ierr)
  !   CALL MPI_BCAST(periodic,1,MPI_INTEGER,0,comm,ierr)
  !   CALL MPI_BCAST(open,1,MPI_INTEGER,0,comm,ierr)
  !
  ! END SUBROUTINE setup_param

END MODULE grid
