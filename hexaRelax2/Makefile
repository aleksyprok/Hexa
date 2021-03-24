FFLAGS =
F90 = mpif90 -O3 --free-line-length-400 -fbounds-check

FFILES1 = relax_var.o relax_main.o

relax:$(FFILES1)
	$(F90) $(FFLAGS1) $(FFILES1) -o relax

relax_var.o:relax_var.f90
	$(F90) -c relax_var.f90 $(FFLAGS1)

relax_main.o:relax_main.f90 relax_var.o
	$(F90) -c  relax_main.f90 $(FFLAGS1)

clean:
	rm *% *.o *.mod core *~ *gtest* *out.dat* relax
