.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

default: code

include Makefile.mac

OBJECTS = \
module_structures.o\
itaux.o\
iters.o\
spmv_kernels.o\
output_for_paraview.o\
mammoth.o

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90
.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f
code:	$(OBJECTS)
	$(F90) $(OPTIONS) $(OBJECTS) $(LIBS) -o mammoth
clean:	
	rm -f *.o
	rm -f *.mod
	rm -f *.dat 
	rm -f *.out 
	rm -f opla
	rm -f fort.*
	rm -f mammoth


