ATI_DISPLAY    = DISPLAY=:0.0
HMPP           = hmpp

SRCS =	initial.f90 sorties.f90
OBJS =	initial.o sorties.o

FC = $(HMPP) gfortran
FC = gfortran
OPTS = -O4 -p -W
FFLAGS = $(OPTS) 

all : vpp.exe 

vpp.exe: $(OBJS) main.F90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

run: $(PROG)
	$(ATI_DISPLAY) ./$(PROG)

.SUFFIXES: $(SUFFIXES) .o .f90 .mod .f

.f90.o:
	$(FC) $(FFLAGS) -c $<
.mod.o:
	$(FC) $(FFLAGS) -c $*.f90

main.o: main.f90 initial.o sorties.o
initial.o: initial.f90
sorties.o: sorties.f90

clean:
	rm -f $(PROG) *.o *.mod fort.* data/* *.dx *.gnu
	rm -f *.exe *.so *.so.* *.cu *.cu.* *.linkinfo *.hmd
