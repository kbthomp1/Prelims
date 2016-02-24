FC = gfortran

#FCFLAGS = -fbacktrace -fcheck=all -g -fdefault-real-8 -fdefault-double-8 
FCFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 
#FCFLAGS += -I/usr/include

PROGRAMS = main

all: $(PROGRAMS)

main.o: grid.o
main: grid.o
main.o: flowsolver.o
main: flowsolver.o
main.o: output.o
main: output.o

# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

