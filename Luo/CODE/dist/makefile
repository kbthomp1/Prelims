#Fortran Compiler
#FC = ifort
FC = gfortran

#Compiler Options
#FCFLAGS = -O0 -traceback -g -check -warn
#FCFLAGS = -O0 -Wall -g -fbacktrace -fbounds-check
FCFLAGS = -O3

PROGRAMS = main

all: $(PROGRAMS)

SRCS = $(wildcard *90)
f90s = $(SRCS:.F90=.f90)
deps = $(f90s:.f90=.d)
OBJS = $(f90s:.f90=.o)

$(PROGRAMS): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

#$(OBJS): $(SRCS)

-include .make.deps

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o %.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o %.mod: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
