#Fortran Compiler
#FC = ifort
FC = gfortran

#Compiler Options
#FCFLAGS = -O0 -traceback -g -check -warn
FCFLAGS = -O0 -Wall -g -fbacktrace -fbounds-check
#FCFLAGS = -O3

PROGRAMS = main

all: $(PROGRAMS)

SRCS = $(wildcard *90)
f90s = $(SRCS:.F90=.f90)
deps = $(f90s:.f90=.d)
OBJS = $(f90s:.f90=.o)

$(PROGRAMS): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

#$(OBJS): $(SRCS)

-include $(deps)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o %.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o %.mod: %.F90
	$(FC) $(FCFLAGS) -c $<

# Auto dependency generation

%.d: %.f90
	./Depend.pl $< . > $@

%.d: %.F90
	./Depend.pl $< . > $@

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *.d

veryclean: clean
	rm -f *~ $(PROGRAMS)
