#Fortran Compiler
FC = ifort
#FC = gfortran

#Compiler Options
#FCFLAGS = -O0 -traceback -g -check all -r8

PROGRAMS = main

all: $(PROGRAMS)

SRCS = $(wildcard *90)
f90s = $(SRCS:.F90=.f90)
OBJS = $(f90s:.f90=.o)

$(PROGRAMS): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJS): $(SRCS)


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

echo_all:
	echo $(SRCS)
	echo $(OBJS)
