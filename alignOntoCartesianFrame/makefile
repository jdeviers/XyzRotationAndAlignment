PROG = $(wildcard prog_*.f90)
MODS = $(wildcard mod_*.f90)
OBJS = $(patsubst %.f90,%.o,$(MODS))

FC = gfortran
FCFLAGS = -Wall -Wextra -fbacktrace -fcheck=all

PROGRAM = orient_xyz

#
default : $(PROGRAM)

$(PROGRAM) : $(OBJS)
	$(FC) -o $@ $(PROG) $^ $(FCFLAGS)

$(OBJS) : %.o : %.f90
	$(FC) -c $< -o $@ $(FCFLAGS)

mod_procedures.o mod_rotations.o mod_translation.o : mod_precision.o
mod_alignment.o : mod_precision.o mod_procedures.o mod_rotations.o

# 
debug :
	@echo $(PROGRAM)
	@echo $(PROG)
	@echo $(MODS)
	@echo $(OBJS)
	
#
clean :
	rm $(PROGRAM) $(OBJS) $(patsubst %.o,%.mod,$(OBJS))

.PHONY = default debug clean
