PROG =	test_matrix

SRCS =	matrix_types.f90 test_matrix.f90 \
	base/base_hooks.f90 base/kinds.f90 base/machine.f90 \
	base/machine_internal.f90 test_utils.f90

OBJS =	matrix_types.o test_matrix.o \
	base/base_hooks.o base/kinds.o base/machine.o \
	base/machine_internal.o test_utils.o

LIBS = -framework Accelerate

CC = gcc
CFLAGS = -cpp -g 
FC = gfortran
FFLAGS =  -cpp -g
F90 = gfortran
F90FLAGS =  -cpp -g
LDFLAGS =

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

$(OBJS): %.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

test_matrix.o: matrix_types.o base/kinds.o test_utils.o
matrix_types.o: base/kinds.o base/machine.o base/base_uses.f90 base/base_hooks.o
negf_env_types.o: base/kinds.o matrix_types.o base/base_uses.f90 base/base_hooks.o
base/base_hooks.o: base/kinds.o base/machine.o
base/machine.o: base/kinds.o base/machine_internal.o
base/machine_internal.o: base/machine_posix.f90
