COMMON_SRCS = base/base_hooks.f90 base/kinds.f90 base/machine.f90 \
	      base/machine_internal.f90 test_utils.f90
COMMON_OBJS = base/base_hooks.o base/kinds.o base/machine.o \
	      base/machine_internal.o test_utils.o

MATRIX_SRCS = matrix_types.f90 test_matrix.f90
MATRIX_OBJS = matrix_types.o test_matrix.o

SANCHO_SRCS = sancho_method.f90 test_sancho.f90 matrix_types.f90
SANCHO_OBJS = sancho_method.o test_sancho.o matrix_types.o

SRCS = $(MATRIX_SRCS) $(SANCHO_SRCS) $(COMMON_SRCS)
OBJS = $(MATRIX_OBJS) $(SANCHO_OBJS) $(COMMON_OBJS)

CC = gcc
CFLAGS = -cpp -g
FC = gfortran
FFLAGS =  -cpp -g
F90 = gfortran
F90FLAGS =  -cpp -g
LDFLAGS =
LIBS = -framework Accelerate

all: test_matrix test_sancho

test_matrix: $(MATRIX_OBJS) $(COMMON_OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(MATRIX_OBJS) $(COMMON_OBJS) $(LIBS)

test_sancho: $(SANCHO_OBJS) $(COMMON_OBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(SANCHO_OBJS) $(COMMON_OBJS) $(LIBS)

clean:
	rm -f test_sancho test_matrix $(OBJS) *.mod

%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

test_matrix.o: matrix_types.o base/kinds.o test_utils.o
test_sancho.o: sancho_method.o matrix_types.o base/kinds.o test_utils.o
sancho_method.o: matrix_types.o base/kinds.o
matrix_types.o: base/kinds.o base/machine.o base/base_uses.f90 base/base_hooks.o
base/base_hooks.o: base/kinds.o base/machine.o
base/machine.o: base/kinds.o base/machine_internal.o
base/machine_internal.o: base/machine_posix.f90
