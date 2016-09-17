# -*- mode: makefile -*-
SRCS =	matrix_types.f90 test_matrix.f90 \
	base/base_hooks.f90 base/kinds.f90 base/machine.f90 \
	base/machine_internal.f90 test_utils.f90 \
        integral.f90 negf_env_types.f90 sancho_method.f90
	common/mathconstants.f90

OBJS_test_matrix =  matrix_types.o test_matrix.o \
                    base/base_hooks.o base/kinds.o base/machine.o \
                    base/machine_internal.o test_utils.o

OBJS_test_quad = integral.o negf_env_types.o sancho_method.o \
		 matrix_types.o test_quad.o\
		 base/base_hooks.o base/kinds.o base/machine.o \
		 base/machine_internal.o test_utils.o common/mathconstants.o

LIBS = -framework Accelerate
# LIBS = /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf

#ln -s /usr/lib/liblapack.so.3gf  /usr/lib/liblapack.so
#ln -s /usr/lib/libblas.so.3gf  /usr/lib/libblas.so
#LDFLAGS = -L/usr/lib -Wl,-rpath=/usr/lib
#LIBS = -lblas llapack

CC = gcc
CFLAGS = -cpp -g 
FC = gfortran
FFLAGS =  -cpp -g
F90 = gfortran
F90FLAGS =  -cpp -g
LDFLAGS =

all: test_matrix test_quad

test_matrix: $(OBJS_test_matrix)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJS_test_matrix) $(LIBS)

test_quad: $(OBJS_test_quad)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJS_test_quad) $(LIBS)

clean:
	rm -f test_matrix test_quad *.o *.mod

$(OBJS_test_matrix): %.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(OBJS_test_quad): %.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

test_matrix.o: matrix_types.o base/kinds.o test_utils.o
matrix_types.o: base/kinds.o base/machine.o base/base_uses.f90 base/base_hooks.o
negf_env_types.o: base/kinds.o matrix_types.o base/base_uses.f90 base/base_hooks.o
                  common/mathconstants.o
base/base_hooks.o: base/kinds.o base/machine.o
base/machine.o: base/kinds.o base/machine_internal.o
base/machine_internal.o: base/machine_posix.f90
test_quad.o: integral.o negf_env_types.o sancho_method.o \
             matrix_types.o base/kinds.o test_utils.o \
             common/mathconstants.o
integral.o: negf_env_types.o matrix_types.o base/kinds.o \
            common/mathconstants.o
sancho_method.o: matrix_types.o base/kinds.o
