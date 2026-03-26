FC = gfortran
MODDIR = src
INCLUDES = -I$(MODDIR) -Isrc/types -Isrc/interfaces -Isrc/integrals -Isrc/scf -Isrc/properties
FFLAGS = -O3 -march=native -fopenmp -Wall -ffree-line-length-none $(INCLUDES)
FFLAGS_DBG = -g -fbacktrace -fbounds-check -O0 -Wall -Wextra -fopenmp -ffree-line-length-none $(INCLUDES)
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
LIBCINT = /Users/sarath/anaconda3/lib/python3.11/site-packages/pyscf/lib/deps/lib/libcint.6.dylib
RPATH   = -Wl,-rpath,/Users/sarath/anaconda3/lib/python3.11/site-packages/pyscf/lib/deps/lib
else ifeq ($(UNAME_S),Linux)
LIBCINT = /home/sarath/anaconda3/lib/python3.12/site-packages/pyscf/lib/deps/lib/libcint.so
RPATH   = -Wl,-rpath,/home/sarath/anaconda3/lib/python3.12/site-packages/pyscf/lib/deps/lib \
          -Wl,-rpath,/usr/lib/x86_64-linux-gnu/openblas-pthread
else
$(error Unsupported OS '$(UNAME_S)')
endif

# Linux links against threaded OpenBLAS; macOS uses Accelerate-backed BLAS
ifeq ($(UNAME_S),Linux)
LIBS    = $(LIBCINT) $(RPATH) -llapack \
          -L/usr/lib/x86_64-linux-gnu/openblas-pthread -lopenblas
else
LIBS    = $(LIBCINT) $(RPATH) -llapack -lblas
endif

# Compilation order matters: modules must be compiled before programs that use them
OBJS = \
	src/types/molecule_t.o \
	src/interfaces/libcint_interface.o \
	src/interfaces/math_utils.o \
	src/integrals/one_electron_integrals.o \
	src/scf/fock_builder.o \
	src/scf/scf_driver.o \
	src/properties/cpks.o \
	src/programs/rhf_main.o

all: rhf_main

rhf_main: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# Explicit dependencies for modules
src/types/molecule_t.o: src/types/molecule_t.f90
src/interfaces/libcint_interface.o: src/interfaces/libcint_interface.f90
src/interfaces/math_utils.o: src/interfaces/math_utils.f90
src/integrals/one_electron_integrals.o: src/integrals/one_electron_integrals.f90 src/types/molecule_t.o src/interfaces/libcint_interface.o src/interfaces/math_utils.o
src/scf/fock_builder.o: src/scf/fock_builder.f90 src/integrals/one_electron_integrals.o src/interfaces/libcint_interface.o src/interfaces/math_utils.o
src/scf/scf_driver.o: src/scf/scf_driver.f90 src/integrals/one_electron_integrals.o src/interfaces/math_utils.o src/scf/fock_builder.o
src/properties/cpks.o: src/properties/cpks.f90 src/interfaces/libcint_interface.o src/integrals/one_electron_integrals.o src/interfaces/math_utils.o
src/programs/rhf_main.o: src/programs/rhf_main.f90 src/integrals/one_electron_integrals.o src/scf/scf_driver.o src/properties/cpks.o

rhf_direct: src/types/molecule_t.o src/interfaces/libcint_interface.o src/interfaces/math_utils.o src/integrals/one_electron_integrals.o src/programs/rhf_direct.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

src/programs/rhf_direct.o: src/programs/rhf_direct.f90 src/types/molecule_t.o src/integrals/one_electron_integrals.o

src/%.o: src/%.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

generate:
	python python/export_cint_env.py geometry/pyridine.xyz sto-6g
	python python/generate_pyscf_integrals.py

run: rhf_main generate
	./rhf_main

clean:
	find src -type f \( -name '*.o' -o -name '*.mod' \) -delete
	rm -f *.mod *.o rhf_main rhf_direct

# Debug build: make debug
debug: FFLAGS = $(FFLAGS_DBG)
debug: clean rhf_main
