FC     = gfortran
FFLAGS = -O3 -march=native -fopenmp -Wall -ffree-line-length-none -Isrc
FFLAGS_DBG = -g -fbacktrace -fbounds-check -O0 -Wall -Wextra -fopenmp -ffree-line-length-none -Isrc
MODDIR = src

UNAME_S := $(shell uname -s)
PYTHON  ?= python3
PYSCF_LIBDIR := $(shell $(PYTHON) -c "import os, pyscf; print(os.path.join(os.path.dirname(pyscf.__file__), 'lib', 'deps', 'lib'))" 2>/dev/null)

ifeq ($(UNAME_S),Darwin)
LIBCINT = $(PYSCF_LIBDIR)/libcint.6.dylib
RPATH   = -Wl,-rpath,$(PYSCF_LIBDIR)
else ifeq ($(UNAME_S),Linux)
LIBCINT = /home/sarath/anaconda3/lib/python3.12/site-packages/pyscf/lib/deps/lib/libcint.so
RPATH   = -Wl,-rpath,/home/sarath/anaconda3/lib/python3.12/site-packages/pyscf/lib/deps/lib \
          -Wl,-rpath,/usr/lib/x86_64-linux-gnu/openblas-pthread
else
$(error Unsupported OS '$(UNAME_S)')
endif

ifeq ($(UNAME_S),Linux)
LIBS = $(LIBCINT) $(RPATH) -llapack \
       -L/usr/lib/x86_64-linux-gnu/openblas-pthread -lopenblas
else
LIBS = $(LIBCINT) $(RPATH) -llapack -lblas
endif

# ---------------------------------------------------------------------------
# Object files — compilation order encodes module dependency
# ---------------------------------------------------------------------------
OBJS = \
    src/interfaces/libcint_interface.o \
    src/interfaces/math_utils.o \
    src/types/molecule_t.o \
    src/integrals/one_electron_integrals.o \
    src/integrals/jk_contraction.o \
    src/integrals/two_electron_cholesky.o \
    src/integrals/two_electron_df.o \
    src/interfaces/transform_sigma.o \
    src/interfaces/davidson.o \
    src/scf/fock_builder.o \
    src/scf/scf_driver.o \
    src/properties/pol_initialisation.o \
    src/programs/rhf_main.o

all: rhf_main

rhf_main: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# ---------------------------------------------------------------------------
# Explicit module dependencies
# ---------------------------------------------------------------------------
src/interfaces/libcint_interface.o: src/interfaces/libcint_interface.f90
src/interfaces/math_utils.o:        src/interfaces/math_utils.f90
src/types/molecule_t.o:             src/types/molecule_t.f90 \
                                    src/interfaces/libcint_interface.o
src/integrals/one_electron_integrals.o: src/integrals/one_electron_integrals.f90 \
                                        src/interfaces/libcint_interface.o \
                                        src/interfaces/math_utils.o
src/integrals/jk_contraction.o:     src/integrals/jk_contraction.f90 \
                                    src/integrals/one_electron_integrals.o \
                                    src/interfaces/libcint_interface.o \
                                    src/interfaces/math_utils.o
src/integrals/two_electron_cholesky.o: src/integrals/two_electron_cholesky.f90 \
                                       src/integrals/one_electron_integrals.o \
                                       src/interfaces/libcint_interface.o \
                                       src/interfaces/math_utils.o
src/integrals/two_electron_df.o:    src/integrals/two_electron_df.f90 \
                                    src/integrals/one_electron_integrals.o \
                                    src/interfaces/libcint_interface.o \
                                    src/interfaces/math_utils.o
src/scf/fock_builder.o:             src/scf/fock_builder.f90 \
                                    src/integrals/jk_contraction.o \
                                    src/integrals/two_electron_cholesky.o \
                                    src/integrals/two_electron_df.o \
                                    src/integrals/one_electron_integrals.o \
                                    src/types/molecule_t.o
src/scf/scf_driver.o:               src/scf/scf_driver.f90 \
                                    src/scf/fock_builder.o \
                                    src/interfaces/math_utils.o
src/properties/cpks.o:              src/properties/cpks.f90 \
                                    src/integrals/jk_contraction.o \
                                    src/integrals/one_electron_integrals.o \
                                    src/interfaces/libcint_interface.o \
                                    src/interfaces/math_utils.o
src/programs/rhf_main.o:            src/programs/rhf_main.f90 \
                                    src/scf/fock_builder.o \
                                    src/scf/scf_driver.o \
                                    src/properties/cpks.o

# ---------------------------------------------------------------------------
# Pattern rule: compile .f90 -> .o, emit .mod files into src/
# ---------------------------------------------------------------------------
src/%.o: src/%.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

generate:
	python python/export_cint_env.py geometry/pyridine.xyz sto-6g
	python python/generate_pyscf_integrals.py

run: rhf_main generate
	./rhf_main

clean:
	rm -f src/*.o src/*.mod *.mod *.o rhf_main rhf_direct

debug: FFLAGS = $(FFLAGS_DBG)
debug: clean rhf_main
