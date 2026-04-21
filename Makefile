FC     = gfortran
FFLAGS = -O3 -march=native -fopenmp -Wall -ffree-line-length-none -Isrc
FFLAGS_DBG = -g -fbacktrace -fbounds-check -O0 -Wall -Wextra -fopenmp -ffree-line-length-none -Isrc
MODDIR = src

UNAME_S := $(shell uname -s)
PYTHON  ?= python3
CONDA   ?= $(shell command -v conda 2>/dev/null)
CONDA_ENV ?= exagrad
PYSCF_LIBDIR := $(shell $(PYTHON) -c "import os, pyscf; print(os.path.join(os.path.dirname(pyscf.__file__), 'lib', 'deps', 'lib'))" 2>/dev/null)

ifeq ($(strip $(PYSCF_LIBDIR)),)
ifneq ($(strip $(CONDA)),)
PYSCF_LIBDIR := $(shell $(CONDA) run -n $(CONDA_ENV) python -c "import os, pyscf; print(os.path.join(os.path.dirname(pyscf.__file__), 'lib', 'deps', 'lib'))" 2>/dev/null)
endif
endif

ifeq ($(strip $(PYSCF_LIBDIR)),)
$(error Could not locate PySCF libcint library. Activate a PySCF-enabled environment or run 'make CONDA_ENV=<env-name>')
endif

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
ROSE_DIR = src/interfaces/rose_interface

COMMON_OBJS = \
    src/interfaces/libcint_interface.o \
    src/interfaces/math_utils.o \
    src/types/molecule_t.o \
    src/types/molecule_loader.o \
    src/integrals/one_electron_integrals.o \
    src/integrals/jk_contraction.o \
    src/integrals/two_electron_cholesky.o \
    src/integrals/two_electron_df.o \
    src/interfaces/transform_sigma.o \
    src/interfaces/davidson.o \
    src/scf/fock_builder.o \
    src/scf/scf_driver.o \
    src/properties/pol_initialisation.o \
    src/properties/cpks.o

ROSE_LIB_OBJS = \
    $(ROSE_DIR)/rose_global.o \
    $(ROSE_DIR)/linear_algebra.o \
    $(ROSE_DIR)/rose_utils_stub.o \
    $(ROSE_DIR)/rose_check_stub.o \
    $(ROSE_DIR)/rose_properties_stub.o \
    $(ROSE_DIR)/make_iao.o \
    $(ROSE_DIR)/make_ibo.o \
    $(ROSE_DIR)/rose_interface.o \
    $(ROSE_DIR)/rose_davidson.o

OBJS = $(COMMON_OBJS) src/programs/rhf_main.o
ROSE_OBJS = $(COMMON_OBJS) $(ROSE_LIB_OBJS) src/interfaces/utils.o src/programs/main_rose.o

all: rhf_main

rhf_main: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

rhf_rose_main: $(ROSE_OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# ---------------------------------------------------------------------------
# Explicit module dependencies
# ---------------------------------------------------------------------------
src/interfaces/libcint_interface.o: src/interfaces/libcint_interface.f90
src/interfaces/math_utils.o:        src/interfaces/math_utils.f90
src/interfaces/utils.o:             src/interfaces/utils.f90
$(ROSE_DIR)/rose_global.o:          $(ROSE_DIR)/rose_global.F90
$(ROSE_DIR)/linear_algebra.o:       $(ROSE_DIR)/linear_algebra.F90 \
                                    $(ROSE_DIR)/rose_global.o
$(ROSE_DIR)/rose_utils_stub.o:      $(ROSE_DIR)/rose_utils_stub.F90 \
                                    $(ROSE_DIR)/rose_global.o
$(ROSE_DIR)/rose_check_stub.o:      $(ROSE_DIR)/rose_check_stub.F90
$(ROSE_DIR)/rose_properties_stub.o: $(ROSE_DIR)/rose_properties_stub.F90 \
                                    $(ROSE_DIR)/rose_global.o
$(ROSE_DIR)/make_iao.o:             $(ROSE_DIR)/make_iao.F90 \
                                    $(ROSE_DIR)/linear_algebra.o \
                                    $(ROSE_DIR)/rose_utils_stub.o
$(ROSE_DIR)/make_ibo.o:             $(ROSE_DIR)/make_ibo.F90 \
                                    $(ROSE_DIR)/linear_algebra.o \
                                    $(ROSE_DIR)/rose_utils_stub.o \
                                    $(ROSE_DIR)/rose_check_stub.o \
                                    $(ROSE_DIR)/rose_properties_stub.o
$(ROSE_DIR)/rose_interface.o:       $(ROSE_DIR)/rose_interface.f90 \
                                    $(ROSE_DIR)/make_iao.o \
                                    $(ROSE_DIR)/make_ibo.o
$(ROSE_DIR)/rose_davidson.o:        $(ROSE_DIR)/rose_davidson.f90 \
                                    src/types/molecule_t.o
src/types/molecule_t.o:             src/types/molecule_t.f90 \
                                    src/interfaces/libcint_interface.o
src/types/molecule_loader.o:        src/types/molecule_loader.f90 \
                                    src/types/molecule_t.o
src/integrals/one_electron_integrals.o: src/integrals/one_electron_integrals.f90 \
                                        src/types/molecule_loader.o \
                                        src/interfaces/libcint_interface.o \
                                        src/interfaces/math_utils.o
src/integrals/jk_contraction.o:     src/integrals/jk_contraction.f90 \
                                                src/types/molecule_loader.o \
                                    src/interfaces/libcint_interface.o \
                                    src/interfaces/math_utils.o
src/integrals/two_electron_cholesky.o: src/integrals/two_electron_cholesky.f90 \
                                                    src/types/molecule_loader.o \
                                       src/interfaces/libcint_interface.o \
                                       src/interfaces/math_utils.o
src/integrals/two_electron_df.o:    src/integrals/two_electron_df.f90 \
                                                src/types/molecule_loader.o \
                                    src/interfaces/libcint_interface.o \
                                    src/interfaces/math_utils.o
src/interfaces/transform_sigma.o:   src/interfaces/transform_sigma.f90 \
                                    src/interfaces/math_utils.o \
                                    src/integrals/jk_contraction.o \
                                    src/interfaces/libcint_interface.o \
                                    src/integrals/two_electron_df.o
src/scf/fock_builder.o:             src/scf/fock_builder.f90 \
                                    src/integrals/jk_contraction.o \
                                    src/integrals/two_electron_cholesky.o \
                                    src/integrals/two_electron_df.o \
                                                src/types/molecule_loader.o \
                                    src/types/molecule_t.o
src/scf/scf_driver.o:               src/scf/scf_driver.f90 \
                                    src/scf/fock_builder.o \
                                    src/interfaces/math_utils.o
src/properties/pol_initialisation.o: src/properties/pol_initialisation.f90 \
                                    src/interfaces/transform_sigma.o \
                                    src/integrals/jk_contraction.o \
                                    src/integrals/two_electron_df.o \
                                    src/interfaces/math_utils.o
src/properties/cpks.o:              src/properties/cpks.f90 \
                                    src/interfaces/davidson.o \
                                    src/interfaces/transform_sigma.o \
                                    src/integrals/two_electron_df.o \
                                    src/integrals/two_electron_cholesky.o \
                                    src/integrals/jk_contraction.o \
                                    src/integrals/one_electron_integrals.o \
                                    src/interfaces/math_utils.o
src/programs/rhf_main.o:            src/programs/rhf_main.f90 \
                                    src/scf/fock_builder.o \
                                    src/scf/scf_driver.o \
                                    $(ROSE_DIR)/rose_interface.o \
                                    src/properties/pol_initialisation.o \
                                    src/properties/cpks.o
src/programs/main_rose.o:           src/programs/main_rose.f90 \
                                    src/types/molecule_loader.o \
                                    src/interfaces/utils.o \
                                    src/scf/fock_builder.o \
                                    src/scf/scf_driver.o \
                                    src/properties/pol_initialisation.o \
                                    src/properties/cpks.o \
                                    $(ROSE_DIR)/rose_interface.o \
                                    $(ROSE_DIR)/rose_davidson.o

# ---------------------------------------------------------------------------
# Pattern rules: compile .f90 / .F90 -> .o, emit .mod files into src/
# ---------------------------------------------------------------------------
src/%.o: src/%.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

src/%.o: src/%.F90
	$(FC) $(FFLAGS) -J$(MODDIR) -o $@ -c $<

generate:
	python python/export_cint_env.py geometry/pyridine.xyz sto-6g
	python python/generate_pyscf_integrals.py

run: rhf_main generate
	./rhf_main

clean:
	find src -type f \( -name '*.o' -o -name '*.mod' \) -delete
	rm -f *.mod *.o rhf_main rhf_rose_main rhf_direct

debug: FFLAGS = $(FFLAGS_DBG)
debug: clean rhf_main
