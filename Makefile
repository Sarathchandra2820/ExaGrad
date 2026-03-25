FC = gfortran
FFLAGS = -O3 -march=native -fopenmp -Wall -Isrc
FFLAGS_DBG = -g -fbacktrace -fbounds-check -O0 -Wall -Wextra -fopenmp -Isrc
LIBCINT = /Users/sarath/anaconda3/lib/python3.11/site-packages/pyscf/lib/deps/lib/libcint.6.dylib
RPATH   = -Wl,-rpath,/Users/sarath/anaconda3/lib/python3.11/site-packages/pyscf/lib/deps/lib
LIBS    = $(LIBCINT) $(RPATH) -llapack -lblas

# Compilation order matters: modules must be compiled before programs that use them
OBJS = src/libcint_interface.o src/math_utils.o src/integrals_gen.o src/scf.o src/rhf_main.o

all: rhf_main

rhf_main: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# Explicit dependencies for modules
src/integrals_gen.o: src/integrals_gen.f90 src/libcint_interface.o src/math_utils.o
src/math_utils.o:    src/math_utils.f90
src/scf.o:           src/scf.f90 src/integrals_gen.o src/math_utils.o
src/rhf_main.o:      src/rhf_main.f90 src/integrals_gen.o src/scf.o

src/%.o: src/%.f90
	$(FC) $(FFLAGS) -Jsrc -o $@ -c $<

generate:
	python python/export_cint_env.py geometry/pyridine.xyz sto-6g
	python python/generate_pyscf_integrals.py

run: rhf_main generate
	./rhf_main

clean:
	rm -f src/*.o src/*.mod *.mod *.o rhf_main rhf_direct

# Debug build: make debug
debug: FFLAGS = $(FFLAGS_DBG)
debug: clean rhf_main
