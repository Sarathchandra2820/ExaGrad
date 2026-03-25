FC = gfortran
FFLAGS = -g -fbacktrace -fbounds-check -O0 -Wall -Wextra -Isrc
LIBCINT = /Users/sarath/anaconda3/lib/python3.11/site-packages/pyscf/lib/deps/lib/libcint.6.dylib
RPATH   = -Wl,-rpath,/Users/sarath/anaconda3/lib/python3.11/site-packages/pyscf/lib/deps/lib
LIBS    = $(LIBCINT) $(RPATH) -llapack -lblas

# Compilation order matters: modules must be compiled before programs that use them
OBJS = src/libcint_interface.o src/integrals_gen.o src/rhf_main.o

all: rhf_main

rhf_main: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

src/%.o: src/%.f90
	$(FC) $(FFLAGS) -Jsrc -o $@ -c $<

clean:
	rm -f src/*.o src/*.mod rhf_main rhf_direct
