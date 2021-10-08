# Makefile

#Intel compiler which supports openMP
FC	= ifort
FFLAGS 	= -O3 -openmp -static

#GNU compiler which supports openMP
FC	= gfortran
FFLAGS 	= -O3 -fopenmp

all: test_omp

test_omp: test_omp.f
	$(FC) $(FFLAGS) -o $@ $@.f $(LIBS) -lm

.f.o: .f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f core test_omp

print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar
