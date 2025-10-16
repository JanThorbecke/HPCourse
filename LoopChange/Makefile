# Makefile

#use this for the GNU compiler
FC	= gfortran
FFLAGS 	= -O1 
CC	= gcc 
CFLAGS = -O3

#use this for the Intel compiler
#FC=ifx
#FFLAGS= -O0 -qopt-report=3 -qopt-report-phase=all -axCORE-AVX2
#CC=icx
#CFLAGS = -O3 -qopt-report=3 -qopt-report-phase=all

#use this for the Nvdia compiler
FC=nvfortran
FFLAGS= -O1 -Minfo=all
CC=nvc
CFLAGS = -O3  -Manno -Minfo=all

# For AMD aocc compilera
FC=flang
FFLAGS= -O3 -fgen-aor
CC=clang
CFLAGS = -O3 -fgen-aor

PRG = loop

OBJ	= wallclock_time.o loop.o

$(PRG):	$(OBJ)
	$(FC) $(FFLAGS) -o $(PRG) $(OBJ) 

loop2c: wallclock_time.o loop2c.o
	$(CC) $(CFLAGS) -o loop2c  $^

install: $(PRG)
	cp $(PRG) $(BIN)

clean:
		rm -f core $(OBJ) $(PRG) loop*.o

print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar

.SUFFIXES: .o .c .f90

.f90.o: 
	$(FC) $(FFLAGS) -c $<

.c.o: .c
	$(CC) $(CFLAGS) -c $<

