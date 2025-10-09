# Makefile

#use this for the GNU compiler
FC	= gfortran
FFLAGS 	= -O3 
CC	= gcc 

#use this for the Intel compiler
FC=ifx
FFLAGS= -O3 -qopt-report=3 -qopt-report-phase=all -axCORE-AVX2
CC=icx
CFLAGS = -O3 -qopt-report=3 -qopt-report-phase=all

#use this for the PGI compiler
#FC=pgf90
#FFLAGS= -O1 -Minfo
#CC=pgcc


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

