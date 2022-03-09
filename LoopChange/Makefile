# Makefile

#use this for the GNU compiler
FC	= gfortran
FFLAGS 	= -O3 
CC	= gcc 

#use this for the Intel compiler
#FC=ifort
#FFLAGS= -O3 -opt-report-phase=hlo -nowarn -vec-report3
#CC=icc
#CFLAGS 	= -O3 -mssse3 -par-report -opt-report 3 -opt-report-phase=hlo
#icc -O1 -mssse3 -no-scalar-rep -opt-report 3 -opt-report-phase=hlo -restrict -c floatCycl.c -S -fcode-asm

#use this for the PGI compiler
#FC=pgf90
#FFLAGS= -O1 -Minfo
#CC=pgcc


PRG = loop

OBJ	= wallclock_time.o loop.o

$(PRG):	$(OBJ)
	$(FC) $(FFLAGS) -o $(PRG) $(OBJ) 

install: $(PRG)
	cp $(PRG) $(BIN)

clean:
		rm -f core $(OBJ) $(PRG)

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

