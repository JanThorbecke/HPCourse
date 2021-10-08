# Makefile

CC	= icc
CFLAGS 	= -O3 -vec-report
#icc -O1 -mssse3 -no-scalar-rep -opt-report 3 -opt-report-phase=hlo -restrict -c floatCycl.c -S -fcode-asm

#CC	= gcc
#CFLAGS 	= -O3 -Ofast -ftree-vectorize -mtune=intel -fopt-info-vec -fopt-info-vec-missed 
#CFLAGS 	= -O3 -ffast-math 

#CC	= pgcc
#CFLAGS 	= -fast -Minfo -Mvect

PRG = vect

SRC	= $(PRG).c \
	wallclock_time.c \

OBJ	= $(SRC:%.c=%.o)

.c.o: .c
	$(CC) $(CFLAGS) -c $<

$(PRG):	$(OBJ)
	$(CC) $(CFLAGS) -o $(PRG) $(OBJ) $(LIBS) $(DEFINES) 

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
