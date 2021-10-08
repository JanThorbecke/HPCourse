# Makefile

#CC	= icc
#CFLAGS 	= -O1 
#CFLAGS 	= -O1 -mssse3 -no-scalar-rep -opt-report 3 -opt-report-phase=hlo -restrict -S -fcode-asm
#CFLAGS 	= -O3 -mssse3 -par-report -opt-report 3 -opt-report-phase=hlo
#icc -O1 -mssse3 -no-scalar-rep -opt-report 3 -opt-report-phase=hlo -restrict -c floatCycl.c -S -fcode-asm

CC	= gcc
CFLAGS 	= -O1
#CFLAGS 	= -O3 -ffast-math 


PRG = floatCycl

SRC	= $(PRG).c \
	wallclock_time.c \

OBJ	= $(SRC:%.c=%.o)

.c.o: .c
	$(CC) $(CFLAGS) -c $<

$(PRG):	$(OBJ)
	$(CC) $(CFLAGS) -o $(PRG) $(OBJ) $(LIBS) $(DEFINES) -lm

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
