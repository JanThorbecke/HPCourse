# Makefile

CC	= gcc
CFLAGS 	= -O1
#CFLAGS 	= -O3 -ffast-math 

#CC	= icc
#CFLAGS 	= -O1 
#CFLAGS 	= -O3 -march=core-avx2 -qmkl -qopt-report -qopt-report-annotate -qopt-report-phase=all
#icc -O1 -mssse3 -no-scalar-rep -opt-report 3 -opt-report-phase=hlo -restrict -c floatCycl.c -S -fcode-asm


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
