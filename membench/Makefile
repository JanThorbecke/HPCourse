# Makefile

#CFLAGS 	= -O3 -mssse3 -par-report -opt-report 3 -opt-report-phase=hlo
#icc -O1 -mssse3 -no-scalar-rep -opt-report 3 -opt-report-phase=hlo -restrict -c floatCycl.c -S -fcode-asm

CC	= gcc
CFLAGS 	= -O0 

PRG = membench

SRC	= $(PRG:%=%.c)
EXE	= $(PRG:%=%.exe
OBJ	= $(SRC:%.c=%.o)

%.o: %.c 
	$(CC) $(CFLAGS) -c $<

%.exe: %.o  $(OBJ)
	#$(CC) $(CFLAGS) -o $@ $(OBJ) 

install: $(PRG)

clean:
		rm -f core $(OBJ) $(PRG)

print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar
