# Makefile

CC	= gcc
CFLAGS 	= -O3 -m64 -mssse3 -ftree-vectorize -ftree-vectorizer-verbose=2

CC	= icc
CFLAGS 	= -O3 -mssse3 -vec_report3
CFLAGS 	= -O3 -mssse3 -opt-report -vec_report
CFLAGS 	= -O3 -mssse3 -vec_report

PRG = Peak
#PRG = lsTests

SRC	= $(PRG).c 

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
