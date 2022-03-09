# Makefile

CC	= gcc
CFLAGS 	= -O0 -m64

PRG = loadStore
#PRG = lsTests

SRC	= $(PRG).c 

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
