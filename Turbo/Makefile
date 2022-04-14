# Makefile


#intel compiler
CC	= icc
CFLAGS 	= -O3 -opt-report
CFLAGS 	= -O3 -par-report -vec-report -Wl,-no_compact_unwind
#CFLAGS 	= -O1 -par-report

#PGI compiler
#CC = pgcc
#CFLAGS  = -fast

#GNU compiler
CC	= gcc
CFLAGS	= -O3 $(DEFINES) $(ALLINC) -fopenmp

#OSX clang compiler
#CC	= clang
#CFLAGS 	= -O3 

#if you have none of the above just try
LIBS= -lblas

#If you still do not have blas and get messages about missing sgemm_ 
#uncomment the following line
#LIBS= -lm

PRG = ClassicMatrixMultiply

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
