# Makefile for corrtest

CC = icc
CFLAGS  = -O3 -vec-report
#CC = gcc
#CFLAGS  = -O3 -msse -ftree-vectorize -ftree-vectorizer-verbose=2 -fstrict-aliasing 
#CC = clang
#CFLAGS  = -Ofast  -fdiagnostics-show-option
#CC = pgcc
#CFLAGS  = -fast -Minfo -Mvect=simd:128

OBJ	= atopkge.o \
	  docpkge.o \
	  getpars.o \
	  wallclock_time.o

MOBJ = corrOpt.o

all: corrOpt

corrOpt: $(OBJ) $(MOBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(MOBJ) -lm

clean:
		rm -f core $(OBJ) $(MOBJ)

print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar

INCLUDE_FILES = segy.h

.SUFFIXES: .o .c .C .F .f

.F.o: $(INCLUDE_FILES)
	$(F77) $(FFLAGS) $(DEFINES) $(ALLINC) -c $<

.f.o: $(INCLUDE_FILES)
	$(F77) $(FFLAGS) $(DEFINES) $(ALLINC) -c $<

.C.o: $(INCLUDE_FILES)
	$(CPLUS) $(CFLAGS) $(DEFINES) $(ALLINC) -c $<

.c.o: $(INCLUDE_FILES)
	$(CC) $(CFLAGS) $(DEFINES) $(ALLINC) -c $<

