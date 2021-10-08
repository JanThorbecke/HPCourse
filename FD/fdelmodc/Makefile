# Makefile

include ../Make_include

########################################################################
# define general include and system library
ALLINC     = -I.

#CFLAGS 	= -O3 -xP -ip -fno-alias -D$(MACHINE) $(ALLINC) -vec-report2
LIBS    = -L../FFTlib/lib -lgenfft -lm

PRG = fdelmodc

SRCC	= $(PRG).c \
		acoustic4.c \
		viscoacoustic4.c \
		elastic4.c \
		viscoelastic4.c \
		defineSource.c  \
		getParameters.c  \
		getWaveletInfo.c  \
		getModelInfo.c  \
		applySource.c  \
		getRecTimes.c  \
		writeSnapTimes.c  \
		writeRec.c  \
		fileOpen.c  \
		recvPar.c  \
		readModel.c  \
		taperEdges.c  \
		verbosepkg.c  \
		SU2su.c  \
		gaussGen.c  \
		spline3.c  \
		wallclock_time.c  \
		name_ext.c  \
		atopkge.c \
		docpkge.c \
		getpars.c

OBJC	= $(SRCC:%.c=%.o)

$(PRG):	$(OBJC) fdelmod.h
	$(CC) $(CFLAGS) -o fdelmodc $(OBJC) $(LIBS)

install: $(PRG)
	cp fdelmodc $(BIN)

clean:
		rm -f core $(OBJC) $(PRG) *.vect

print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar



