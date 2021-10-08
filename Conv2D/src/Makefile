# Makefile for puls3D

ROOT=$(HOME)/src/HPCourse/Conv2D/

include $(ROOT)/Make_include

ALLINC = -I$(ROOT)/FFTlib/include/

CFLAGS += $(ALLINC) 

LIBS    = -L../FFTlib/lib -lgenfft -lm

PRG	= puls3D

OBJ	= puls3D.o \
	  extr3D.o \
	  tablecalc_2D.o \
	  conv2D.o \
	  getpars.o \
	  docpkge.o \
	  atopkge.o

$(PRG):	$(OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $(ALLINC) $(OBJ) $(LIBS)

clean:
		rm -f core $(OBJ)

print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar

