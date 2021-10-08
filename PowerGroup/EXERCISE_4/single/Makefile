#
#  Makefile for Jacobi solver - Peter Michielse, SGI
#  March 30, 1999

.SUFFIXES: .o

FC      = /usr/bin/f77
CC      = /usr/bin/cc

EXE     = jacobi

OPT_O3  = -c -O3

LIB     = 

LINK    =

OBJ     =  driver.o error_check.o initialize.o jacobi.o main.o

all: $(EXE)

$(OBJ):	param.h

$(EXE): $(OBJ)
	$(FC) $(LINK) -o $(EXE) $(OBJ) $(LIB)

.f.o:
	$(FC) -o $*.o $(OPT_O3) $(PAR) $<

clean:
	@rm -f $(EXE) $(OBJ)
