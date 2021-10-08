#
#  Makefile for Poisson solver - Peter Michielse, SGI
#  Code obtained from Roel Verstappen
#  Jan. 18, 1999
#

.SUFFIXES: .o1 .o2 .o3

FC	= /usr/bin/f77
CC	= /usr/bin/cc

EXE3	= poisson_O3

OPT_O3	= -c -O3 -r8

PAR     = -pfa keep

LIB	=

LINK_O3	= -mp

OBJ_O3 = \
bitr.o3 ifft.o3 rhs.o3 fft.o3 inifft.o3 solu.o3 grid.o3 poissn.o3 \
poisson.o3 iccgp.o3 precon.o3 verify.o3 wclock.o

all: $(EXE3)

.IGNORE:

$(OBJ_O3) :

$(EXE3): $(OBJ_O3)
	$(FC) $(LINK_O3) -o $(EXE3) $(OBJ_O3) $(LIB)

.f.o3 :
	$(FC) -o $*.o3 $(OPT_O3) $(PAR) $<

clean:
	@rm -f $(EXE3) $(OBJ_O3)
