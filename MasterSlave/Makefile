# Makefile

#MPI compiler 
MPICC	= mpicc
MPIF90	= mpif90
#MPICC	= gcc
MPICFLAGS 	= -O3 -mpe=mpilog
MPICFLAGS 	= -O3 -I/opt/mpich2/include  -m32
MPICFLAGS 	= -O3  -static
MPIFFLAGS 	= -O3  -static

#MPILIBS=-L/opt/mpich2/lib -lfmpich /u/thorbec/src/mpe2-1.0.6p1/lib/liblmpe.a /u/thorbec/src/mpe2-1.0.6p1/lib/libmpe.a -lmpich -lpthread -lrt
#MPILIBS=-L/u/thorbec/src/fpmpi-2.1f -lfpmpi

all:  MPIprobe MasterSlave

MPIprobe: MPIprobe.c
	$(MPICC) $(MPICFLAGS) -o $@ $@.c $(MPILIBS) -lm

MasterSlave: MasterSlave.c
	$(MPICC) $(MPICFLAGS) $@.c -o $@ $(MPILIBS) -lm

.c.o: .c
	$(MPICC) $(MPICFLAGS) -c $<

clean:
	rm -f core MPIprobe MasterSlave

print:	Makefile $(SRC)
	$(PRINT) $?
	@touch print

count:
	@wc $(SRC)

tar:
	@tar cf $(PRG).tar Makefile $(SRC) && compress $(PRG).tar
