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
#CC	= gcc
#CFLAGS	= -O3 $(DEFINES) $(ALLINC) 

#OSX clang compiler
CC	= clang
CFLAGS 	= -O3 

#BLAS library
#marken (TA cluster)
MKL=/opt/intel/mkl/10.0.3.020/
LIBS = -L$(MKL)/lib/em64t/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm

#linux (TNW cluster)
MKL=/opt/mkl-8.0.1/lib/32/
LIBS = -L$(MKL) -lmkl_solver -lmkl_lapack -lmkl_ia32 -lguide -lpthread -lm

#Linux with Atlas installed using latest GNU compilers
LIBS = -L/usr/local/atlas/lib -latlas -lf77blas -lgfortran -lm
# using older GNU compilers
#LIBS = -L/usr/local/atlas/lib -latlas -lf77blas -lg2c

#Mac OSX with intel compiler
MKL=/opt/intel/Compiler/11.0/059/Frameworks/mkl/lib/em64t
MKL=/opt/intel/Compiler/11.1/046/Frameworks/mkl/lib/universal
LIBS = -L$(MKL) -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -lm

#Mac OSX with M1 aarch64 chips
LIBS= /Library/Developer/CommandLineTools/SDKs/MacOSX12.1.sdk//usr/lib/liblapack.tbd /Library/Developer/CommandLineTools/SDKs/MacOSX12.1.sdk//usr/lib/libblas.tbd


#if you have none of the above just try
#LIBS= -lblas

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
