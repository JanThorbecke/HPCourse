#
#
# 
#

LIB = lib_timer.a

CC   = gcc
CCFLAGS= -O2 

FC   = gfortran

AR  = ar -cr
#

$(LIB): wallclock_time.o 
	$(AR) $(LIB) wallclock_time.o 

wallclock_time.o: wallclock_time.c
	$(CC) -c $(CCFLAGS) wallclock_time.c

clean:
	rm -f wallclock_time.o $(LIB)
