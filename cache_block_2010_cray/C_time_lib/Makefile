#
# source options -DCRAY_USER_WTIME="function id" -DCRAY_NO_INIT_CALL
#                -DCRAY_USER_WTIME_F="fortran id" -DCRAY_LIST_CLOCK_RATE
#
SRCOPT = -DCRAY_NO_INIT_CALL -DCRAY_USER_WTIME_F=secondr_ -DCRAY_USER_WTIME=secondr

LIB = libCray_timer.a

CC   = gcc
CCFLAGS= -O2 $(SRCOPT)

FC   = ftn

AR  = ar -cr
#

$(LIB): rtc.o cray_rt_clock.o
	$(AR) $(LIB) rtc.o cray_rt_clock.o

rtc.o: rtc.c
	$(CC) -c $(CCFLAGS) rtc.c

cray_rt_clock.o: cray_rt_clock.c
	$(CC) -c $(CCFLAGS) cray_rt_clock.c

clean:
	rm -f rtc.o cray_rt_clock.o $(LIB)
