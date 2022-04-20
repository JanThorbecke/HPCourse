SHELL=/bin/bash

CC=icc
MAKE=make
AR=ar
ARCREATE=cr
CFLAGS=-DHAVE_uint -O

COMPILE=$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)

INCS =	bench.h lib_mem.h stats.h timing.h

SRCS =  lat_mem_rd.c

lat_mem_rd: getopt.o lib_mem.o lib_timing.o 
	$(COMPILE) lat_mem_rd.c -o lat_mem_rd $?

lib_timing.o : lib_timing.c $(INCS)
	$(COMPILE) -c lib_timing.c -o lib_timing.o
lib_mem.o : lib_mem.c $(INCS)
	$(COMPILE) -c lib_mem.c -o lib_mem.o
lib_unix.o : lib_unix.c $(INCS)
	$(COMPILE) -c lib_unix.c -o lib_unix.o
lib_debug.o : lib_debug.c $(INCS)
	$(COMPILE) -c lib_debug.c -o lib_debug.o
lib_stats.o : lib_stats.c $(INCS)
	$(COMPILE) -c lib_stats.c -o lib_stats.o
getopt.o : getopt.c $(INCS)
	$(COMPILE) -c getopt.c -o getopt.o


clean:
	/bin/rm -f *.[oas]

