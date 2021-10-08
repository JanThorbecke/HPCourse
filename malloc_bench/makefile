CC	= cc
COPT	= -O3
LIBS	= -L~/lus/tcmalloc/lib -ltcmalloc_minimal

OBJS	= malloc_bench.o
BINS	= malloc_bench tcmalloc_bench

.c.o:
	$(CC) $(COPT) -c -o $@ $*.c

default: $(BINS)

malloc_bench:   $(OBJS)
	$(CC) $(COPT) -o $@ $(OBJS)
tcmalloc_bench: $(OBJS)
	$(CC) $(COPT) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f *.o *.oo core* $(BINS)
