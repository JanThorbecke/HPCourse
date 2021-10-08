/*
 * malloc_bench.c
 * (c)2009 Seiji Nishimura
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#define DEFAULT_NTIMES	1000000
#define DEFAULT_NBYTES	1048576

#define MIN(x,y)	(((x)<(y))?(x):(y))
#define MAX(x,y)	(((x)>(y))?(x):(y))

typedef struct {
    int ntimes;
    size_t nbytes;
} prms_t;

// function prototypes
void prms_init(prms_t * prms, int argc, char **argv);
double malloc_test(int ntimes, size_t nbytes);
inline double wtime(void);

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

int main(int argc, char **argv)
{
    prms_t prms;

    prms_init(&prms, argc, argv);

    printf("MALLOC Bench: %d x %ld [Bytes], TIME=%g[sec.]\n",
	   prms.ntimes, prms.nbytes,
	   malloc_test(prms.ntimes, prms.nbytes));

    return EXIT_SUCCESS;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

void prms_init(prms_t * prms, int argc, char **argv)
{
    int option;
    const char *optstring = "m:n:";

    // zero clear all elements of prms
    memset(prms, 0x00, sizeof(prms_t));

    // defaults
    prms->ntimes = DEFAULT_NTIMES;
    prms->nbytes = DEFAULT_NBYTES;

    // command line analysis
    while ((option = getopt(argc, argv, optstring)) != -1) {
	switch (option) {
	case ('m'):		// memory size
	    prms->nbytes = MAX(1, labs(atol(optarg)));
	    break;
	case ('n'):		// number of iteration
	    prms->ntimes = MAX(1, abs(atoi(optarg)));
	    break;
	case ('?'):
	default:		// unknown option
	    exit(EXIT_FAILURE);
	    break;
	}
    }

    return;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

double malloc_test(int ntimes, size_t nbytes)
{
    int i;
    void *ptr;
    double tst, ted;

    tst = wtime();

#pragma omp parallel for private(i,ptr)
    for (i = 0; i < ntimes; i++) {
	// allocation
	if ((ptr = malloc(nbytes)) == NULL) {
	    perror(__FUNCTION__);
#ifndef _OPENMP
	    exit(EXIT_FAILURE);
#endif
	}
	// touch all pages of the allocated memory
	memset(ptr, 0x00, nbytes);
	// deallocation
	free(ptr);
    }

    ted = wtime();

    return ted - tst;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

inline double wtime(void)
{
    struct timeval tv;

    gettimeofday(&tv, NULL);

    return ((double) tv.tv_sec + 1.E-6 * tv.tv_usec);
}
