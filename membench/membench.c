//////////////////////////////////////////////////////////////////////
// $Id: membench.c,v 1.2 1998/02/01 23:13:27 dmartin Exp $
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <sys/times.h>
#include <sys/types.h>

/* time function for Apple OSX */
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#define SAMPLE    10
#define CACHE_MIN (64*1024)
#define CACHE_MAX (32*1024*1024)

#define TIME_DIF_TO_NS(s,f) \
    ((f.tv_sec-s.tv_sec)*1000000000.0 + (f.tv_nsec-s.tv_nsec))

int x[CACHE_MAX];

void current_utc_time(struct timespec *ts) {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME, ts);
#endif
 
}


int main()
{
    int register i, index, stride, limit, temp;
    double sample_ns, sample_sec, sec, ns;
    struct timespec start,finish;
    int steps, tsteps, csize;
  
    for (csize = CACHE_MIN; csize <= CACHE_MAX; csize *= 2){
        for (stride = 1; stride <= csize/2; stride *= 2){
	    sec = 0;
	    ns = 0.0;
	    limit = csize-stride+1;
	    steps = 0;
	    do {
		current_utc_time(&start);
		for (i = SAMPLE*stride; i != 0; i--)
		    for (index = 0; index < limit; index += stride)
			x[index]++;
		current_utc_time(&finish);
		sample_ns = TIME_DIF_TO_NS(start,finish);
		sample_sec = sample_ns/1000000000.0;
		steps++;
		sec += sample_sec;
		ns += sample_ns;
	    } while (sec < 1.0);
	    tsteps=0;
	    do {
		current_utc_time(&start);
		for (i = SAMPLE*stride; i != 0; i--)
		    for (index = 0; index < limit; index += stride)
			temp += index;
		current_utc_time(&finish);
		sample_ns = TIME_DIF_TO_NS(start,finish);
		sample_sec = sample_ns/1000000000.0;
		tsteps++;
		sec -= sample_sec;
		ns -= sample_ns;
	    } while (tsteps < steps);
	    printf("Size: %7ld Stride: %7ld read+write: %14.0f ns\n",
		   csize*sizeof(int),
		   stride*sizeof(int), 
		   (double) (sec*1000000000.0)/(steps*SAMPLE*stride*((limit-1.0)/stride+1.0)));
	}
	printf ("\n");
    }
    return (0);
}

//////////////////////////////////////////////////////////////////////
// $Log: membench.c,v $
// Revision 1.2  1998/02/01 23:13:27  dmartin
// added exit(0)
//
// Revision 1.1  1998/01/19 00:47:39  dmartin
// Initial revision
//
//////////////////////////////////////////////////////////////////////

