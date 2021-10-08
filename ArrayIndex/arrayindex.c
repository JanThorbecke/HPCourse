#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

double wallclock_time(void);

int main ()
{
	FILE *fp;
	int mhz=2500, l2_size=6291456, l1_size=32768, l2_line=64;
	int j, *index;
	int N, ix, iy;
	int power, size, loops;
	double t0, t1, t_bias, time;
	double *a, *c, **d;
	char *filename;

	t0  = wallclock_time();
	t1  = wallclock_time();
	t_bias = t1-t0;
	fprintf(stderr,"t_bias = %e\n", t_bias);
	
/************************************/

	N = 2048;
	size = N*N;
	loops = 100;

	d = (double **)malloc(N*sizeof(double*));
	assert(d != NULL);
	d[0] = (double *)malloc(size*sizeof(double));
	assert(d[0] != NULL);

	for (ix=1; ix<N; ix++) {
		d[ix] = d[0]+N*ix;
//		printf("&d[%d] = %p\n", ix, &d[ix][0]); 
	}

	c = &d[0][0];

	index = (int *)malloc(size*sizeof(int));
	assert(index!= NULL);
	for (ix=0; ix<size; ix++) {
		index[ix] = ix;
	}

	t0  = wallclock_time();
	for (j=0; j<loops; j++) {
		d[0][0] = j;
		for (iy=0; iy<N; iy++) {
			for (ix=0; ix<N; ix++) {
				d[iy][ix] = j+1.1;
			}
		}
	}
	t1  = wallclock_time();
	time = (t1-t0-t_bias)/loops;
	fprintf(stderr, "Direct index time = %e [%d] = %d cycles = %f MB/s\n", 
		time, loops, (int)((time)*mhz*1e+6), d[129][313]);
	
	
	t0  = wallclock_time();
	for (j=0; j<loops; j++) {
		d[0][0] = j;
		for (iy=0; iy<N; iy++) {
			for (ix=0; ix<N; ix++) {
				c[iy*N+ix] = j+1.1;
			}
		}
	}
	t1  = wallclock_time();
	time = (t1-t0-t_bias)/loops;
	fprintf(stderr, "Explicit index time = %e [%d] = %d cycles = %f MB/s\n", 
		time, loops, (int)((time)*mhz*1e+6), d[129][313]);

	t0  = wallclock_time();
	for (j=0; j<loops; j++) {
		d[0][0] = j;
		iy=0;
		for (ix=0; ix<size; ix++) {
			c[iy] = j+1.1;
            iy = iy+1;
		}
	}
	t1  = wallclock_time();
	time = (t1-t0-t_bias)/loops;
	fprintf(stderr, "Loop carried index time = %e [%d] = %d cycles = %f MB/s\n", 
		time, loops, (int)((time)*mhz*1e+6), d[129][313]);
	
	t0  = wallclock_time();
	for (j=0; j<loops; j++) {
		d[0][0] = j;
		for (ix=0; ix<size; ix++) {
			c[index[ix]] = j+1.1;
//			c[ix] = ix*M_PI;
		}
	}
	t1  = wallclock_time();
	time = (t1-t0-t_bias)/loops;
	fprintf(stderr, "Indirect time = %e [%d] = %d cycles = %f MB/s\n", 
		time, loops, (int)((time)*mhz*1e+6), d[129][313]);
	
	free(d[0]);

	return 0;
}

double wallclock_time(void )
{
    struct timeval s_val;
    static struct timeval b_val;
    double time;
    static int base=0;

    gettimeofday(&s_val,0);

    if (!base) {
        b_val = s_val;
        base = 1;
        return 0.0;
    }

    time = (double)(s_val.tv_sec-b_val.tv_sec) +
           (double)(1e-6*(s_val.tv_usec-b_val.tv_usec));

    return (double)time;
}

