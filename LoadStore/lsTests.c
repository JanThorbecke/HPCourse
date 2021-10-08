#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <blas.h>

float wallclock_time(void);

typedef struct _complexStruct { /* complex number */
	float r,i;
} complex;

#define MIN(x,y) ((x) < (y) ? (x) : (y))

void main ()
{
	int ix, iy, j, l2_size=1048576, l1_size=32768;
	int loops, power, size, l2_line=256;
	float t0, t1, t_bias, time;
	float *a, *c, *d, zero=0;


	t0  = wallclock_time();
	t1  = wallclock_time();
	t_bias = t1-t0;
	a = (float *) malloc(l2_size);

/************************************/

	size = 0;
	power = 8;
	while (size < pow(2,27) ) {
		size = pow(2,power);
		fprintf(stderr,"size = %ld power = %d\n", size*sizeof(float), power);

		posix_memalign(&c, l2_line, size*sizeof(float)+l2_line);
		assert(c!= NULL);
		posix_memalign(&d, l2_line, size*sizeof(float)+l2_line);
		assert(d!= NULL);

		for (ix=0; ix<size; ix++) d[ix] = ix*M_PI;
		/* clean L2-cache */
		for (ix=0; ix<l2_size/(sizeof(float)); ix++) a[ix] = 0.0;
		
		loops = l2_size/size;
		if (!loops) loops=1;

		t0  = wallclock_time();
		for (j=0; j<loops; j++) {
		#pragma ivdep
		for (ix=0; ix<size; ix+=7) {
			c[ix] = d[ix];
			c[ix+1] = d[ix+1];
			c[ix+2] = d[ix+2];
			c[ix+3] = d[ix+3];
			c[ix+4] = d[ix+4];
			c[ix+5] = d[ix+5];
			c[ix+6] = d[ix+6];
		}
		}
		t1  = wallclock_time();
		time = (t1-t0-t_bias)/loops;

		fprintf(stderr, "float copy time = %e = %d cycles = %f MB/s\n", 
	time, (int)((time)*250e+6), size*sizeof(float)/((1024*1024)*(time)) );

	fprintf(stdout, "%ld %e\n", (size*sizeof(float)), size*sizeof(float)/((1024*1024)*(time)) );

		free(c);
		free(d);
		power++;
	}

fprintf(stdout, "\n\n");

/************************************/
/************************************/

	size = 0;
	power = 8;
	while (size < pow(2,27) ) {
	size = pow(2,power);
	fprintf(stderr,"size = %ld power = %d\n", size*sizeof(float), power);

		posix_memalign(&c, l2_line, size*sizeof(float)+l2_line);
		assert(c!= NULL);
		posix_memalign(&d, l2_line, size*sizeof(float)+l2_line);
		assert(d!= NULL);

		for (ix=0; ix<size; ix++) d[ix] = ix*M_PI;
		for (ix=0; ix<l2_size/(sizeof(float)); ix++) a[ix] = 0.0;

		loops = l2_size/size;
		if (!loops) loops=1;

		t0  = wallclock_time();
		for (j=0; j<loops; j++) {
			memcpy(c, d, size*sizeof(float));
		}
		t1  = wallclock_time();
		time = (t1-t0-t_bias)/loops;

	fprintf(stderr, "float memcpy time = %e = %d cycles = %f MB/s\n", 
	time, (int)((time)*250e+6), size*sizeof(float)/((1024*1024)*(time)) );

	fprintf(stdout, "%ld %e\n", (size*sizeof(float)), size*sizeof(float)/((1024*1024)*(time)) );

	free(c);
	free(d);
	power++;
	}

/************************************/

fprintf(stdout, "\n\n");

	size = 0;
	power = 8;
	while (size < pow(2,27) ) {
		size = pow(2,power);
		fprintf(stderr,"size = %ld power = %d\n", size*sizeof(float), power);

		posix_memalign(&c, l2_line, size*sizeof(float)+l2_line);
		assert(c!= NULL);
		posix_memalign(&d, l2_line, size*sizeof(float)+l2_line);
		assert(d!= NULL);

		for (ix=0; ix<size; ix++) d[ix] = ix*M_PI;
		for (ix=0; ix<l2_size/(sizeof(float)); ix++) a[ix] = 0.0;
	
		loops = l2_size/size;
		if (!loops) loops=1;

		t0  = wallclock_time();
		for (j=0; j<loops; j++) {
			bcopy(d, c, size*sizeof(float));
		}
		t1  = wallclock_time();
		time = (t1-t0-t_bias)/loops;

		fprintf(stderr, "bcopy copy time = %e = %d cycles = %f MB/s\n", 
	time, (int)((time)*250e+6), size*sizeof(float)/((1024*1024)*(time)) );

	fprintf(stdout, "%ld %e\n", (size*sizeof(float)), size*sizeof(float)/((1024*1024)*(time)) );

		free(c);
		free(d);
		power++;
	}

/************************************/

fprintf(stdout, "\n\n");

	size = 0;
	power = 8;
	while (size < pow(2,27) ) {
		size = pow(2,power);
		fprintf(stderr,"size = %ld power = %d\n", size*sizeof(float), power);

		posix_memalign(&c, l2_line, size*sizeof(float)+l2_line);
		assert(c!= NULL);
		posix_memalign(&d, l2_line, size*sizeof(float)+l2_line);
		assert(d!= NULL);

		for (ix=0; ix<size; ix++) d[ix] = ix*M_PI;
		for (ix=0; ix<l2_size/(sizeof(float)); ix++) a[ix] = 0.0;
	
		loops = l2_size/size;
		if (!loops) loops=1;

		t0  = wallclock_time();
		for (j=0; j<loops; j++) {
			scopy(size, c, 1, d, 1);
		}
		t1  = wallclock_time();
		time = (t1-t0-t_bias)/loops;

		fprintf(stderr, "blas scopy copy time = %e = %d cycles = %f MB/s\n", 
	time, (int)((time)*250e+6), size*sizeof(float)/((1024*1024)*(time)) );

	fprintf(stdout, "%ld %e\n", (size*sizeof(float)), size*sizeof(float)/((1024*1024)*(time)) );

		free(c);
		free(d);
		power++;
	}

	return;
}



/*################################################*/
/*################################################*/
/*################################################*/
/*################################################*/
/*################################################*/
/*################################################*/

float wallclock_time(void )
{
    struct timeval s_val;
    static struct timeval b_val;
    float time;
    static int base=0;

    gettimeofday(&s_val,0);

    if (!base) {
        b_val = s_val;
        base = 1;
        return 0.0;
    }

    time = (float)(s_val.tv_sec-b_val.tv_sec) +
           (float)(1e-6*(s_val.tv_usec-b_val.tv_usec));

    return (float)time;
}

