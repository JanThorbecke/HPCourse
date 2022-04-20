#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

double wallclock_time(void);

#define MIN(x,y) ((x) < (y) ? (x) : (y))

int main ()
{
	FILE *fp;
	int mhz=3000, l3_size=67108864, l2_size=1048576, l1_size=32768, l2_line=64;
	int ix, iy, j;
	size_t power, size, loops;
	double t0, t1, t_bias, time;
	double *a, *c, *d, zero=0, *f, *t;
	char *filename;

	filename="DelftBlue";
	t0  = wallclock_time();
	t1  = wallclock_time();
	t_bias = t1-t0;
	fprintf(stderr,"t_bias = %e\n", t_bias);
	a  = (double *) malloc(l2_size);
	fp = fopen(filename, "w+");
	
/************************************/

	size = 0;
	power = 1;
	while (power < 32) {
		size = pow(2,power);

		c = (double *)malloc(size*sizeof(double)+l2_line);
		assert(c!= NULL);
		d = (double *)malloc(size*sizeof(double)+l2_line);
		assert(d!= NULL);

		for (ix=0; ix<size; ix++) d[ix] = ix*M_PI;

		for (ix=0; ix<MIN(l2_size/(sizeof(double)),size); ix++ ) c[ix] = d[ix];
		for (ix=0; ix<MIN(l1_size/(sizeof(double)),size); ix++ ) c[ix] = d[ix];

		loops = l2_size/size;
		if (loops<=2) loops=3;
		fprintf(stderr,"size = %ld power = %ld number of loops=%ld\n", size*sizeof(double), power, loops);
		t0  = wallclock_time();
		for (j=0; j<loops; j++) {
			for (ix=0; ix<size; ix+=1) {
				c[ix] = d[ix];
//				c[ix+1] = d[ix+1];
//				c[ix+2] = d[ix+2];
//				c[ix+3] = d[ix+3];
//				c[ix+4] = d[ix+4];
//				c[ix+5] = d[ix+5];
//				c[ix+6] = d[ix+6];
//				c[ix+7] = d[ix+7];
			}
		}
		t1  = wallclock_time();
	
		time = (t1-t0-t_bias)/loops;
		fprintf(stderr, "double copy time = %e [%ld] = %ld cycles = %f MB/s\n", 
			time, loops, (size_t)((time)*mhz*1e+6), size*sizeof(double)/((1024*1024)*(time)) );
	
		fprintf(fp, "%ld %e\n", (size*sizeof(double)), size*sizeof(double)/((1024*1024)*(time)) );
	
		free(c);
		free(d);
		power++;
	}
	fclose(fp);

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

