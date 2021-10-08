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
	int mhz=2500;
	int i, loops;
	double t0, t1, t_bias, time;
	float a[2][2], b[2][2], c[2][2];
	char *filename;

	filename="MacCore2Duo";
	t0  = wallclock_time();
	t1  = wallclock_time();
	t_bias = t1-t0;
	fprintf(stderr,"t_bias = %e\n", t_bias);
	fp = fopen(filename, "w+");

	a[0][0] = -1.00000001; a[0][1] = -1.0;
	a[1][0] =  1.0; a[1][1] =  1.0;
	
	b[0][0] =  2.0; b[0][1] =  2.0;
	b[1][0] = -2.0; b[1][1] = -2.0;

	c[0][0] =  5.0; c[0][1] =  5.0;
	c[1][0] = -5.0; c[1][1] = -5.0;

/************************************/

	t0  = wallclock_time();

	loops = 400000000;
	for (i=0; i<loops; i++) {
/* C = C + AB */
	c[0][0] = c[0][0] + a[0][0]*b[0][0]
					  + a[1][0]*b[0][1];
	c[0][1] = c[0][1] + a[0][1]*b[0][0]
					  + a[1][1]*b[0][1];
	c[1][0] = c[1][0] + a[0][0]*b[1][0]
					  + a[1][0]*b[1][1];
	c[1][0] = c[1][0] + a[0][1]*b[1][0]
					  + a[1][1]*b[1][1];
	
/* A = A + BC */
	a[0][0] = a[0][0] + b[0][0]*c[0][0]
					  + b[1][0]*c[0][1];
	a[0][1] = a[0][1] + b[0][1]*c[0][0]
					  + b[1][1]*c[0][1];
	a[1][0] = a[1][0] + b[0][0]*c[1][0]
					  + b[1][0]*c[1][1];
	a[1][0] = a[1][0] + b[0][1]*c[1][0]
					  + b[1][1]*c[1][1];
	
/* B = B + CA */
	b[0][0] = b[0][0] + c[0][0]*a[0][0]
					  + c[1][0]*a[0][1];
	b[0][1] = b[0][1] + c[0][1]*a[0][0]
					  + c[1][1]*a[0][1];
	b[1][0] = b[1][0] + c[0][0]*a[1][0]
					  + c[1][0]*a[1][1];
	b[1][0] = b[1][0] + c[0][1]*a[1][0]
					  + c[1][1]*a[1][1];
	
	}
	t1  = wallclock_time();

	fprintf(stderr,"A = | %f %f |\n    | %f %f | \n", a[0][0], a[1][0], a[0][1], a[1][1]);
	fprintf(stderr,"B = | %f %f |\n    | %f %f | \n", b[0][0], b[1][0], b[0][1], b[1][1]);
	fprintf(stderr,"C = | %f %f |\n    | %f %f | \n", c[0][0], c[1][0], c[0][1], c[1][1]);

	time = (t1-t0-t_bias);
	fprintf(stderr, "performance time = %e  mflop = %f Mflop/s\n", 
			time, (48.0*(loops/1e6))/(time) );
	
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

