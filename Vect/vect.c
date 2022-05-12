#include<stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include<math.h>
#include<stdlib.h>

double wallclock_time(void);

int main(int argc, char *argv[]) {
    size_t k, i, j, N, Loop, mhz, len, cycles;
    float *a, *b;
    double t0, t1, fcycles, s;

    mhz = 3200;
    Loop = 1000000;
    N=4096;
    a=(float *)calloc(N,sizeof(float));
    b=(float *)calloc(N,sizeof(float));

    for (i=0; i<N; i++) {
        a[i] = (float) rand()/RAND_MAX;
        b[i] = (float) rand()/RAND_MAX;
    }

	t0=wallclock_time();
	s = 0.0;
#pragma novector
	for (k=0; k<Loop; k++){
/* this loop will not vectorize */
#pragma novector
#pragma clang loop vectorize(disable)
		for (j=0; j<N; j++) {
			s += a[j]*b[j];
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"non-vector time=%f gives %f Mflop/s\n", t1-t0, ((2.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(2*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(2.0*(Loop/1000000)*N);
	fprintf(stderr,"loop used %d cycles per instruction %f\n", (int)cycles, fcycles);
	fprintf(stderr,"s %f\n", s);

	t0=wallclock_time();
	s = 0.0;
	for (k=0; k<Loop; k++){
/* this loop will vectorize */
#pragma ivdep
#pragma vector
		for (j=0; j<N; j++) {
			s += a[j]*b[j];
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"vector time=%f gives %f Mflop/s\n", t1-t0, ((2.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(2*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(2.0*(Loop/1000000)*N);
	fprintf(stderr,"loop used %d cycles per instruction %f\n", (int)cycles, fcycles);
	fprintf(stderr,"s %f\n", s);

/* a not correct vectorizing example  */
	t0=wallclock_time();
	s = 0.0;
	for (k=0; k<Loop; k++){
/* this loop will vectorize */
#pragma ivdep
#pragma vector
#pragma clang loop vectorize(enable)
		for (j=0; j<N; j++) {
			a[j+1] = b[j];
			b[j+1] = a[j];
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"vector time=%f gives %f Mflop/s\n", t1-t0, ((2.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(2*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(2.0*(Loop/1000000)*N);
	fprintf(stderr,"loop used %d cycles per instruction %f\n", (int)cycles, fcycles);
	fprintf(stderr,"a[2] %e\n", a[2]);

	return 0;
}

