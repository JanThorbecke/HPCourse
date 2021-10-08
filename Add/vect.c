// TAB -- The Ampersand BUG function
// Returns a pointer to an int
#include<stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include<math.h>
#include<stdlib.h>

double wallclock_time(void);

int main(int argc, char *argv[]) {
    size_t k, i, j, N, Loop, mhz, len, cycles;
    float *a, *b, s;
    double t0, t1, t2, fcycles;

    mhz = 2500;
    Loop = 10000000;
    N=1024;
    a=(float *)calloc(N,sizeof(float));
    b=(float *)calloc(N,sizeof(float));

    for (i=0; i<N; i++) {
        a[i] = (float) rand()/RAND_MAX;
        b[i] = (float) rand()/RAND_MAX;
    }

	t0=wallclock_time();
	s = 0.0;
	for (k=0; k<Loop; k++){
#pragma novector
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

	return 0;
}

