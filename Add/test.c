// TAB -- The Ampersand BUG function
// Returns a pointer to an int
#include<stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
double wallclock_time(void);

int* TAB() {
   int temp[10];
   return(&temp[11]); // returns a pointer to the local int
}

void Victim() {
   int* ptr;
   ptr = TAB();
   *ptr = 42;
	fprintf(stderr,"ptr =%d\n", *ptr);
}

int main(int argc, char *argv[]) {
    size_t k, i, j, N, Loop, mhz, len, cycles;
    float *a, *b, *c;
    double t0, t1, t2, fcycles;

    mhz = 2500;
    Loop = 100000000;
    N=32;
    a=(float *)calloc(N,sizeof(float));
    b=(float *)calloc(N,sizeof(float));
    c=(float *)calloc(N,sizeof(float));

    for (i=0; i<N; i++) {
        a[i] = (float) rand()/RAND_MAX;
        b[i] = (float) rand()/RAND_MAX;
        c[i] = (float) rand()/RAND_MAX;
    }


	t0=wallclock_time();
	for (k=0; k<Loop; k++){
		for (j=0; j<N; j++) {
			c[j] = pow(b[j],1.5);
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"float pow time=%f gives %f Mflop/s\n", t1-t0, ((1.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(1*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(1.0*(Loop/1000000)*N);
	fprintf(stderr,"float pow used %d cycles per instruction %f\n", (int)cycles, fcycles);

	return 0;
}

double wallclock_time(void)
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
		   (double)(1e-6*((double)s_val.tv_usec-(double)b_val.tv_usec));

	return (double)time;
}

double wallclock_time_(void)
{
	return (double)wallclock_time();
}

