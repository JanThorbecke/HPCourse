#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double wallclock_time(void);

int main (int argc, char *argv[])
{
	size_t k, i, j, N, Loop, mhz, len, cycles;
	float *a, *b, *c;
	double t0, t1, t2, fcycles;

	mhz = 2000;
	Loop = 10000000;
	N=128;
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
			c[j] = b[j]*a[j];
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"float mul time=%f gives %f Mflop/s\n", t1-t0, ((1.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(1*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(1.0*(Loop/1000000)*N);
	fprintf(stderr,"float mul used %d cycles per instruction (%f)\n", (int)cycles, fcycles);


	t0=wallclock_time();
	for (k=0; k<Loop; k++){
		for (j=0; j<N; j++) {
			c[j] = c[j]+b[j]*a[j];
		}
	}
	
	t1=wallclock_time();
	fprintf(stderr,"float add mul time=%f gives %f Mflop/s\n", t1-t0, ((2.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(2*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(2.0*(Loop/1000000)*N);
	fprintf(stderr,"float add mul used %d cycles per instruction %f\n", (int)cycles, fcycles);

	t0=wallclock_time();
	for (k=0; k<Loop; k++){
		for (j=0; j<N; j+=4) {
			c[j] = c[j]+b[j]*a[j];
			c[j+1] = c[j+1]+b[j+1]*a[j+1];
			c[j+2] = c[j+2]+b[j+2]*a[j+2];
			c[j+3] = c[j+3]+b[j+3]*a[j+3];
//			c[j+4] = c[j+4]+b[j+4]*a[j+4];
//			c[j+5] = c[j+5]+b[j+5]*a[j+5];
//			c[j+6] = c[j+6]+b[j+6]*a[j+6];
//			c[j+7] = c[j+7]+b[j+7]*a[j+7];
		}
	}
	t1=wallclock_time();
//	fprintf(stderr,"c=%f %f %f %f\n",c[0], c[N/2], c[N/4], c[12]);
	fprintf(stderr,"unrolled float add mul time=%f gives %f Mflop/s\n", t1-t0, ((2.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(2*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(2.0*(Loop/1000000)*N);
	fprintf(stderr,"unrolled float add mul used %d cycles per instruction %f\n", (int)cycles, fcycles);

	t0=wallclock_time();
	t0=wallclock_time();
	for (k=0; k<Loop; k++){
		for (j=0; j<N; j++) {
			c[j] = b[j]/a[j];
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"float div time=%f gives %f Mflop/s\n", t1-t0, ((1.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(1*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(1.0*(Loop/1000000)*N);
	fprintf(stderr,"float div used %d cycles per instruction %f\n", (int)cycles, fcycles);

	t0=wallclock_time();
	for (k=0; k<Loop; k++){
		for (j=0; j<N; j++) {
			c[j] = pow(b[j],2.1);
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"float pow time=%f gives %f Mflop/s\n", t1-t0, ((1.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(1*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(1.0*(Loop/1000000)*N);
	fprintf(stderr,"float pow used %d cycles per instruction %f\n", (int)cycles, fcycles);

	t0=wallclock_time();
	for (k=0; k<Loop; k++){
		for (j=0; j<N; j++) {
			c[j] = b[j]*b[j];
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"float mul^2 time=%f gives %f Mflop/s\n", t1-t0, ((1.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(1*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(1.0*(Loop/1000000)*N);
	fprintf(stderr,"float mul used %d cycles per instruction %f\n", (int)cycles, fcycles);

	Loop = Loop/10;
	t0=wallclock_time();
	for (k=0; k<Loop; k++){
		for (j=0; j<N; j++) {
			c[j] = sin(b[j]);
		}
	}
	t1=wallclock_time();
	fprintf(stderr,"float sin time=%f gives %f Mop/s\n", t1-t0, ((1.0*(Loop/1000000)*N))/(t1-t0));
	cycles  = ((t1-t0)*mhz)/(1*(Loop/1000000)*N);
	fcycles = ((t1-t0)*mhz)/(1.0*(Loop/1000000)*N);
	fprintf(stderr,"float sin used %d cycles per instruction %f\n", (int)cycles, fcycles);


	return 0;
}

