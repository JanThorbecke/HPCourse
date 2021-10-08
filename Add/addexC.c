#include <math.h>
#include<stdio.h>
#include <time.h>
#include <sys/time.h>

void addC(const double a[], const double b[], double c[], double f, int n);
double wallclock_time(void);

int main(int argc, char *argv[])
{	
	double a[1024], b[1024], c[1024], f;
	double t0, t1;
	int n, i;

	n = 1024;
	f = 1.0;
	for (i = 0; i < n; i++) {
		a[i] = i;
		b[i] = i;
		c[i] = 1.0;
	}
	t0=wallclock_time();
	for (i = 0; i < n; i++) {
		addC(a, b, c, f, n);
	}
	t1=wallclock_time();
	fprintf(stderr,"C-function %e seconds\n", t1-t0);
	fprintf(stderr,"C[0]= %f c[1] = %f\n",c[0], c[1]);

	for (i = 0; i < n; i++) {
		c[i] = 1.0;
	}

	addC(a, &c[0], &c[1], f, n-1);
	fprintf(stderr,"C[0]= %f c[1] = %f\n",c[0], c[1]);
	return 0;
}

void addC(const double a[], const double b[], double c[], double f, int n)
{
	int k;
	for(k = 0; k < n; k++) c[k] = a[k] + f * b[k];
}

