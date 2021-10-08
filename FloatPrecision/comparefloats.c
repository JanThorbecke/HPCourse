#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char *argv[]) 
{
	float a, b, c, d;
	int i, j;

	a=0.1;
	b=0.0;
	fprintf(stderr,"a=%e b=%e\n",a, b);

	for (i=0; i<10; i++) {
		b = b + a;
	}
	if ( b == 1.0 ) {
		fprintf(stderr,"Summation is correct\n");
	}
	else {
		fprintf(stderr,"What is wrong here?\n");
	}
	fprintf(stderr,"difference b-1.0=%g\n",b-1.0);

/* adding small and big */

	a = 1.2000e8;
	b = a + 1.123e-6;
	fprintf(stderr," b=%lf\n",b);

/* adding numbers in order */

	a = 1.200000;
	b = 0.0000001;
	c = 0.00000006;
	d = 0.0;
	for (i=0; i<10; i++) {
		d = d+(a+b)+c;
	}
	fprintf(stderr," d=%lf\n",d);
	d = 0.0;
	for (i=0; i<10; i++) {
		d = d+a+(b+c);
	}
	fprintf(stderr," d=%lf\n",d);

/* inversion of numbers */

	j = 0;
	for (d=1.0; d<1000.0; d+=1.0) {
		a = 1.0 / d;
		b = a * d;
		if ( b != 1.0 ) {
			j++;
		}
	}
	fprintf(stderr," out of 1000 found %d problems\n",j);

	return 0;
}

