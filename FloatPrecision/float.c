#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char *argv[]) 
{
	float a, b;
	double c, d;
	int i, j, k;

	a=1.0;
	d=1.0;
	c=1.0;
	fprintf(stderr,"a=%f c=%g\n",a, c);

	for (i=0; i<10; i++) {
		c = 0.1*c;
		a = a + c;
		d = d + c;
		fprintf(stderr,"i=%d a=%16.14f c=%10.14f d=%16.14f \n",i, a, c, d);
	}

	return 0;
}

