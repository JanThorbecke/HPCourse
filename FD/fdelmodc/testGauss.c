#include<stdio.h>
#include<assert.h>
#include<math.h>
#include<stdlib.h>

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

float gaussGen();

int main (int argc, char *argv[]) 
{
	FILE *fp;
	float g1, *gauss;
	int i, n;
	size_t nwrite;

	srand48(10);
	n=1024;
	gauss = (float *)calloc(n,sizeof(float));	
	for (i=0;i<n;i++) {
		g1 = 0.5*n+100*gaussGen();
		gauss[(int)MIN(MAX(0,g1),n-1)] += 1;
	}

	fp = fopen("test.bin", "w");
    nwrite = fwrite( gauss, sizeof(float), n, fp);
    assert(nwrite == n);
	fclose(fp);
	return 0;
}


float gaussGen()
{
	double x1, x2, w, y1, y2;
 
	do {
		x1 = 2.0 * drand48() - 1.0;
		x2 = 2.0 * drand48() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;
	y2 = x2 * w;

	return (float)y1;
}
