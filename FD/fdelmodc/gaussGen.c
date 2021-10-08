#include<math.h>
#include<stdlib.h>

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
//	y2 = x2 * w;

	return (float) y1;
}
