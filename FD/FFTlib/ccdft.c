#include <genfft.h>

void ccdft(complex *cdata, int n, int sign)
{
	int i, j, k; 
	float scl, sumr, sumi;
	complex *tmp;
	static float *csval;
	static int nprev=0;

	if (nprev != n) {
		scl = 2.0*PI/n;
		if (csval) free(csval);
		csval = (float *) malloc(2*n*sizeof(float));
		for (i=0; i<n; i++) {
			csval[2*i] = cos(scl*i);
			csval[2*i+1] = sin(scl*i);
		}
		nprev = n;
	}

	tmp = (complex *) malloc(n*sizeof(complex));

	for (i=0; i<n; i++) {
		sumr = sumi = 0.0;
		for (j=0; j<n; j++) {
			k = 2*((i*j)%n);
			sumr += cdata[j].r*csval[k]-sign*cdata[j].i*csval[k+1];
			sumi += cdata[j].i*csval[k]+sign*cdata[j].r*csval[k+1];
		}
		tmp[i].r = sumr;
		tmp[i].i = sumi;
	}

	for (i=0; i<n; i++) cdata[i] = tmp[i];
	free(tmp);

	return;
}

