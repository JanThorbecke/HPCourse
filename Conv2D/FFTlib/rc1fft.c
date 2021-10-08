#include "genfft.h"

/**
*   NAME:     rc1fft
*
*   DESCRIPTION: real to complex FFT
*
*   USAGE:
*	      void rc1fft(float *rdata, complex *cdata, int n, int sign)
*
*   INPUT:  - *rdata: real 1D input vector 
*           -      n: number of (real) samples in input vector data
*           -   sign: sign of the Fourier kernel 
*
*   OUTPUT: - *cdata: complex 1D output vector unscaled 
*
*   Notice in the preceding formula that there are n real input values,
*     and n/2 + 1 complex output values.  This property is characteristic of
*     real-to-complex FFTs.
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - CRAY T3D and T3E
*          - CRAY T90
*          - CRAY J90
*          - SGI/CRAY ORIGIN 2000 (scsl)
*          - SGI Power Challenge (complib)
*          - inplace FFT from Mayer and SU (see file fft_mayer.c)
*
*   AUTHOR:
*	        Jan Thorbecke (jant@demeern.sgi.com)
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*	        The Netherlands
*
*
*----------------------------------------------------------------------
*  REVISION HISTORY:
*  VERSION        AUTHOR          DATE         COMMENT
*    1.0       Jan Thorbecke    Feb  '94    Initial version (TU Delft)
*    1.1       Jan Thorbecke    June '94    faster in-place FFT 
*    2.0       Jan Thorbecke    July '97    added Cray SGI calls 
*
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*	        The Netherlands
*
----------------------------------------------------------------------*/

void rc1fft(float *rdata, complex *cdata, int n, int sign)
{
#if defined(CRAY_MPP)
	static int nprev=0;
	int   ntable, nwork, zero=0;
	static int isys;
	static float *work, *table, scale=1.0;
#elif defined(CRAY_MPP_64)
	static int nprev=0;
	int   ntable, nwork, zero=0, i;
	static int isys;
	static double *ddata, *work, *table, scale=1.0;
#elif defined(CRAY_PVP)
	static int nprev=0;
	int   ntable, nwork, zero=0;
	static int isys;
	static float *work, *table, scale=1.0;
#elif defined(HAVE_LIBSCS)
	static int nprev=0;
	int   ntable, nwork, zero=0;
	static int isys;
	static float *work, *table, scale=1.0;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	static int nprev=0;
	float *r, *c;
	int j;
	static float *coeff;
#elif defined(CONVEX)
	static int nprev=0;
	int   iopt, ier;
	static float *work;
#endif

#if defined(CRAY_MPP)
	if (n != nprev) {
		isys   = 0;
		ntable = 2*n;
		nwork  = 2*n;

		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		HGFFT(&zero, &n, &scale, rdata, cdata, table, work, &isys);
		nprev = n;
	}
	
	HGFFT(&sign, &n, &scale, rdata, cdata, table, work, &isys);
#elif defined(CRAY_MPP_64)
	if (n != nprev) {
		isys   = 0;
		ntable = 2*n;
		nwork  = 2*n;

		if (work) free(work);
		work = (double *)malloc(nwork*sizeof(double));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		if (table) free(table);
		table = (double *)malloc(ntable*sizeof(double));
		if (table == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		if (ddata) free(ddata);
		ddata = (double *)malloc((n+2)*sizeof(double));
		if (ddata == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		SCFFT(&zero, &n, &scale, ddata, ddata, table, work, &isys);
		nprev = n;
	}
	for (i=0; i<n; i++) {
		ddata[i] = (double) rdata[i];
	}
	SCFFT(&sign, &n, &scale, ddata, ddata, table, work, &isys);
	for (i=0; i<n/2+1; i++) {
		cdata[i].r = (float) ddata[2*i];
		cdata[i].i = (float) ddata[2*i+1];
	}
#elif defined(CRAY_PVP)
	if (n != nprev) {
		isys   = 0;
		ntable = 4*n + 100;
		nwork  = 4*n + 4;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		SCFFT(&zero, &n, &scale, rdata, cdata, table, work, &isys);
		nprev = n;
	}
	SCFFT(&sign, &n, &scale, rdata, cdata, table, work, &isys);
#elif defined(HAVE_LIBSCS)
	if (n != nprev) {
		isys   = 0;
		ntable = n + 15;
		nwork  = n + 2;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		scfft_(&zero, &n, &scale, rdata, cdata, table, work, &isys);
		nprev = n;
	}
	scfft_(&sign, &n, &scale, rdata, cdata, table, work, &isys);
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	if (n != nprev) {
		if (coeff) free(coeff);
		coeff = (float *)scfft1dui(n, NULL);
		nprev = n;
	}

	r = (float *)&rdata[0];
	c = (float *)&cdata[0];
	for (j = 0; j < n; j++) *c++ = *r++;
	*c++ = 0.0; *c++ = 0.0;

	scfft1du(sign, n, (float *)&cdata[0], 1, coeff);
#else
	rc1_fft(rdata, cdata, n, sign);
#endif

	return;
}


/****************** NO COMPLEX DEFINED ******************/

void Rrc1fft(float *rdata, float *cdata, int n, int sign)
{
	rc1fft(rdata, (complex *)cdata, n , sign);
	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nrc1fft	FNAME(RC1FFTF)
#else
#define nrc1fft	FNAME(rc1fftf)
#endif

void nrc1fft(float *rdata, complex *cdata, int *n, int *sign)
{
	rc1fft(rdata, cdata, *n, *sign);

	return;
}

