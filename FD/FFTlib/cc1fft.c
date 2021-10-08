#include "genfft.h"

/**
*   NAME:     cc1fft
*
*   DESCRIPTION: complex to complex FFT
*
*   USAGE:
*	      void cc1fft(complex *data, int n, int sign)
*
*   INPUT:  - *data: complex 1D input vector 
*           -     n: number of samples in input vector data
*           -  sign: sign of the Fourier kernel 
*
*   OUTPUT: - *data: complex 1D output vector unscaled 
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - CONVEX
*          - CRAY T3D and T3E
*          - CRAY T90
*          - CRAY J90
*          - SGI/CRAY ORIGIN 2000 (scsl)
*          - inplace FFT from Mayer and SU (see file fft_mayer.c)
*
*   AUTHOR:
*           Jan Thorbecke (jant@demeern.sgi.com)
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*           The Netherlands
*
*
*----------------------------------------------------------------------
*  REVISION HISTORY:
*  VERSION        AUTHOR          DATE         COMMENT
*    1.0       Jan Thorbecke    Feb  '94    Initial version (TU Delft)
*    1.1       Jan Thorbecke    June '94    faster in-place FFT 
*    2.0       Jan Thorbecke    July '97    added Cray SGI calls 
*    2.1       Alexander Koek   June '98    updated SCS for use inside
*                                           parallel loops
*
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*           The Netherlands
*
----------------------------------------------------------------------*/

void cc1fft(complex *data, int n, int sign)
{
#if defined(CRAY_MPP)
	int   ntable, nwork, zero=0;
	static int isys, nprev=0;
	static float *work, *table, scale=1.0;
#elif defined(CRAY_MPP_64)
	int   ntable, nwork, zero=0, i;
	static int isys, nprev=0;
	static double *ddata, *work, *table, scale=1.0;
#elif defined(CRAY_PVP)
	int   ntable, nwork, zero=0;
	static int isys, nprev=0;
	static float *work, *table, scale=1.0;
#elif defined(HAVE_LIBSCS)
	int    ntable, nwork, zero=0;
	static int isys, nprev[MAX_NUMTHREADS];
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS], scale=1.0;
	int    pe, i;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	static int nprev[MAX_NUMTHREADS];
	static complex *coeff[MAX_NUMTHREADS];
	int    pe;
#elif defined(CONVEX)
	static int nprev=0;
	int   iopt, ier;
	static float *work;
#endif

#if defined(CRAY_MPP)
	if (n != nprev) {
		isys   = 0;
		ntable = 2*n;
		nwork  = 4*n;

		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");
		GGFFT(&zero, &n, &scale, data, data, table, work, &isys);
		nprev = n;
	}
	
	GGFFT(&sign, &n, &scale, data, data, table, work, &isys);
#elif defined(CRAY_MPP_64)
	if (n != nprev) {
		if (factorized(n)) {
			isys   = 0;
			ntable = 2*n;
			nwork  = 4*n;
		}
		else {
			isys   = 1;
			ntable = 12*n;
			nwork  = 8*n;
		}
		if (work) free(work);
		work = (double *)malloc(nwork*sizeof(double));
		if (work == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");
		if (table) free(table);
		table = (double *)malloc(ntable*sizeof(double));
		if (table == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");
		if (ddata) free(ddata);
		ddata = (double *)malloc(2*n*sizeof(double));
		if (ddata == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");
		CCFFT(&zero, &n, &scale, ddata, ddata, table, work, &isys);
		nprev = n;
	}
	for (i=0; i<n; i++) {
		ddata[2*i] = (double) data[i].r;
		ddata[2*i+1] = (double) data[i].i;
	}
	CCFFT(&sign, &n, &scale, ddata, ddata, table, work, &isys);
	for (i=0; i<n; i++) {
		data[i].r = (float) ddata[2*i];
		data[i].i = (float) ddata[2*i+1];
	}
#elif defined(CRAY_PVP)
	if (n != nprev) {
		isys   = 0;
		ntable = 8*n + 100;
		nwork  = 8*n;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");
		CCFFT(&zero, &n, &scale, data, data, table, work, &isys);
		nprev = n;
	}
	CCFFT(&sign, &n, &scale, data, data, table, work, &isys);
#elif defined(HAVE_LIBSCS)
	pe = mp_my_threadnum();
	assert ( pe <= MAX_NUMTHREADS );
	if (n != nprev[pe]) {
		isys   = 0;
		ntable = 2*n + 30;
		nwork  = 2*n;
		/* allocate memory on each processor locally for speed */
		if (work[pe]) free(work[pe]);
		work[pe] = (float *)malloc(nwork*sizeof(float));
		if (work[pe] == NULL) 
			fprintf(stderr,"cc1fft: memory allocation error\n");
		if (table[pe]) free(table[pe]);
		table[pe] = (float *)malloc(ntable*sizeof(float));
		if (table[pe] == NULL) 
			fprintf(stderr,"cc1fft: memory allocation error\n");
		ccfft_(&zero, &n, &scale, data, data, table[pe], work[pe], &isys);
		nprev[pe] = n;
	}
	ccfft_(&sign, &n, &scale, data, data, table[pe], work[pe], &isys);
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	pe = mp_my_threadnum();
	assert ( pe <= MAX_NUMTHREADS );
	if (n != nprev[pe]) {
		if (coeff[pe]) free(coeff[pe]);
		coeff[pe] = (complex *)cfft1di(n, NULL);
		nprev[pe] = n;
	}
	cfft1d(sign, n, (complex *) data, 1, coeff[pe]);
#elif defined(CONVEX)
	if (n != nprev) {
		if (work) free(work);
		work = (float *)malloc(5*n/2*sizeof(float));
		if (work == NULL) fprintf(stderr,"cc1fft: memory allocation error\n");

		iopt = -3;
		c1dfft_(data, &n, work, &iopt, &ier);
		if (ier != 0) fprintf(stderr,"cc1fft: Error in fft ier = %d\n", ier);
		nprev = n;
	}

	if (sign > 0) iopt = -2;
	else iopt = 1;
	c1dfft_(data, &n, work, &iopt, &ier);
	if (ier != 0) fprintf(stderr,"cc1fft: Error in fft ier = %d\n", ier);
#else
	cc1_fft(data, n, sign);
#endif

	return;
}

/****************** NO COMPLEX DEFINED ******************/

void Rcc1fft(float *data, int n, int sign)
{
    cc1fft((complex *)data, n , sign);
    return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncc1fft	FNAME(CC1FFTF)
#else
#define ncc1fft	FNAME(cc1fftf)
#endif

void ncc1fft(complex *data, int *n, int *sign)
{
	cc1fft(data, *n, *sign);

	return;
}

