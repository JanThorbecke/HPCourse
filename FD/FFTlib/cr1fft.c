#include "genfft.h"

/**
*   NAME:     cr1fft
*
*   DESCRIPTION: complex to real FFT
*
*   USAGE:
*	      void cr1fft(complex *cdata, float *rdata, int n, int sign)
*
*   INPUT:  - *cdata: complex 1D input vector 
*           -      n: number of (real) samples in input vector data
*           -   sign: sign of the Fourier kernel 
*
*   OUTPUT: - *rdata: real 1D output vector unscaled 
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

void cr1fft(complex *cdata, float *rdata, int n, int sign)
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
	float *r;
	int j;
	static float *coeff;
#elif defined(CONVEX)
	static int nprev=0;
	int   iopt, ier;
	static float *work;
#elif defined(FFTW3)
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
		if (work == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		GHFFT(&zero, &n, &scale, cdata, rdata, table, work, &isys);
		nprev = n;
	}
	
	GHFFT(&sign, &n, &scale, cdata, rdata, table, work, &isys);
#elif defined(CRAY_MPP_64)
	if (n != nprev) {
		isys   = 0;
		ntable = 2*n;
		nwork  = 2*n;

		if (work) free(work);
		work = (double *)malloc(nwork*sizeof(double));
		if (work == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		if (table) free(table);
		table = (double *)malloc(ntable*sizeof(double));
		if (table == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		if (ddata) free(ddata);
		ddata = (double *)malloc((n+2)*sizeof(double));
		if (ddata == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		CSFFT(&zero, &n, &scale, ddata, ddata, table, work, &isys);
		nprev = n;
	}
	for (i=0; i<(n+2)/2; i++) {
		ddata[2*i] = (double) cdata[i].r;
		ddata[2*i+1] = (double) cdata[i].i;
	}
	CSFFT(&sign, &n, &scale, ddata, ddata, table, work, &isys);
	for (i=0; i<n; i++) {
		rdata[i] = (float) ddata[i];
	}
#elif defined(CRAY_PVP)
	if (n != nprev) {
		isys   = 0;
		ntable = 4*n + 100;
		nwork  = 4*n + 4;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		CSFFT(&zero, &n, &scale, cdata, rdata, table, work, &isys);
		nprev = n;
	}
	CSFFT(&sign, &n, &scale, cdata, rdata, table, work, &isys);
#elif defined(HAVE_LIBSCS)
	if (n != nprev) {
		isys   = 0;
		ntable = n + 15;
		nwork  = n+1;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		csfft_(&zero, &n, &scale, cdata, rdata, table, work, &isys);
		nprev = n;
	}
	csfft_(&sign, &n, &scale, cdata, rdata, table, work, &isys);
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	if (n != nprev) {
		if (coeff) free(coeff);
		coeff = (float *)scfft1dui(n, NULL);
		nprev = n;
	}

	r = (float *)calloc( (n+2),sizeof(float) );
	memcpy((float *)&r[0], (float *)&cdata[0], (n+2)*sizeof(float));

	csfft1du(sign, n, r, 1, coeff);

	memcpy((float *)&rdata[0], (float *)&r[0], n*sizeof(float));
	free(r);
#else
	cr1_fft(cdata, rdata, n, sign);
#endif

	return;
}

/****************** NO COMPLEX DEFINED ******************/

void Rcr1fft(float *cdata, float *rdata, int n, int sign)
{
    cr1fft((complex *)cdata, rdata, n, sign);
    return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncr1fft	FNAME(CR1FFTF)
#else
#define ncr1fft	FNAME(cr1fftf)
#endif

void ncr1fft(complex *cdata, float *rdata, int *n, int *sign)
{
	cr1fft(cdata, rdata, *n, *sign);

	return;
}

