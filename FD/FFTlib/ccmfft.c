#include "genfft.h"

/**
*   NAME:     ccmfft
*
*   DESCRIPTION: Multipple vector complex to complex FFT
*
*   USAGE:
*	      void ccmfft(complex *data, int n1, int n2, int ld1, int sign)
*
*   INPUT:  - *data: complex 2D input array 
*           -    n1: number of samples to be transformed
*           -    n2: number of vectors to be transformed
*           -   ld1: leading dimension
*           -  sign: sign of the Fourier kernel 
*
*   OUTPUT: - *data: complex 2D output array unscaled 
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - CRAY T3D and T3E
*          - CRAY T90
*          - CRAY J90
*          - SGI/CRAY ORIGIN 2000 (scsl)
*          - SGI Power Challenge (complib.sgimath)
*          - CONVEX
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

void ccmfft(complex *data, int n1, int n2, int ld1, int sign)
{
#if defined(CRAY_MPP)
	int   ntable, nwork, zero=0, j;
	static int isys, nprev=0;
	static float *work, *table, scale=1.0;
#elif defined(CRAY_MPP_64)
	int   ntable, nwork, zero=0, i, j;
	static int isys, nprev=0;
	static double *ddata, *work, *table, scale=1.0;
#elif defined(CRAY_PVP)
	int   ntable, nwork, zero=0, i, j, ld0;
	complex *tmp;
	static int isys, n2prev=0, nprev=0;
	static float *work, *table, scale=1.0;
#elif defined(HAVE_LIBSCS)
	int   ntable, nwork, zero=0;
	static int isys, nprev=0;
	static float *work, *table, scale=1.0;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	static int nprev=0;
	static complex *coeff;
#elif defined(CONVEX)
	static int nprev=0;
	int   iopt, ier;
	static float *work;
#endif

#if defined(CRAY_MPP)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1;
		nwork  = 4*n1;

		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		GGFFT(&zero, &n1, &scale, data, data, table, work, &isys);
		nprev = n1;
	}
	
	for (j=0; j<n2; j++ ) {
		GGFFT(&sign, &n1, &scale, &data[j*ld1], &data[j*ld1], table, 
			work, &isys);
	}
#elif defined(CRAY_MPP_64)
	if (n1 != nprev) {
		if (factorized(n1)) {
			isys   = 0;
			ntable = 2*n1;
			nwork  = 4*n1;
		}
		else {
			isys   = 1;
			ntable = 12*n1;
			nwork  = 8*n1;
		}
		if (work) free(work);
		work = (double *)malloc(nwork*sizeof(double));
		if (work == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		if (table) free(table);
		table = (double *)malloc(ntable*sizeof(double));
		if (table == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		if (ddata) free(ddata);
		ddata = (double *)malloc(2*n1*sizeof(double));
		if (ddata == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		CCFFT(&zero, &n1, &scale, ddata, ddata, table, work, &isys);
		nprev = n1;
	}
	for (j=0; j<n2; j++) {
		for (i=0; i<n1; i++) {
			ddata[2*i] = (double) data[j*ld1+i].r;
			ddata[2*i+1] = (double) data[j*ld1+i].i;
		}
		CCFFT(&sign, &n1, &scale, ddata, ddata, table, work, &isys);
		for (i=0; i<n1; i++) {
			data[j*ld1+i].r = (float) ddata[2*i];
			data[j*ld1+i].i = (float) ddata[2*i+1];
		}
	}
#elif defined(CRAY_PVP)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1 + 100;
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		CCFFTM(&zero, &n1, &n2, &scale, data, &ld1, data, &ld1, table, work, &isys);
		nprev  = n1;
	}
	if(n2 != n2prev || n1 != nprev) {
		nwork  = 4*n1*n2;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		n2prev = n2;
	}
	if (!(ld1 & 01)) { /* leading dimension is not odd */
		ld0 = n1+1;
		tmp = (complex *) malloc(ld0*n2*sizeof(complex));
		if (tmp == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		for (j=0; j<n2; j++) {
			for (i=0; i<n1; i++) tmp[j*ld0+i] = data[j*ld1+i];
		}
		CCFFTM(&sign, &n1, &n2, &scale, tmp, &ld0, tmp, &ld0, table, work, &isys);
		for (j=0; j<n2; j++) {
			for (i=0; i<n1; i++) data[j*ld1+i] = tmp[j*ld0+i];
		}
		free(tmp);
	}
	else {	
		CCFFTM(&sign, &n1, &n2, &scale, data, &ld1, data, &ld1, table, work, &isys);
	}
#elif defined(HAVE_LIBSCS)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1 + 30;
		nwork  = 2*n1;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		ccfftm_(&zero, &n1, &n2, &scale, data, &ld1, data, &ld1, table, work, &isys);
		nprev = n1;
	}
	ccfftm_(&sign, &n1, &n2, &scale, data, &ld1, data, &ld1, table, work, &isys);
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	if (n1 != nprev) {
		if (coeff) free(coeff);
		coeff = (complex *)cfftm1di(n1, NULL);
		nprev = n1;
	}
	cfftm1d(sign, n1, n2, data, 1, ld1, coeff);
#elif defined(CONVEX)
	if (n1 != nprev) {
		if (work) free(work);
		work = (float *)malloc(5*n1/2*sizeof(float));
		if (work == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");

		iopt = -3;
		c1dfft_(data, &n1, work, &iopt, &ier);
		if (ier != 0) fprintf(stderr,"ccmfft: Error in fft ier = %d\n", ier);
		nprev = n1;
	}

	if (sign > 0) iopt = -2;
	else iopt = 1;
	for (j=0; j<n2; j++ {
		c1dfft_(data[j*ld1], &n1, work, &iopt, &ier);
		if (ier != 0) fprintf(stderr,"ccmfft: Error in fft ier = %d\n", ier);
	}
#else
	ccm_fft(data, n1, n2, ld1, sign);
#endif

	return;
}


/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nccmfft	FNAME(CCMFFTF)
#else
#define nccmfft	FNAME(ccmfftf)
#endif

void nccmfft(complex *data, int *n1, int *n2, int *ld1, int *sign)
{
	ccmfft(data, *n1, *n2, *ld1, *sign);

	return;
}

