#include "genfft.h"

/**
*   NAME:     rcmfft
*
*   DESCRIPTION: Multiple vector real to complex FFT
*
*   USAGE:
*	      void rcmfft(float *rdata, complex *cdata, int n1, int n2, 
*                     int ldr, int ldc, int sign)
*
*   INPUT:  - *rdata: real 2D input array 
*           -     n1: number of (real) samples to be transformed
*           -     n2: number of vectors to be transformed
*           -    ldr: leading dimension (number of real samples)
*           -    ldc: leading dimension (number of complex samples)
*           -   sign: sign of the Fourier kernel 
*
*   OUTPUT: - *cdata: complex 2D output array unscaled 
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - CRAY T3D and T3E
*          - CRAY T90
*          - CRAY J90
*          - SGI/CRAY ORIGIN 2000 (scsl)
*          - SGI Power Challenge (complib.sgimath)
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
*    2.1       Alexander Koek   June '98    updated SCS for use inside
*                                           parallel loops
*
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*	        The Netherlands
*
----------------------------------------------------------------------*/

void rcmfft(float *rdata, complex *cdata, int n1, int n2, int ldr, int ldc, int sign)
{
#if defined(CRAY_MPP)
	static int nprev=0;
	int   ntable, nwork, zero=0, j;
	static int isys;
	static float *work, *table, scale=1.0;
#elif defined(CRAY_MPP_64)
	static int nprev=0;
	int   ntable, nwork, zero=0, i, j;
	static int isys;
	static double *ddata, *work, *table, scale=1.0;
#elif defined(CRAY_PVP)
	static int nprev=0;
	int   ntable, nwork, zero=0, i, j, ld0;
	float *tmp;
	static int isys, n2prev=0;
	static float *work, *table, scale=1.0;
#elif defined(HAVE_LIBSCS)
	static int nprev[MAX_NUMTHREADS];
	int    nmp, ntable, nwork, zero=0;
	static int isys;
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS], scale=1.0;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	static int nprev=0;
	int j;
	static float *coeff;
#elif defined(CONVEX)
	static int nprev=0;
	int   iopt, ier;
	static float *work;
#endif

#if defined(CRAY_MPP)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1;
		nwork  = 2*n1;

		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		HGFFT(&zero, &n1, &scale, rdata, cdata, table, work, &isys);
		nprev = n1;
	}
	
	for (j=0; j<n2; j++ ) {
		HGFFT(&sign, &n1, &scale, &rdata[j*ldr], &cdata[j*ldc], table, 
			work, &isys);
	}
#elif defined(CRAY_MPP_64)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1;
		nwork  = 2*n1;
		if (work) free(work);
		work = (double *)malloc(nwork*sizeof(double));
		if (work == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		if (table) free(table);
		table = (double *)malloc(ntable*sizeof(double));
		if (table == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		if (ddata) free(ddata);
		ddata = (double *)malloc(2*n1*sizeof(double));
		if (ddata == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		SCFFT(&zero, &n1, &scale, ddata, ddata, table, work, &isys);
		nprev = n1;
	}
	for (j=0; j<n2; j++) {
		for (i=0; i<n1; i++) {
			ddata[i] = (double) rdata[j*ldr+i];
		}
		SCFFT(&sign, &n1, &scale, ddata, ddata, table, work, &isys);
		for (i=0; i<(n1+2)/2; i++) {
			cdata[j*ldc+i].r = (float) ddata[2*i];
			cdata[j*ldc+i].i = (float) ddata[2*i+1];
		}
	}
#elif defined(CRAY_PVP)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1 + 100;
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		SCFFTM(&zero, &n1, &n2, &scale, rdata, &ldr, cdata, &ldc, table, work, &isys);
		nprev  = n1;
	}
	if(n2 != n2prev || n1 != nprev) {
		nwork  = (2*n1+4)*n2;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		n2prev = n2;
	}
	if (!(ldr & 01)) { /* leading dimension is not odd */
		ld0 = ldc*2+1;
		tmp = (float *) malloc(ld0*n2*sizeof(float));
		if (tmp == NULL) fprintf(stderr,"rcmfft: memory allocation error\n");
		for (j=0; j<n2; j++) {
			memcpy((float *)&tmp[j*ld0], (float *)&rdata[j*ldr], sizeof(float)*n1);
		}
		SCFFTM(&sign, &n1, &n2, &scale, tmp, &ld0, cdata, &ldc, table, work, &isys);
		free(tmp);
	}
	else {	
		SCFFTM(&sign, &n1, &n2, &scale, rdata, &ldr, cdata, &ldc, table, work, &isys);
	}
#elif defined(HAVE_LIBSCS)
	nmp = mp_my_threadnum();
	if(nmp>=MAX_NUMTHREADS) {
	   fprintf(stderr,"rcmfft: cannot handle more than %d processors\n",MAX_NUMTHREADS);
		exit(1);
	}

	if (n1 != nprev[nmp]) {
		isys   = 0;
		ntable = n1 + 15;
		nwork  = n1+2;
		if (work[nmp]) free(work[nmp]);
		work[nmp] = (float *)malloc(nwork*sizeof(float));
		if (work[nmp] == NULL) {
			fprintf(stderr,"rcmfft: memory allocation error in work[%d]\n",nmp);
			exit(1);
		}

		if (table[nmp]) free(table[nmp]);
		table[nmp] = (float *)malloc(ntable*sizeof(float));
		if (table[nmp] == NULL) {
			fprintf(stderr,"rcmfft: memory allocation error in table[%d]\n",nmp);
			exit(1);
		}
		scfftm_(&zero, &n1, &n2, &scale, rdata, &ldr, cdata, &ldc, table[nmp], work[nmp], &isys);
		nprev[nmp] = n1;
	}
	scfftm_(&sign, &n1, &n2, &scale, rdata, &ldr, cdata, &ldc, table[nmp], work[nmp], &isys);
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	if (n1 != nprev) {
		if (coeff) free(coeff);
		coeff = scfftm1dui(n1, NULL);
		nprev = n1;
	}

	for (j = 0; j < n2; j++) {
		memcpy((float *)&cdata[j*ldc], (float *)&rdata[j*ldr], sizeof(float)*n1);
	}
	scfftm1du(sign, n1, n2, (float *)&cdata[0], 1, (2*ldc), coeff);
#else
	rcm_fft(rdata, cdata, n1, n2, ldr, ldc, sign);
#endif

	return;
}


/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nrcmfft	FNAME(RCMFFTF)
#else
#define nrcmfft	FNAME(rcmfftf)
#endif

void nrcmfft(float *rdata, complex *cdata, int *n1, int *n2, int *ldr, int *ldc, int *sign)
{
	rcmfft(rdata, cdata, *n1, *n2, *ldr, *ldc, *sign);

	return;
}

