#include "genfft.h"

/**
*   NAME:     crmfft
*
*   DESCRIPTION: Multiple vector real to complex FFT
*
*   USAGE:
*	      void crmfft(complex *cdata, float *rdata, int n1, int n2, 
*                     int ldc, int ldr, int sign)
*
*   INPUT:  - *cdata: complex 2D input array 
*           -     n1: number of (real) samples to be transformed
*           -     n2: number of vectors to be transformed
*           -    ldc: leading dimension (number of complex samples)
*           -    ldr: leading dimension (number of real samples)
*           -   sign: sign of the Fourier kernel 
*
*   OUTPUT: - *rdata: real 2D output array unscaled 
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
*    2.1       Alexander Koek   Feb. '98    updated complib version
*    2.2       Alexander Koek   June '98    updated SCS for use inside
*                                           parallel loops
*
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*	        The Netherlands
*
----------------------------------------------------------------------*/

void crmfft(complex *cdata, float *rdata, int n1, int n2, int ldc, int ldr, int sign)
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
	int    ntable, nwork, zero=0;
	int    ld1, j, nmp;
	static int isys;
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS], scale=1.0;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	static int nprev=0;
	int j;
	static float *coeff;
	float  *ctmp;
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
		if (work == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		GHFFT(&zero, &n1, &scale, cdata, rdata, table, work, &isys);
		nprev = n1;
	}
	
	for (j=0; j<n2; j++ ) {
		GHFFT(&sign, &n1, &scale, &cdata[j*ldc], &rdata[j*ldr], table, 
			work, &isys);
	}
#elif defined(CRAY_MPP_64)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1;
		nwork  = 2*n1;
		if (work) free(work);
		work = (double *)malloc(nwork*sizeof(double));
		if (work == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		if (table) free(table);
		table = (double *)malloc(ntable*sizeof(double));
		if (table == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		if (ddata) free(ddata);
		ddata = (double *)malloc(4*ldc*sizeof(double));
		if (ddata == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		CSFFT(&zero, &n1, &scale, ddata, ddata, table, work, &isys);
		nprev = n1;
	}
	for (j=0; j<n2; j++) {
		for (i=0; i<ldc; i++) {
			ddata[2*i] = (double) cdata[j*ldc+i].r;
			ddata[2*i+1] = (double) cdata[j*ldc+i].i;
		}
		CSFFT(&sign, &n1, &scale, ddata, ddata, table, work, &isys);
		for (i=0; i<n1; i++) {
			rdata[j*ldr+i] = (float) ddata[i];
		}
	}
#elif defined(CRAY_PVP)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1 + 100;
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		CSFFTM(&zero, &n1, &n2, &scale, cdata, &ldc, rdata, &ldr, table, work, &isys);
		nprev  = n1;
	}
	if(n2 != n2prev || n1 != nprev) {
		nwork  = (2*n1+4)*n2;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		n2prev = n2;
	}
	if (!(ldr & 01)) { /* leading dimension is not odd */
		ld0 = 2*ldc+1;
		tmp = (float *) malloc(ld0*n2*sizeof(float));
		if (tmp == NULL) fprintf(stderr,"crmfft: memory allocation error\n");
		CSFFTM(&sign, &n1, &n2, &scale, cdata, &ldc, tmp, &ld0, table, work, &isys);
		for (j=0; j<n2; j++) {
			memcpy((float *)&rdata[j*ldr], (float *)&tmp[j*ld0], sizeof(float)*n1);
		}
		free(tmp);
	}
	else {	
		CSFFTM(&sign, &n1, &n2, &scale, cdata, &ldc, rdata, &ldr, table, work, &isys);
	}
#elif defined(HAVE_LIBSCS)
	nmp = mp_my_threadnum();
	if(nmp>=MAX_NUMTHREADS) {
		fprintf(stderr,"crmfft: cannot handle more than %d processors\n",
			MAX_NUMTHREADS);
		exit(1);
	}

	if (n1 != nprev[nmp]) {
		isys   = 0;
		ntable = n1 + 15;
		nwork  = n1+1;
		if (work[nmp]) free(work[nmp]);
		work[nmp] = (float *)malloc(nwork*sizeof(float));
		if (work[nmp] == NULL) {
		 	fprintf(stderr,"crmfft: memory allocation error in work[%d]\n",
				nmp);
			exit(1);
		}

		if (table[nmp]) free(table[nmp]);
		table[nmp] = (float *)malloc(ntable*sizeof(float));
		if (table[nmp] == NULL) {
		 	fprintf(stderr,"crmfft: memory allocation error in table[%d]\n",
				nmp);
			exit(1);
		}
		csfftm_(&zero, &n1, &n2, &scale, cdata, &ldc, rdata, &ldr, table[nmp], work[nmp], &isys);
		nprev[nmp] = n1;
	}
	ld1 = 2*ldc;
	csfftm_(&sign, &n1, &n2, &scale, cdata, &ldc, cdata, &ld1, table[nmp], work[nmp], &isys);

	for (j = 0; j < n2; j++) {
		memcpy((float *)&rdata[j*ldr], (float *)&cdata[j*ldc], sizeof(float)*n1);
	}
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	if (n1 != nprev) {
		if (coeff) free(coeff);
		coeff = (float *)scfftm1dui(n1, NULL);
		nprev = n1;
	}

	ctmp = (float *)malloc( 2*ldc*n2*sizeof(float) );
	memcpy( ctmp, &cdata[0], 2*ldc*n2*sizeof(float));

	csfftm1du(sign, n1, n2, ctmp, 1, (2*ldc), coeff);

	for (j = 0; j < n2; j++) {
		memcpy(&rdata[j*ldr], (float *)&ctmp[j*2*ldc], sizeof(float)*n1);
	}
	free( ctmp );
#else
	crm_fft(cdata, rdata, n1, n2, ldc, ldr, sign);
#endif

	return;
}


/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncrmfft	FNAME(CRMFFTF)
#else
#define ncrmfft	FNAME(crmfftf)
#endif

void ncrmfft(complex *cdata, float *rdata, int *n1, int *n2, int *ldc, int *ldr, int *sign)
{
	crmfft(cdata, rdata, *n1, *n2, *ldc, *ldr, *sign);

	return;
}

