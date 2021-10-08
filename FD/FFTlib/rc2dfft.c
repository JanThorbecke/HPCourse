#include "genfft.h"

/**
*   NAME:     rc2dfft
*
*   DESCRIPTION: 2 Dimensional real to complex FFT
*
*   USAGE:
*         void rc2dfft(float *rdata, complex *cdata, int nr, int nc, 
*                      int ldr, int ldc, int sign)
*
*   INPUT:  - *rdata: real 2D input array
*           -     nr: number of real (fast) samples to be transformed
*           -     nc: number of complex (slow) samples to be transformed
*           -    ldr: leading dimension (number of real samples)
*           -    ldc: leading dimension (number of complex samples)
*           -   sign: sign of the Fourier kernel
*
*   OUTPUT: - *data: complex 2D output array unscaled
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
*    1.0       Jan Thorbecke    July  '97    Initial version
*
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*           The Netherlands
*
----------------------------------------------------------------------*/

void rc2dfft(float *rdata, complex *cdata, int nr, int nc, int ldr, int ldc, int sign)
{
#if defined(CRAY_MPP)
	static int nr_prev=0, nc_prev=0;
	int nwork, ntable, i, j, nf;
	float scale=1.0;
	int isys=0, zero=0;
	static float *work, *table;
	static float *work2, *table2;
	static complex *tmp;
#elif defined(CRAY_MPP_64)
	static int nr_prev=0, nc_prev=0;
	static int isys[3], ldd, lddc;
	int ntable, nwork, zero=0, i, j;
	static double *work, *table, scale=1.0, *ddata;
#elif defined(CRAY_PVP)
	static int nr_prev=0, nc_prev=0;
	int  isys=0, ntable, nwork, zero=0, i, j;
	int ld1;
	complex *tmp;
	static float *work, *table;
	float scale=1.0;
#elif defined(HAVE_LIBSCS)
	static int nr_prev=0, nc_prev=0;
	int   isys=0, ntable, nwork, zero=0;
	static float *work, *table;
	float scale=1.0;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	static int nr_prev=0, nc_prev=0;
	int j;
	static float *coeff;
#else
	int i, j, nf;
	complex *tmp;
#endif

#if defined(CRAY_MPP)
	if (nr != nr_prev) {
		nwork = 2*nr;
		ntable = 2*nr;
		if (work) free (work);
		work  = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"rc2dfft: memory allocation error\n");
		if (table) free (table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"rc2dfft: memory allocation error\n");
		HGFFT(&isys,&nr,&scale,rdata,cdata,table,work,&isys);
		nr_prev = nr;
	}
	for (j=0; j<nc; j++) {
		HGFFT(&sign,&nr,&scale,&rdata[j*ldr],&cdata[j*ldc],table,work,&isys);
	}

	if (nc != nc_prev) {
		nwork = 4*nc;
		ntable = 2*nc;
		if (work2) free (work2);
		work2  = (float *)malloc(nwork*sizeof(float));
		if (work2 == NULL) fprintf(stderr,"rc2dfft: memory allocation error\n");
		if (table2) free (table2);
		table2 = (float *)malloc(ntable*sizeof(float));
		if (table2 == NULL) fprintf(stderr,"rc2dfft: memory allocation error\n");
		if (tmp) free (tmp);
		tmp = (complex *)malloc(nc*sizeof(complex));
		if (tmp == NULL) fprintf(stderr,"rc2dfft: memory allocation error\n");
		GGFFT(&zero,&nc,&scale,tmp,tmp,table2,work2,&isys);
		nc_prev = nc;
	}

	nf = (nr+2)/2;
	for (i=0; i<nf; i++) {
		for (j=0; j<nc; j++) tmp[j] = cdata[j*ldc+i];
	 	GGFFT(&sign,&nc,&scale,tmp,tmp,table2,work2,&isys);
		for (j=0; j<nc; j++) cdata[j*ldc+i] = tmp[j];
	}
#elif defined(CRAY_MPP_64)
	if (nr != nr_prev || nc != nc_prev) {
		if (factorized(nr) && factorized(nc)) {
			isys[0] = 2;
			isys[1] = 0;
			isys[2] = 0;
		 	ntable  = 2*(nr+nc);
		}
		else {
			isys[0] = 2;
			isys[1] = 1-factorized(nr);
			isys[2] = 1-factorized(nc);
		 	ntable  = 12*(nr+nc);
		}
		nwork = (nr+nc)*nc;
		if (work) free (work);
		work  = (double *)malloc(nwork*sizeof(double));
		if (table) free (table);
		table = (double *)malloc(ntable*sizeof(double));
		ldd = nr+2;
		lddc = ldd/2;
		if (ddata) free(ddata);
		ddata = (double *)malloc(ldd*nc*sizeof(double));
		if (ddata == NULL) fprintf(stderr,"rc2dfft: memory allocation error\n");
		SCFFT2D(&zero,&nr,&nc,&scale,ddata,&ldd,ddata,&lddc,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}
	for (j=0; j<nc; j++) {
		for (i=0; i<nr; i++) {
			ddata[j*ldd+i] = (double) rdata[j*ldr+i];
		}
	}
	SCFFT2D(&sign,&nr,&nc,&scale,ddata,&ldd,ddata,&lddc,table,work,&isys);
	for (j=0; j<nc; j++) {
		for (i=0; i<lddc; i++) {
			cdata[j*ldc+i].r = (float) ddata[j*2*lddc+2*i];
			cdata[j*ldc+i].i = (float) ddata[j*2*lddc+2*i+1];
		}
	}
#elif defined(CRAY_PVP)
	if (nr != nr_prev || nc != nc_prev) {
		nwork = 512*MAX(nr,nc);
		ntable = 100+2*(nr+nc);
		if (work) free (work);
		work  = (float *)malloc(nwork*sizeof(float));
		if (table) free (table);
		table = (float *)malloc(ntable*sizeof(float));
		SCFFT2D(&zero,&nr,&nc,&scale,rdata,&ldr,cdata,&ldc,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}

	SCFFT2D(&sign,&nr,&nc,&scale,rdata,&ldr,cdata,&ldc,table,work,&isys);

#elif defined(HAVE_LIBSCS)
	if (nr != nr_prev || nc != nc_prev) {
		nwork = nr*nc;
		ntable = 15+nr+2*(15+nc);
		if (work) free (work);
		if (table) free (table);
		work  = (float *)malloc(nwork*sizeof(float));
		table = (float *)malloc(ntable*sizeof(float));
		scfft2d_(&zero,&nr,&nc,&scale,rdata,&ldr,cdata,&ldc,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}

	scfft2d_(sign,&nr,&nc,&scale,rdata,&ldr,cdata,&ldc,table,work,&isys);

#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	if (nr != nr_prev || nc != nc_prev) {
		if (coeff) free(coeff);
		coeff = (float *)scfft2dui(nr, nc, NULL);
		nr_prev = nr;
		nc_prev = nc;
	}

	for (j=0; j<nc; j++) {
		memcpy((float *)&cdata[j*ldc], (float *)&rdata[j*ldr], sizeof(float)*nr);
	}

	scfft2du(sign, nr, nc, (float *)&cdata[0], 2*ldc, coeff);

#else 
	rcmfft(rdata, cdata, nr, nc, ldr, ldc, sign);
	tmp = (complex *)malloc(nc*sizeof(complex));
	nf = (nr+2)/2;
	for (i=0; i<nf; i++) {
		for (j=0; j<nc; j++) tmp[j] = cdata[j*ldc+i];
		cc1fft(tmp, nc, sign);
		for (j=0; j<nc; j++) cdata[j*ldc+i] = tmp[j];
	}
	free (tmp);
#endif

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nrc2dfft FNAME(RC2DFFTF)
#else
#define nrc2dfft FNAME(rc2dfftf)
#endif

void nrc2dfft(float *rdata, complex *cdata, int *nr, int *nc, int *ldr, int *ldc, int *sign)
{
	rc2dfft(rdata, cdata, *nr, *nc, *ldr, *ldc, *sign);

	return;
}

