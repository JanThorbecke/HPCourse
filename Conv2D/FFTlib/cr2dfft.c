#include "genfft.h"

/**
*   NAME:     cr2dfft
*
*   DESCRIPTION: 2 Dimensional complex to real FFT
*
*   USAGE:
*         void cr2dfft(complex *cdata, float *rdata, int nr, int nc, 
*                      int ldc, int ldr, int sign)
*
*   INPUT:  - *cdata: complex 2D input array [nc][nr]
*           -     nr: number of real (fast) samples to be transformed
*           -     nc: number of complex (slow) samples to be transformed
*           -    ldc: leading dimension (number of complex samples)
*           -    ldr: leading dimension (number of real samples)
*           -   sign: sign of the Fourier kernel
*
*   OUTPUT: - *rdata: real 2D output array unscaled [nc][nr]
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

void cr2dfft(complex *cdata, float *rdata, int nr, int nc, int ldc, int ldr, int sign)
{
#if defined(CRAY_MPP)
	static int nr_prev=0, nc_prev=0;
	int nwork, ntable, i, j, nf;
	static float scale=1.0;
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
	int ld1, j;
	static float *work, *table;
	float scale=1.0;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	static int nr_prev=0, nc_prev=0;
	int j;
	static float *coeff;
	float *ctmp;
#else
	int i, j, nf;
	complex *tmp;
#endif

#if defined(CRAY_MPP)

	if (nc != nc_prev) {
		nwork = 4*nc;
		ntable = 2*nc;
		if (work2) free (work2);
		work2  = (float *)malloc(nwork*sizeof(float));
		if (work2 == NULL) fprintf(stderr,"cr2dfft: memory allocation error\n");
		if (table2) free (table2);
		table2 = (float *)malloc(ntable*sizeof(float));
		if (table2 == NULL) fprintf(stderr,"cr2dfft: memory allocation error\n");
		if (tmp) free (tmp);
		tmp = (complex *)malloc(nc*sizeof(complex));
		if (tmp == NULL) fprintf(stderr,"cr2dfft: memory allocation error\n");
		GGFFT(&zero,&nc,&scale,tmp,tmp,table2,work2,&isys);
		nc_prev = nc;
	}

	nf = (nr+2)/2;
	for (i=0; i<nf; i++) {
		for (j=0; j<nc; j++) tmp[j] = cdata[j*ldc+i];
	 	GGFFT(&sign,&nc,&scale,tmp,tmp,table2,work2,&isys);
		for (j=0; j<nc; j++) cdata[j*ldc+i] = tmp[j];
	}

	if (nr != nr_prev) {
		nwork = 2*nr;
		ntable = 2*nr;
		if (work) free (work);
		work  = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"cr2dfft: memory allocation error\n");
		if (table) free (table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"cr2dfft: memory allocation error\n");
		GHFFT(&isys,&nr,&scale,cdata,rdata,table,work,&isys);
		nr_prev = nr;
	}
	for (j=0; j<nc; j++) {
		GHFFT(&sign,&nr,&scale,&cdata[j*ldc],&rdata[j*ldr],table,work,&isys);
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
		if (ddata == NULL) fprintf(stderr,"cr2dfft: memory allocation error\n");
		CSFFT2D(&zero,&nr,&nc,&scale,ddata,&lddc,ddata,&ldd,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}
	for (j=0; j<nc; j++) {
		for (i=0; i<lddc; i++) {
			ddata[j*2*lddc+2*i] = (double) cdata[j*ldc+i].r;
			ddata[j*2*lddc+2*i+1] = (double) cdata[j*ldc+i].i;
		}
	}
	CSFFT2D(&sign,&nr,&nc,&scale,ddata,&lddc,ddata,&ldd,table,work,&isys);
	for (j=0; j<nc; j++) {
		for (i=0; i<nr; i++) {
			rdata[j*ldr+i] = (float) ddata[j*ldd+i];
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
		CSFFT2D(&zero,&nr,&nc,&scale,cdata,&ldc,rdata,&ldr,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}

	ld1 = 2*ldc;
	CSFFT2D(&sign,&nr,&nc,&scale,cdata,&ldc,cdata,&ld1,table,work,&isys);

	for (j=0; j<nc; j++) {
		memcpy((float *)&rdata[j*ldr], (float *)&cdata[j*ldc].r, sizeof(float)*nr);
	}

#elif defined(HAVE_LIBSCS)
	if (nr != nr_prev || nc != nc_prev) {
		nwork = nr*nc;
		ntable = 15+nr+2*(15+nc);
		if (work) free (work);
		if (table) free (table);
		work  = (float *)malloc(nwork*sizeof(float));
		table = (float *)malloc(ntable*sizeof(float));
		csfft2d_(&zero,&nr,&nc,&scale,cdata,&ldc,rdata,&ldr,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}
	ld1 = 2*ldc;
	csfft2d_(sign,&nr,&nc,&scale,cdata,&ldc,cdata,&ld1,table,work,&isys);

	for (j=0; j<nc; j++) {
		memcpy((float *)&rdata[j*ldr], (float *)&cdata[j*ldc].r, sizeof(float)*nr);
	}

#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	if (nr != nr_prev || nc != nc_prev) {
		if (coeff) free(coeff);
		coeff = (float *)scfft2dui(nr, nc, NULL);
		nr_prev = nr;
		nc_prev = nc;
	}

	/* make copy to avoid overwriting the input */
	ctmp = (float *)malloc(2*ldc*nc*sizeof(float));
	memcpy( ctmp, (float *)&cdata[0], 2*ldc*nc*sizeof(float) );

	csfft2du(sign, nr, nc, (float *)&ctmp[0], 2*ldc, coeff);

	for (j=0; j<nc; j++) {
		memcpy((float *)&rdata[j*ldr], (float *)&ctmp[j*2*ldc], sizeof(float)*nr);
	}

	free(ctmp);
#else 
	tmp = (complex *)malloc(nc*sizeof(complex));
	nf = (nr+2)/2;
	for (i=0; i<nf; i++) {
		for (j=0; j<nc; j++) tmp[j] = cdata[j*ldc+i];
		cc1fft(tmp, nc, sign);
		for (j=0; j<nc; j++) cdata[j*ldc+i] = tmp[j];
	}
	free (tmp);
	crmfft(cdata, rdata, nr, nc, ldc, ldr, sign);
#endif

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncr2dfft FNAME(CR2DFFTF)
#else
#define ncr2dfft FNAME(cr2dfftf)
#endif

void ncr2dfft(complex *cdata, float *rdata, int *nr, int *nc, int *ldc, int *ldr, int *sign)
{
	cr2dfft(cdata, rdata, *nr, *nc, *ldc, *ldr, *sign);

	return;
}

