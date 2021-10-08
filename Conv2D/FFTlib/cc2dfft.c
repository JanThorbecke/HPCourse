#include "genfft.h"

/**
*   NAME:     cc2dfft
*
*   DESCRIPTION: 2 Dimensional complex to complex FFT
*
*   USAGE:
*         void cc2dfft(complex *data, int nx, int ldx, 
*                      int ny, int sign)
*
*   INPUT:  - *data: complex 2D input array
*           -    nx: number of x (fast) samples to be transformed
*           -   ldx: leading dimension in fast direction
*           -    ny: number of y (slow) samples to be transformed
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

void cc2dfft(complex *data, int nx, int ny, int ldx, int sign)
{
#if defined(CRAY_MPP)
	static int isys=0;
	static int nx_prev=0, ny_prev=0;
	static float *work, *table;
	static float *work2, *table2;
	static complex *tmp;
	int nwork, ntable, i, j;
	float scale=1.0;
#elif defined(CRAY_MPP_64)
	static int nx_prev=0, ny_prev=0;
	static int isys[3];
	int ntable, nwork, zero=0, i, j;
	static double *work, *table, scale=1.0, *ddata;
#elif defined(CRAY_PVP)
	static int nx_prev=0, ny_prev=0;
	int  isys=0, ntable, nwork, zero=0, i, j;
	int ld1;
	complex *tmp;
	static float *work, *table;
	float scale=1.0;
#elif defined(HAVE_LIBSCS)
	int pe;
	static int nx_prev[MAX_NUMTHREADS], ny_prev[MAX_NUMTHREADS];
	int   isys=0, ntable, nwork, zero=0;
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS];
	float scale=1.0;
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	int pe;
	static int nx_prev[MAX_NUMTHREADS], ny_prev[MAX_NUMTHREADS];
	static complex *coeff[MAX_NUMTHREADS];
#else
	int i, j;
	complex *tmp;
#endif

#if defined(CRAY_MPP)
	if (ny != 1) {
		if (nx != nx_prev) {
			nwork = 4*nx;
			ntable = 2*nx;
			if (work) free (work);
			work  = (float *)malloc(nwork*sizeof(float));
			if (work == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			if (table) free (table);
			table = (float *)malloc(ntable*sizeof(float));
			if (table == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			GGFFT(&isys,&nx,&scale,data,data,table,work,&isys);
			nx_prev = nx;
		}
		for (j=0; j<ny; j++) {
	 		GGFFT(&sign,&nx,&scale,&data[j*ldx].r,&data[j*ldx].r,table,work,&isys);
		}

		if (ny != ny_prev) {
			nwork = 4*ny;
			ntable = 2*ny;
			if (work2) free (work2);
			work2  = (float *)malloc(nwork*sizeof(float));
			if (work2 == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			if (table2) free (table2);
			table2 = (float *)malloc(ntable*sizeof(float));
			if (table2 == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			if (tmp) free (tmp);
			tmp = (complex *)malloc(ny*sizeof(complex));
			if (tmp == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			GGFFT(&isys,&ny,&scale,tmp,tmp,table2,work2,&isys);
			ny_prev = ny;
		}

		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) tmp[j] = data[j*ldx+i];
	 		GGFFT(&sign,&ny,&scale,tmp,tmp,table2,work2,&isys);
			for (j=0; j<ny; j++) data[j*ldx+i] = tmp[j];
		}
	}
	else {
		cc1fft(data, nx, sign);
	}
#elif defined(CRAY_MPP_64)
	if (ny != 1) {
		if (nx != nx_prev || ny != ny_prev) {
			if (factorized(nx) && factorized(ny)) {
				isys[0] = 2;
				isys[1] = 0;
				isys[2] = 0;
			 	ntable  = 2*(nx+ny);
			}
			else {
				isys[0] = 2;
				isys[1] = 1-factorized(nx);
				isys[2] = 1-factorized(ny);
			 	ntable  = 12*(nx+ny);
			}
			nwork = 2*nx*ny;
			if (work) free (work);
			work  = (double *)malloc(nwork*sizeof(double));
			if (work == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			if (table) free (table);
			table = (double *)malloc(ntable*sizeof(double));
			if (table == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			if (ddata) free(ddata);
			ddata = (double *)malloc(2*nx*ny*sizeof(double));
			if (ddata == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			CCFFT2D(&zero,&nx,&ny,&scale,ddata,&nx,ddata,&nx,table,work,&isys);
			nx_prev = nx;
			ny_prev = ny;
		}
		for (j=0; j<ny; j++) {
			for (i=0; i<nx; i++) {
				ddata[j*2*nx+2*i]   = (double) data[j*ldx+i].r;
				ddata[j*2*nx+2*i+1] = (double) data[j*ldx+i].i;
			}
		}
		CCFFT2D(&sign,&nx,&ny,&scale,ddata,&nx,ddata,&nx,table,work,&isys);
		for (j=0; j<ny; j++) {
			for (i=0; i<nx; i++) {
				data[j*ldx+i].r = (float) ddata[j*2*nx+2*i];
				data[j*ldx+i].i = (float) ddata[j*2*nx+2*i+1];
			}
		}
	}
	else {
		cc1fft(data, nx, sign);
	}
#elif defined(CRAY_PVP)
	if (ny != 1) {
		if (nx != nx_prev || ny != ny_prev) {
			nwork = 512*MAX(nx,ny);
			ntable = 100+2*(nx+ny);
			if (work) free (work);
			work  = (float *)malloc(nwork*sizeof(float));
			if (work == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			if (table) free (table);
			table = (float *)malloc(ntable*sizeof(float));
			if (table == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			CCFFT2D(&zero,&nx,&ny,&scale,data,&ldx,data,&ldx,table,work,&isys);
			nx_prev = nx;
			ny_prev = ny;
		}
		if ( !(ldx & 01) ) { /* leading dimension is not odd */
			ld1 = nx+1;
			tmp = (complex *) malloc(ld1*ny*sizeof(complex));
			if (tmp == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			for (j=0; j<ny; j++) {
				for (i=0; i<nx; i++) tmp[j*ld1+i] = data[j*ldx+i];
			}
			CCFFT2D(&sign,&nx,&ny,&scale,tmp,&ld1,tmp,&ld1,table,work,&isys);
			for (j=0; j<ny; j++) {
				for (i=0; i<nx; i++) data[j*ldx+i] = tmp[j*ld1+i];
			}
			free(tmp);
		}
		else {
			CCFFT2D(&sign,&nx,&ny,&scale,data,&ldx,data,&ldx,table,work,&isys);
		}
	}
	else {
		cc1fft(data, nx, sign);
	}
#elif defined(HAVE_LIBSCS)
	pe = mp_my_threadnum();
	if (ny != 1) {
		if (nx != nx_prev[pe] || ny != ny_prev[pe]) {
			nwork = 2*nx*ny;
			ntable = 60+2*(nx+ny);
			if (work[pe]) free (work[pe]);
			if (table[pe]) free (table[pe]);
			work[pe]  = (float *)malloc(nwork*sizeof(float));
			if (work[pe] == NULL) 
				fprintf(stderr,"cc2dfft: memory allocation error\n");
			table[pe] = (float *)malloc(ntable*sizeof(float));
			if (table[pe] == NULL) 
				fprintf(stderr,"cc2dfft: memory allocation error\n");
			ccfft2d_(&zero,&nx,&ny,&scale,data,&ldx,data,&ldx,
				table[pe],work[pe],&isys);
			nx_prev[pe] = nx;
			ny_prev[pe] = ny;
		}
		ccfft2d_(sign,&nx,&ny,&scale,data,&ldx,data,&ldx,
			table[pe],work[pe],&isys);
	}
	else {
		cc1fft(data, nx, sign);
	}
#elif defined(HAVE_LIBCOMPLIB_SGIMATH)
	pe = mp_my_threadnum();
	if (ny != 1) {
		if (nx != nx_prev[pe] || ny != ny_prev[pe]) {
			if (coeff[pe]) free(coeff[pe]);
			coeff[pe] = (complex *)cfft2di(nx, ny, NULL);
			nx_prev[pe] = nx;
			ny_prev[pe] = ny;
		}
		cfft2d(sign, nx, ny, data, ldx, coeff[pe]);
	}
	else {
		cc1fft(data, nx, sign);
	}
#else 
	if (ny != 1) {
		ccmfft(data, nx, ny, ldx, sign);
		tmp = (complex *)malloc(nx*ny*sizeof(complex));
		if (tmp == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) tmp[j+i*ny] = data[j*ldx+i];
		}
		ccmfft(tmp, ny, nx, ny, sign);
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) data[j*ldx+i] = tmp[j+i*ny];
		}
/*
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) tmp[j] = data[j*ldx+i];
			cc1fft(tmp, ny, sign);
			for (j=0; j<ny; j++) data[j*ldx+i] = tmp[j];
		}
*/
		free (tmp);
	}
	else {
		cc1fft(data, nx, sign);
	}
#endif

	return;
}

void free_cc2dfft()
{
#if defined(T3D) || defined(T3E)
	if (work2) free (work2);
	if (table2) free (table2);
	if (tmp) free (tmp);
	nx_prev = 0;
	ny_prev = 0;
#endif

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncc2dfft FNAME(CC2DFFTF)
#else
#define ncc2dfft FNAME(cc2dfftf)
#endif

void ncc2dfft(complex *data, int *nx, int *ny, int *ldx, int *sign)
{
	cc2dfft(data, *nx, *ny, *ldx, *sign);

	return;
}

