#include <genfft.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void fft(int n, float *real, float *imag);
void ifft(int n, float *real, float *imag);
void realifft(int n, float *real);
void realfft(int n, float *real);
int npfa (int nmin);
void pfacc (int isign, int n, complex z[]);
void ccdft(complex *cdata, int n, int sign);
void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);


/**
*
*   DESCRIPTION: Local FFT implementation
*
*   USAGE:
*	     void rc1_fft(float *data, complex *cdata, int n, int sign)
*        void rcm_fft(float *data, complex *cdata, int n1, int n2, int sign)
*        void cr1_fft(complex *cdata, float *data, int n, int sign)
*        void crm_fft(complex *cdata, float *data, int n1, int n2, int sign)
*        void cc1_fft(complex *cdata, int n, int sign)
*        void ccm_fft(complex *cdata, int n1, int n2, int sign)
*
*   INPUT:  see the documentation in the files without '_'
*
*   OUTPUT: see the documentation in the files without '_'
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
*    1.0       Jan Thorbecke    Feb  '94    Initial version with NumRec
*    1.1       Jan Thorbecke    June '94    faster Mayer local FFT 
*    2.0       Jan Thorbecke    July '97    make it more general
*
*           Silicon Graphics / Cray Research
*           Veldzigt 2a
*           3454 PW De Meern
*           The Netherlands
*
***********************************************************************/

void rc1_fft(float *data, complex *cdata, int n, int sign)
{
	int    j;
	float *datft;

	datft = (float *)malloc(n*sizeof(float));
	if (datft == NULL) fprintf(stderr,"rc1_fft: memory allocation error\n");

	for (j = 0; j < n; j++) datft[j] = (float)data[j];
	realfft(n, datft);
	cdata[0].i = 0.0;
	for (j = 0; j < n/2; j++) {
		cdata[j].r = (float)datft[j]; 
		cdata[j+1].i = sign*(float)datft[n-j-1]; 
	}
	cdata[n/2].r = datft[n/2];
	cdata[n/2].i = 0.0; 

	free(datft);

	return;
}

void rcm_fft(float *data, complex *cdata, int n1, int n2, int ldr, int ldc, int sign)
{
	int    j, i;
	float *datft;

	datft = (float *)malloc(n1*sizeof(float));
	if (datft == NULL) fprintf(stderr,"rc2_fft: memory allocation error\n");

	for (i = 0; i < n2; i++) {
		for (j = 0; j < n1; j++) datft[j] = (float)data[i*ldr+j];
		realfft(n1, datft);
		cdata[i*ldc].i = 0.0;
		for (j = 0; j < n1/2; j++) {
			cdata[i*ldc+j].r = (float)datft[j]; 
			cdata[i*ldc+j+1].i = sign*(float)datft[n1-j-1]; 
		}
		cdata[i*ldc+n1/2].r = (float)datft[n1/2]; 
		cdata[i*ldc+n1/2].i = 0.0; 
	}

	free(datft);

	return;
}

void cr1_fft(complex *cdata, float *data, int n, int sign)
{
	int    j;
	float *datft;

	datft = (float *)malloc(n*sizeof(float));
	if (datft == NULL) fprintf(stderr,"cr1_fft: memory allocation error\n");

	for (j = 0; j < n/2; j++) {
		datft[j] = (float)cdata[j].r;
		datft[n-1-j] = (float)cdata[j+1].i;
	}
	datft[n/2] = (float)cdata[n/2].r;

	realifft(n, datft);

	if (sign == -1) {
		for (j = 0; j < n; j++) 
			data[j] = (float)datft[j];
	}
	else if (sign == 1) {
		for (j = 1; j < n; j++) 
			data[j] = (float)datft[n-j];
		data[0] = (float)datft[0];
	}

	free(datft);

	return;
}

void crm_fft(complex *cdata, float *data, int n1, int n2, int ldc, int ldr, int sign)
{
	int    j, i;
	float *datft;

	datft = (float *)malloc(n1*sizeof(float));
	if (datft == NULL) fprintf(stderr,"cr2_fft: memory allocation error\n");

	for (i = 0; i < n2; i++) {
		for (j = 0; j < n1/2; j++) {
			datft[j] = (float)cdata[i*ldc+j].r;
			datft[n1-1-j] = (float)cdata[i*ldc+j+1].i;
		}
		datft[n1/2] = (float)cdata[i*ldc+n1/2].r;

		realifft(n1, datft);

		if (sign == -1) {
			for (j = 0; j < n1; j++) 
				data[i*ldr+j] = (float)datft[j];
		}
		else if (sign == 1) {
			for (j = 1; j < n1; j++) 
				data[i*ldr+j] = (float)datft[n1-j];
			data[i*ldr] = (float)datft[0];
		}
	}

	free(datft);

	return;
}


void cc1_fft(complex *cdata, int n, int sign)
{
	int    j;
	float  *real, *imag;

	if (NINT(pow(2.0, (double)NINT(log((double)n)/log(2.0)))) != n) {
		if (npfa(n) == n) pfacc(sign, n, cdata);
		else ccdft(cdata,n,sign);
	}
	else {
		real = (float *)malloc(n*sizeof(float));
		if (real == NULL) fprintf(stderr,"cc1_fft: memory allocation error\n");
		imag = (float *)malloc(n*sizeof(float));
		if (imag == NULL) fprintf(stderr,"cc1_fft: memory allocation error\n");
	
		for (j = 0; j < n; j++) {
			real[j] = (float)cdata[j].r;
			imag[j] = (float)cdata[j].i;
		}

		if (sign < 0) fft(n, real, imag);
		else ifft(n, real, imag);

		for (j = 0; j < n; j++) {
			cdata[j].r = (float)real[j];
			cdata[j].i = (float)imag[j];
		}

		free(real);
		free(imag);
	}

	return;
}

void ccm_fft(complex *cdata, int n1, int n2, int ld1, int sign)
{
	int    i, j;
	float  *real, *imag;

	if (NINT(pow(2.0, (double)NINT(log((double)n1)/log(2.0)))) != n1) {
		if (npfa(n1) == n1) pfamcc(sign, n1, n2, 1, ld1, cdata);
		else {
			for (i = 0; i < n2; i++) {
				ccdft(&cdata[i*ld1],n1,sign);
			}
		}
	}
	else {
		real = (float *)malloc(n1*sizeof(float));
		if (real == NULL) fprintf(stderr,"ccm_fft: memory allocation error\n");
		imag = (float *)malloc(n1*sizeof(float));
		if (imag == NULL) fprintf(stderr,"ccm_fft: memory allocation error\n");
	
		for (i = 0; i < n2; i++) {
			for (j = 0; j < n1; j++) {
				real[j] = cdata[i*ld1+j].r;
				imag[j] = cdata[i*ld1+j].i;
			}
	
			if (sign < 0) fft(n1, real, imag);
			else ifft(n1, real, imag);
	
			for (j = 0; j < n1; j++) {
				cdata[i*ld1+j].r = real[j];
				cdata[i*ld1+j].i = imag[j];
			}

		}

		free(real);
		free(imag);
	}

	return;
}

