#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "fdelmod.h"
#include "segyhdr.h"
#include "par.h"

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#endif/* complex */

int optncr(int n);
void rc1fft(float *rdata, complex *cdata, int n, int sign);
void cr1fft(complex *cdata, float *rdata, int n, int sign);

int writesufile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2);
float gaussGen();
void spline3(float x1, float x2, float z1, float z2, float dzdx1,
         float dzdx2, float *a, float *b, float *c, float *d);

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int defineSource(wavPar wav, srcPar src, float *src_nwav, int verbose)
{
    FILE    *fp;
    size_t  nread;
    int optn, nfreq, i, j, k, iwmax, tracesToDo;
	int iw, n1, itbeg, itmax, nsmth;
    float scl, d1, df, deltom, om;
	float amp1, amp2, amp3;
	float *trace, maxampl;
	float   x1, x2, z1, z2, dzdx1, dzdx2, a, b, c, d, t;
    complex *ctrace, tmp;
    segyhdr hdr;
    
	n1 = wav.nt;
	if (wav.random) {
		srand48(10);
	}
	else {

/* read first header and last byte to get file size */

    	fp = fopen( wav.file_src, "r" );
    	assert( fp != NULL);
    	nread = fread( &hdr, 1, TRCBYTES, fp );
    	assert(nread == TRCBYTES);
	
/* read all traces */

		tracesToDo = wav.nx;
		i = 0;
		while (tracesToDo) {
			memset(&src_nwav[i*n1],0,n1*sizeof(float));
        	nread = fread(&src_nwav[i*n1], sizeof(float), hdr.ns, fp);
        	assert (nread == hdr.ns);
	
        	nread = fread( &hdr, 1, TRCBYTES, fp );
        	if (nread==0) break;
			tracesToDo--;
			i++;
		}
    	fclose(fp);
	}

/*=========== Scale source with 1/jw and set maximum amplitude to 1 =========*/

	optn = optncr(2*n1);
	nfreq = optn/2 + 1;
	ctrace = (complex *)calloc(nfreq,sizeof(complex));
	trace = (float *)calloc(optn,sizeof(float));

	df     = 1.0/(optn*wav.dt);
    deltom = 2.*M_PI*df;
	scl    = 1.0/optn;
	maxampl=0.0;
    
	if (wav.random) iwmax = MIN(NINT(wav.fmax/df),nfreq);
	else iwmax = nfreq;
//	fprintf(stderr,"fmax=%f iwmax = %d nfreq=%d dt=%f df=%f src_n=%d wav.nx=%d\n", wav.fmax, iwmax, nfreq, wav.dt, deltom, src.n, wav.nx);

	for (i=0; i<wav.nx; i++) {
		if (wav.random) {
			for (iw=1;iw<iwmax;iw++) {
				ctrace[iw].r = (float)(drand48()-0.5);
				ctrace[iw].i = (float)(drand48()-0.5);
			}
		}
		else {
			memset(&ctrace[0].r,0,nfreq*sizeof(complex));
			memset(&trace[0],0,optn*sizeof(float));
			memcpy(&trace[0],&src_nwav[i*n1],n1*sizeof(float));
			rc1fft(trace,ctrace,optn,-1);
			for (iw=1;iw<iwmax;iw++) {
				om = 1.0/(deltom*iw);
				tmp.r = om*ctrace[iw].i;
				tmp.i = -om*ctrace[iw].r;
				ctrace[iw].r = tmp.r;
				ctrace[iw].i = tmp.i;
			}

			/* zero frequency iw=0 */
			amp1 = sqrt(ctrace[1].r*ctrace[1].r+ctrace[1].i*ctrace[1].i);
			if (amp1 == 0.0) {
				ctrace[0].r = ctrace[0].i = 0.0;
			}
			else {
				amp2 = sqrt(ctrace[2].r*ctrace[2].r+ctrace[2].i*ctrace[2].i);
				amp3 = sqrt(ctrace[3].r*ctrace[3].r+ctrace[3].i*ctrace[3].i);
				ctrace[0].r = amp1+(2.0*(amp1-amp2)-(amp2-amp3));
				ctrace[0].i = 0.0;
				if (ctrace[1].r < 0.0) {
					ctrace[0].r *= -1.0;
				}
//				fprintf(stderr,"amp1=%f amp2=%f amp3=%f\n", amp1, amp2, amp3);
//				fprintf(stderr,"ctrace[0]=%f =%f \n", ctrace[0].r, ctrace[0].i);
//				fprintf(stderr,"ctrace[1]=%e =%e \n", ctrace[1].r, ctrace[1].i);
//				fprintf(stderr,"ctrace[5]=%e =%e \n", ctrace[5].r, ctrace[5].i);
//				fprintf(stderr,"ctrace[50]=%e =%e \n", ctrace[50].r, ctrace[50].i);
			}
		}
		for (iw=iwmax;iw<nfreq;iw++) {
			ctrace[iw].r = 0.0;
			ctrace[iw].i = 0.0;
		}
		cr1fft(ctrace,trace,optn,1);

		if (wav.random) {
			/* find first zero crossing in wavelet */
			amp1 = trace[0];
			j = 1;
			if (amp1 < 0.0) {
				while (trace[j] < 0.0) j++;
			}
			else {
				while (trace[j] > 0.0) j++;
			}
			itbeg = j;

			/* find last zero crossing in wavelet */
			itmax = itbeg+MIN(NINT((src.tend[i]-src.tbeg[i])/wav.dt),n1);
			amp1 = trace[itmax-1];
			j = itmax;
			if (amp1 < 0.0) {
				while (trace[j] < 0.0) j++;
			}
			else {
				while (trace[j] > 0.0) j++;
			}
			itmax = j;

			/* make smooth transitions to zero aamplitude */
			nsmth=10;
			x1 = 0.0;
			z1 = 0.0;
			dzdx1 = 0.0;
			x2 = nsmth;
			z2 = trace[itbeg+nsmth];
//			dzdx2 = (trace[itbeg+(nsmth+1)]-trace[itbeg+(nsmth-1)])/(2.0);
			dzdx2 = (trace[itbeg+nsmth-2]-8.0*trace[itbeg+nsmth-1]+
					 8.0*trace[itbeg+nsmth+1]-trace[itbeg+nsmth+2])/(12.0);
			spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);
			for (j=0; j<nsmth; j++) {
				t = j;
				trace[itbeg+j] = a*t*t*t+b*t*t+c*t+d;
			}

			x1 = 0.0;
			z1 = trace[itmax-nsmth];
//			dzdx1 = (trace[itmax-(nsmth-1)]-trace[itmax-(nsmth+1)])/(2.0);
			dzdx1 = (trace[itmax-nsmth-2]-8.0*trace[itmax-nsmth-1]+
					 8.0*trace[itmax-nsmth+1]-trace[itmax-nsmth+2])/(12.0);
			x2 = nsmth;
			z2 = 0.0;
			dzdx2 = 0.0;

//			fprintf(stderr,"x1=%f z1=%f d=%f x2=%f, z2=%f d=%f\n",x1,z1,dzdx1,x2,z2,dzdx2);
			spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);
			for (j=0; j<nsmth; j++) {
				t = j;
				trace[itmax-nsmth+j] = a*t*t*t+b*t*t+c*t+d;
//			fprintf(stderr,"a=%f b=%f c=%f d=%f trace%d=%f \n",a,b,c,d,j,trace[itmax-nsmth+j]);
			}


			for (j=itbeg; j<itmax; j++) src_nwav[i*n1+j-itbeg] = trace[j];
		}
		else {
			for (j=0; j<n1; j++) src_nwav[i*n1+j] = scl*(trace[j]-trace[optn-1]);
		}

//		for (j=0; j<n1; j++) {
//			maxampl = MAX(maxampl, src_nwav[i*n1+j]);
//		}
    }
    free(ctrace);
    free(trace);

	if (src.amplitude > 0.0) {
		trace = (float *)calloc(100,sizeof(float));
		for (i=0; i<wav.nx; i++) {
			if (src.distribution) scl = gaussGen()*src.amplitude;
			else scl = (float)(drand48()-0.5)*src.amplitude;
			k = (int)MAX(MIN(100*(scl+2.5*src.amplitude)/(5*src.amplitude),99),0);
			trace[k] += 1.0;
			for (j=0; j<n1; j++) {
				src_nwav[i*n1+j] *= scl;
			}
    	}
		d1 = 5*src.amplitude/99;
		if (verbose>2) writesufile("src_ampl.su", trace, 100, 1, -2.5*src.amplitude, 0.0, d1, 1);
		free(trace);
	}

	if (verbose>3) writesufile("src_nwav.su", src_nwav, n1, wav.nx, 0.0, 0.0, wav.dt, 1);

/* set maximum amplitude in source file to 1.0 */
/*
	assert(maxampl != 0.0);
	scl = wav.dt/(maxampl);
	scl = 1.0/(maxampl);
	for (i=0; i<wav.nx; i++) {
		for (j=0; j<n1; j++) {
			src_nwav[i*n1+j] *= scl;
		}
    }
*/

    return 0;
}

