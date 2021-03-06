#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "segyhdr.h"
#include <math.h>
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

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int getWaveletInfo(char *file_src, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *fmax, int *nxm, int verbose)
{
    FILE    *fp;
    size_t  nread, data_sz;
	off_t bytes, ret, trace_sz, ntraces;
    int sx_shot, gx_start, one_shot;
    int optn, nfreq, i, iwmax;
    float *trace;
	float ampl, amplmax; 
    complex *ctrace;
    segyhdr hdr;
    
    if (file_src == NULL) return 0; /* Input pipe can not be handled */
    else fp = fopen( file_src, "r" );
    assert( fp != NULL);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
	if (ret<0) perror("fseeko");
    bytes = ftello( fp );

    *n1 = hdr.ns;
    if (hdr.trid == 1 || hdr.dt != 0) {
        *d1 = ((float) hdr.dt)*1.e-6;
        *f1 = ((float) hdr.delrt)/1000.;
		if (*d1 == 0.0) *d1 = hdr.d1;
    }
    else {
        *d1 = hdr.d1;
        *f1 = hdr.f1;
    }
    *f2 = hdr.f2;

    data_sz = sizeof(float)*(*n1);
    trace_sz = sizeof(float)*(*n1)+TRCBYTES;
    ntraces  = (int) (bytes/trace_sz);
	*n2 = ntraces;
//	fprintf(stderr,"data_sz %ld trace_sz %lld  bytes = %lld\n", data_sz, trace_sz, bytes);

    /* check to find out number of traces in shot gather */

	optn = optncr(hdr.ns);
	nfreq = optn/2 + 1;
	ctrace = (complex *)malloc(nfreq*sizeof(complex));
    one_shot = 1;
    sx_shot  = hdr.sx;
    gx_start = hdr.gx;
    trace = (float *)malloc(optn*sizeof(float));
    fseeko( fp, TRCBYTES, SEEK_SET );

    while (one_shot) {
		memset(trace,0,optn*sizeof(float));
        nread = fread( trace, sizeof(float), hdr.ns, fp );
        assert (nread == hdr.ns);
		rc1fft(trace,ctrace,optn,1);
		amplmax = 0.0;
		for (i=0;i<nfreq;i++) {
			ampl = sqrt(ctrace[i].r*ctrace[i].r+ctrace[i].i*ctrace[i].i);
			if (ampl > amplmax) {
				amplmax = ampl;
				iwmax = i;
			}
		}
		for (i=iwmax;i<nfreq;i++) {
			ampl = sqrt(ctrace[i].r*ctrace[i].r+ctrace[i].i*ctrace[i].i);
			if (400*ampl < amplmax) {
				*fmax = (i-1)*(1.0/(optn*(*d1)));
				break;
			}
		}

        nread = fread( &hdr, 1, TRCBYTES, fp );
        if (nread==0) break;
    }
	*nxm = (int)ntraces;

	if (verbose>2) {
		vmess("For file %s", file_src);
		vmess("nt=%d nx=%d", *n1, *n2);
		vmess("dt=%f dx=%f", *d1, *d2);
		vmess("fmax=%f", *fmax);
		vmess("tstart=%f", *f1);
	}

    fclose(fp);
    free(trace);
    free(ctrace);

    return 0;
}
