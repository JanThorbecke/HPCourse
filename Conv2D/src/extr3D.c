#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define ISODD(n) ((n) & 01)

double wallclock_time(void);

/***** Operator tables *****/
void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, float dz, float alpha, float fmin, float fmax, float cmin, float cmax, float df, float weight, int fine, int method, int verbose);

/***** Direct convolution *****/
void conv2D(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode);

void conv2D_q8(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode);

void conv2D_q4(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode);

void conv2D_q1(int nx, int ny, complex *data, float *velmod, int oplx, int oply, float om, int mode);

void extr3D(complex *wavelet, float df, float fmin, float fmax, int nfreq, int nx, float dx, int ny, float dy, int oplx, int oply, float dz, int nz, float alpha, float cp, float wfacto, float wfacts, float *data3d, int order, int McC, int opt, int d1op, int oper, int verbose)
{
	int    	iomin, iomax, iom, ix, iy, d, j;
	int     i, fine, method=1, mode=1;
	float  	deltom, om, omax, omin, *p;
	double  t0, t2, top;
	float   *velmod, cmin, cmax;
	complex	*data;

	deltom	= 2.*M_PI*df;
	omin	= 2.*M_PI*fmin;
	omax	= 2.*M_PI*fmax;
	iomin	= (int)MIN((omin/deltom), (nfreq-1));
	iomin	= MAX(iomin, 1);
	iomax	= MIN((int)(omax/deltom), (nfreq-1));

	fine  = 1;

	data  = (complex *)malloc(nx*ny*sizeof(complex));
	velmod  = (float *)malloc(nx*ny*nz*sizeof(float));

/* determine velocity model */

	for (d=0; d<nz; d++) {
		for (i=0; i<ny; i++) {
			for (j=0; j<nx; j++) {
				velmod[d*nx*ny+i*nx+j] = cp;
			}
		}
	}

/* find minimum and maximum value in velocity model */

	cmin = cmax = velmod[0];
	for (i=0; i<nx*ny*nz; i++) {
		if (velmod[i] > cmax) cmax = velmod[i];
		if (velmod[i] < cmin) cmin = velmod[i];
	}
	if (verbose) fprintf(stderr,"cmin = %.2f cmax = %.2f\n", cmin, cmax);
	
/* Calculate operator table */

	t0 = wallclock_time();
	method = opt;
	tablecalc_2D(oplx, oply, nx, dx, ny, dy, dz, alpha, fmin, fmax, 
		cmin, cmax, df, wfacto, fine, method, verbose);
	top = wallclock_time() - t0;

	if (verbose) {
		switch ( oper ) {
			case 1 :
				fprintf(stderr," using simple convolution: conv2D_q1\n");
				break;
			case 4 :
				fprintf(stderr," using symmetric convolution: conv2D_q4\n");
				break;
			case 8 :
				fprintf(stderr," using symmetric all convolution: conv2D_q8\n");
				break;
			case 0 :
				fprintf(stderr," using vector convolution: conv2D\n");
				break;
		}
	}

/* Extrapolation for all frequencies of interest */

	t0 = wallclock_time();
	for (iom = iomin; iom <= iomax; iom++) {

		p = &data[0].r;
		for (j = 0; j < 2*nx*ny; j++) *p++ = 0.0;

		if (ISODD(nx)) ix = (nx-1)/2;
		else ix = nx/2;
		if (ISODD(ny)) iy = (ny-1)/2;
		else iy =  ny/2;

		data[iy*nx+ix].r = wavelet[iom].r;
		data[iy*nx+ix].i = wavelet[iom].i;

		om = iom*deltom;

		if (verbose) fprintf(stderr,"Extr iom = %d to do %d freq \n",iom,iomax-iom+1);

		for (d = 0; d < nz; d++) {

			switch ( oper ) {
				case 1 :
					conv2D_q1(nx,ny,data,&velmod[d*nx*ny],oplx,oply,om,mode);
					break;
				case 4 :
					conv2D_q4(nx,ny,data,&velmod[d*nx*ny],oplx,oply,om,mode);
					break;
				case 8 :
					conv2D_q8(nx,ny,data,&velmod[d*nx*ny],oplx,oply,om,mode);
					break;
				case 0 :
					conv2D(nx,ny,data,&velmod[d*nx*ny],oplx,oply,om,mode);
					break;
			}

			/* Imaging condition */
			for (iy = 0; iy < ny; iy++) {
				for (ix = 0; ix < nx; ix++) {
					data3d[d*nx*ny+iy*nx+ix] += 2.0*data[iy*nx+ix].r;
				}
			}
		}
	}

	t2 = wallclock_time();
	fprintf(stderr, "- Operator Computation .......... = %f s.\n",top);
	fprintf(stderr, "- Extrapolation Computation ..... = %f s.\n",t2-t0);

	free(data);
	free(velmod);

	return;
}

