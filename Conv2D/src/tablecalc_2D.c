#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct _complexStruct { 
    float r,i;
} complex;

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

static float dkx, kmin;
static complex *table;

void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, float dz, float alpha, float fmin, float fmax, float cmin, float cmax, float df, float weight, int fine, int method, int verbose)
{
	int 	ikx, nkx, nky, nkx2, nky2, hoplx, hoply, ntable, ix, iy;
	int     size;
	float	k, kmax, wfact, fampl, w_start, limit;
	float   ampl, phase;
	complex *hopkx, *hopx;

	kmin   = 2.0*M_PI*(MAX(fmin-df,0))/cmax;
	kmax   = 2.0*M_PI*(fmax+df)/cmin;
	dkx    = 2.0*M_PI*df/(float)(cmax*fine);
	ntable = (int)((kmax - kmin)/dkx)+1;
	nkx	   = pow(2.0, ceil(log(nx)/log(2.0)));
	nky	   = pow(2.0, ceil(log(ny)/log(2.0)));
	nkx = nky = 512;
	nkx2   = nkx/2+1;
	nky2   = nky/2+1;
	hoplx  = (oplx+1)/2;
	hoply  = (oply+1)/2;
	size   = hoplx*hoply;
    w_start= 1e-6;
	limit  = 1.002;

	table = (complex *)malloc(hoplx*hoply*ntable*sizeof(complex));
	assert (table != NULL);

	if (verbose) {
		fprintf(stderr,"Number of operators to calculate = %d\n", ntable);
		fprintf(stderr,"Size of operator table = %d bytes \n", (int)sizeof(complex)*ntable*hoplx*hoply);
	}

/* WLSQ operator */


	k = kmin;
	for (ikx = 0; ikx < ntable; ikx++) {
		for (iy = 0; iy < hoply; iy++) {
			for (ix = 0; ix < hoplx; ix++) {
				ampl  = MIN(((float) rand()/RAND_MAX),0.99);
				ampl  = 1.001;
				phase = 2.0*M_PI*((float) rand()/RAND_MAX);
				phase = M_PI*0.5;
				table[ikx*size+iy*hoplx+ix].r = ampl*cos(phase);
				table[ikx*size+iy*hoplx+ix].i = ampl*sin(phase);
			}
		}
/*
			table[ikx*size].r = 0.25*hopx[0].r;
			table[ikx*size].i = 0.25*hopx[0].i;
			for (ix = 1; ix < hoplx; ix++) {
				table[ikx*size+ix].r = 0.5*hopx[ix].r;
				table[ikx*size+ix].i = 0.5*hopx[ix].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				table[ikx*size+iy*hoplx].r = 0.5*hopx[iy*hoplx].r;
				table[ikx*size+iy*hoplx].i = 0.5*hopx[iy*hoplx].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				for (ix = 1; ix < hoplx; ix++) {
					table[ikx*size+iy*hoplx+ix] = hopx[iy*hoplx+ix];
				}
			}
*/
	}

	return;
}

void readtable2D(complex *oper, float k, int hoplx, int hoply, int mode)
{
	int p1, p2, i, size;
	float linscale;

	size = hoplx*hoply;
	p1 = (int)((k-kmin)/dkx);
	p2 = p1+1;
	linscale = (k-kmin)/dkx - (float)p1;

	for (i = 0; i < hoplx*hoply; i++) {
		oper[i].r = (table[p1*size+i].r + linscale*(table[p2*size+i].r - 
					table[p1*size+i].r));
		oper[i].i = (table[p1*size+i].i + linscale*(table[p2*size+i].i - 
					table[p1*size+i].i))*mode;
	}

	return;
}


