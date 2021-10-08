#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;

#define ISODD(n) ((n) & 01)

double powj(double b, double i);

void cSolveAxb(complex *A, int nrow, int ncol, complex *b, complex *x);

void coeffCosKr(float *coeff, int McC, int nkx, int nky, float dx, float dy, float dz, float k, float alpha, complex *x, int order, float wfact)
{
	int 	i, j, l, nk, iky, ikx, nrow, ncol;
	float	kx, ky, k2, kz2, ky2, kx2, kr, *weight, dkx, dky, kmax, *oper;
	complex *A, *b;
	
	nk	   = nkx*nky;
	A      = (complex *)malloc(nk*order*sizeof(complex));
	b      = (complex *)malloc(nk*sizeof(complex));
	weight = (float *)malloc(nk*sizeof(float));
	oper   = (float *)malloc(nk*sizeof(float));

	k2 	  = k*k;
	dkx  = M_PI/((nkx-1)*dx);
	dky  = M_PI/((nky-1)*dy);
	kmax = k*sin(alpha*M_PI/180.0);
	if (kmax > 0.85*M_PI/dx) kmax = 0.85*M_PI/dx;

	if (McC == 1) {
		for (iky = 0; iky < nky; iky++) {
			ky   = iky*dkx;
			for (ikx = 0; ikx < nkx; ikx++) {
				kx = ikx*dkx;
				oper[ikx*nkx+iky] = coeff[0] + coeff[1]*cos(kx*dx) +
					coeff[2]*cos(ky*dy) + coeff[3]*cos(kx*dx)*cos(ky*dy);
			}
		}
	}
	else if (McC == 2) {
		for (iky = 0; iky < nky; iky++) {
			ky   = iky*dkx;
			for (ikx = 0; ikx < nkx; ikx++) {
				kx = ikx*dkx;
				oper[ikx*nkx+iky] = coeff[0] + coeff[1]*cos(kx*dx) +
					coeff[3]*cos(ky*dy) + coeff[4]*cos(kx*dx)*cos(ky*dy) + 
					coeff[2]*cos(2.0*kx*dx) + coeff[6]*cos(2.0*ky*dy) +
					coeff[5]*cos(2.0*kx*dx)*cos(ky*dy) +
					coeff[7]*cos(kx*dx)*cos(2.0*ky*dy) +
					coeff[8]*cos(2.0*kx*dx)*cos(2.0*ky*dy);
			}
		}
	}


/*  Define matrix A, vector b and the circular weighting function */

	for (iky = 0; iky < nky; iky++) {
		ky   = iky*dky;
		ky2  = ky*ky;
		for (ikx = 0; ikx < nkx; ikx++) {
			kx = ikx*dkx;
			kx2 = kx*kx;
			kr = sqrt(kx2+ky2);
			if (kr <= kmax) {
				A[ikx+iky*nkx].r = 1.0;
				A[ikx+iky*nkx].i = 0.0;
				weight[ikx+iky*nkx] = 1.0;
			}
			else{
				A[ikx+iky*nkx].r = wfact;
				A[ikx+iky*nkx].i = 0.0;
				weight[ikx+iky*nkx] = wfact;
			}
		}
	}

	for(i = 1; i < order; i++) {
		for(l = 0; l < nky; l++) {
			for(j = 0; j < nkx; j++) {
				A[i*nk+j+l*nkx].r = powj(oper[l*nkx+j], i);
				A[i*nk+j+l*nkx].r *= weight[j+l*nkx];
				A[i*nk+j+l*nkx].i = 0.0;
			}
		}
	}

/*  Defining known function b */

	for(l = 0; l < nky; l++) {
		ky = l*dky;
		for(j = 0; j < nkx; j++) {
			kx = j*dkx;
			kz2 = k2 - (ky*ky+kx*kx);
			if (kz2 > 0) {
				b[j+l*nkx].r = cos(sqrt(kz2)*dz);
				b[j+l*nkx].i = sin(sqrt(kz2)*dz);
			}
			else {
				b[j+l*nkx].r = exp(-sqrt(-kz2)*dz);
				b[j+l*nkx].i = 0.0;
			} 
			b[j+l*nkx].r *= weight[j+l*nkx];
			b[j+l*nkx].i *= weight[j+l*nkx];
		}
	}

	ncol = order;
	nrow = nk;

	cSolveAxb(A, nrow, ncol, b, x);

	free(A);
	free(b);
	free(weight);
	free(oper);

	return;
}

double powj(double b, double i)
{
	if (b >= 0 )
		return pow(b, i);
	else {
		if (ISODD((int)i))
			return -1.0*pow(fabs(b), i);
		else
			return pow(fabs(b), i);
	}
}

