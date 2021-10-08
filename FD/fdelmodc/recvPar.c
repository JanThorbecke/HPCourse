#include <stdio.h>
#include <assert.h>
#include "par.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void name_ext(char *filename, char *extension);

int recvPar(int *rec_x, int *rec_z, int *nrcv, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz)
{
	float *xrcv1, *xrcv2, *zrcv1, *zrcv2;
	int ix0, ix1, iz0, iz1, ix, iz, isign, verbose;
	float dxrcv, dzrcv, r;
	int Nx1, Nx2, Nz1, Nz2, iarray, nskip, nrec;

	Nx1 = countparval("xrcv1");
	Nx2 = countparval("xrcv2");
	Nz1 = countparval("zrcv1");
	Nz2 = countparval("zrcv2");
	if (Nx1 != 1 || Nz1 != 1) {
		assert(Nx1==Nx2);
		assert(Nz1==Nz2);
		assert(Nx1==Nz1);
	}
	if (Nx1==0) Nx1=1;

	xrcv1 = (float *)malloc(Nx1*sizeof(float));
	xrcv2 = (float *)malloc(Nx1*sizeof(float));
	zrcv1 = (float *)malloc(Nx1*sizeof(float));
	zrcv2 = (float *)malloc(Nx1*sizeof(float));

	if(!getparint("verbose", &verbose)) verbose = 0;
	if(!getparfloat("xrcv1", xrcv1)) xrcv1[0]=sub_x0;
	if(!getparfloat("xrcv2", xrcv2)) xrcv2[0]=(nx-1)*dx+sub_x0;
	if(!getparfloat("dxrcv", &dxrcv)) dxrcv=dx;
	if(!getparfloat("zrcv1", zrcv1)) zrcv1[0]=sub_z0;
	if(!getparfloat("zrcv2", zrcv2)) zrcv2[0]=zrcv1[0];
	if(!getparfloat("dzrcv", &dzrcv)) dzrcv=0;

	nrec=0;
	for (iarray=0; iarray<Nx1; iarray++) {

		ix0=MAX(0,NINT((xrcv1[iarray]-sub_x0)/dx));
		ix0=MIN(ix0,nx-1);
		ix1=MAX(0,NINT((xrcv2[iarray]-sub_x0)/dx));
		ix1=MIN(ix1,nx-1);
	
		iz0=MAX(0,NINT((zrcv1[iarray]-sub_z0)/dz));
		iz0=MIN(iz0,nz-1);
		iz1=MAX(0,NINT((zrcv2[iarray]-sub_z0)/dz));
		iz1=MIN(iz1,nz-1);

		if (abs(ix1-ix0) <= abs(iz1-iz0)) {
			nskip = MAX(1,NINT(dzrcv/dz));
		}
		else if (abs(ix1-ix0) > abs(iz1-iz0)) {
			nskip = MAX(1,NINT(dxrcv/dx));
		}

// calculate receiver array and store into rec_x,z

		if ((ix1==ix0) && (iz1==iz0)) {
			rec_x[nrec]=ix0;
			rec_z[nrec]=iz0;
			nrec++;
		}
		else if (abs(ix1-ix0) <= abs(iz1-iz0)) {	
// compute derivative
			if (abs(iz1-iz0) != 0) {
				r=(1.*ix1-ix0)/(iz1-iz0);
				isign=abs(iz1-iz0)/(iz1-iz0);
			}
				else {
				r=0.0;
				isign=1;
			}
// calculate coordinates
			for (iz=iz0; iz<=iz1; iz+=nskip*isign) {
				rec_x[nrec]=NINT( (iz-iz0)*r+ix0 );
				rec_z[nrec]=iz;
				nrec++;
			}
		}
		else if (abs(ix1-ix0) > abs(iz1-iz0)) {	
// compute derivative
			if (abs(ix1-ix0)!=0) {
				r=(1.*iz1-iz0)/(ix1-ix0);
				isign=abs(ix1-ix0)/(ix1-ix0);
			}
			else {
				r = 0.0;
				isign = 1;
			}
// calculate coordinates
			for (ix=ix0; ix<=ix1; ix+=nskip*isign) {
				rec_x[nrec]=ix;
				rec_z[nrec]=NINT( (ix-ix0)*r+iz0 );
				nrec++;
			}
		}
	}

	*nrcv = nrec;
	free(xrcv1);
	free(xrcv2);
	free(zrcv1);
	free(zrcv2);

	return 0;
}

