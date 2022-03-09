#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"par.h"

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int main(int argc, char *argv[]) 
{
	float *xr, *zr, dx, dz, sub_x0, sub_z0, dxrcv, dzrcv;
	int *x, *z, ir, nrec, nrecv;

    initargs(argc, argv);

    nrecv=501;

	getparfloat("dx", &dx);
	if(!getparfloat("dz", &dz)) dz=dx;
	getparfloat("dxrcv", &dxrcv);
	if(!getparfloat("dzrcv", &dzrcv)) dzrcv=0.0;
	getparfloat("x0", &sub_x0);
	if(!getparfloat("z0", &sub_z0)) sub_z0=0.0;

	fprintf(stderr,"dx=%f dz=%f dxrcv=%f dzrcv=%f\n", dx, dz, dxrcv, dzrcv);

	//sub_x0=sub_z0=dx*0.5;

	xr = (float *)malloc(nrecv*sizeof(float));
	zr = (float *)malloc(nrecv*sizeof(float));
	x  = (int *)malloc(nrecv*sizeof(int));
	z  = (int *)malloc(nrecv*sizeof(int));

	nrec=0;
	for (ir=0; ir<nrecv; ir++) {
		xr[nrec]=sub_x0+ir*dxrcv;
		zr[nrec]=sub_z0+ir*dzrcv;
	
		x[nrec]=NINT((xr[nrec])/dx);
		z[nrec]=NINT((zr[nrec])/dz);
		nrec++;
	}

	for (ir=1; ir<nrecv; ir++) {
		if (x[ir]-x[ir-1] != 4) fprintf(stderr,"xr=%f zr=%f at grid coordinates ix=%d iz=%d\n",xr[ir], zr[ir], x[ir], z[ir]);
	}


	return 0;
}

