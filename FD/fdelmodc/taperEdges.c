#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmod.h"

int taperEdges(modPar mod, bndPar bnd, float *vx, float *vz, int verbose)
{
	int   ix, iz, ibnd, ib, ntaper;
	int   nx, nz, n1;

	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;
	ibnd = mod.iorder/2-1;

	/* top */
	if (bnd.tap[0] > 0) {
		ntaper = bnd.tap[0];
		ib = (ntaper+ibnd-1);
		for (ix=ibnd; ix<nx+ibnd; ix++) {
#pragma ivdep
			for (iz=ibnd; iz<ibnd+ntaper; iz++) {
				vx[ix*n1+iz] *= bnd.tapx[ib-iz];
				vz[ix*n1+iz+1] *= bnd.tapz[ib-iz];
			}
		}
	}
	/* right */
	if (bnd.tap[1] > 0) {
		ntaper = bnd.tap[1];
		ib = (nx+ibnd-ntaper);
		for (ix=nx+ibnd-ntaper; ix<nx+ibnd; ix++) {
#pragma ivdep
			for (iz=ibnd; iz<nz+ibnd; iz++) {
				vx[ix*n1+iz] *= bnd.tapx[ix-ib];
				vz[ix*n1+iz] *= bnd.tapz[ix-ib];
			}
		}
	}
	/* bottom */
	if (bnd.tap[2] > 0) {
		ntaper = bnd.tap[2];
		ib = (nz+ibnd-ntaper);
		for (ix=ibnd; ix<nx+ibnd; ix++) {
#pragma ivdep
			for (iz=nz+ibnd-ntaper; iz<nz+ibnd; iz++) {
				vx[ix*n1+iz]   *= bnd.tapx[iz-ib];
				vz[ix*n1+iz+1] *= bnd.tapz[iz-ib];
			}
		}
	}
	/* left */
	if (bnd.tap[3] > 0) {
		ntaper = bnd.tap[3];
		ib = (ntaper+ibnd-1);
		for (ix=ibnd; ix<ntaper+ibnd; ix++) {
#pragma ivdep
			for (iz=ibnd; iz<nz+ibnd; iz++) {
				vx[ix*n1+iz] *= bnd.tapx[ib-ix];
				vz[ix*n1+iz] *= bnd.tapz[ib-ix];
			}
		}
	}

	return 0;
}
